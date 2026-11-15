# This file is part of chemsuite.
# Copyright (C) 2026 Hunter Picard
#
# chemsuite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# chemsuite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with chemsuite.  If not, see <https://www.gnu.org/licenses/>.

"""Prepare a structure in one call: inspect, repair, waters, hydrogens, flips and states.

The result is a new structure (old atoms keep their uids, so the workbench can apply it as one
undoable edit) and a report: what inspection found, what each site chose and how sure it is
(:data:`AMBIGUOUS`), and what could not be done.

:class:`Options` is how the preparation is run; :class:`Decisions` is the chemist overruling
it case by case — which residues to repair, which waters to keep, which location of an
alternate, and which state a given site must take. A pinned site is held to that state and
**the rest of the network is solved around it** (``solve.solve``'s ``allowed``), so the
coupling between sites is kept rather than one answer being edited in afterwards. Reviewing
and overruling is :mod:`.review`.

Every titrating group (:attr:`Options.titrate`) takes its protonation state at the pH from
:mod:`.populations`: the structure carries each group's most populated form (the best state
within it), and the result each group's pKa and protonated share. A group between
:data:`MIXED` and 1 − :data:`MIXED` is marked ⚠: both forms are present at this pH.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

from ..structure.components import CcdCache
from ..structure.model import LIGAND, WATER, Structure
from . import energy, hydrogens, ligand_sites, populations, repair, sites, solve, waters
from .energy import Weights
from .inspect import Finding, inspect
from .templates import Library

ALL, NEAR_LIGANDS, NONE = "all", "near ligands", "none"
WATER_REACH = 5.0  # Å from a ligand atom, for NEAR_LIGANDS
POCKET_REACH = 6.0  # Å from a ligand atom: the residues that titrate with titrate=POCKET
AMBIGUOUS = 1.0  # a site whose next best state is closer than this is marked ⚠
#: a rotor whose chosen orientation scores better than this makes a hydrogen bond; one that
#: does not is free to spin, so its near ties say nothing (as a water's do not) and it is
#: counted in one line instead of a row of its own
BONDED = -0.1
#: a rotor whose next best orientation is this close makes the same bonds (an ammonium's three
#: hydrogens, a hydroxyl stepping past its acceptor): either answer is the same answer, so it is
#: not reported. Between this and :data:`AMBIGUOUS` the orientations differ and the user's eye
#: is worth having; an amide or His tie is reported however small, the flip being a real question
DEGENERATE = 0.05
MIXED = 0.1  # a titrating group with more than this share of each form is marked ⚠
#: the forms a report names (the others are the textbook ones at pH 7)
UNCOMMON = {"ASH": "Asp", "GLH": "Glu", "LYN": "Lys", "TYM": "Tyr", "CYM": "Cys",
            "NH2": "N-terminus", "HIP": "His"}  # fmt: skip
SITE = "site"
NOTE = "note"


@dataclass(frozen=True)
class Decisions:
    """What the chemist chose, against what the preparation would have done on its own.

    Everything is optional: the default decides nothing, so ``prepare`` runs as it always has.
    """

    #: residues whose missing atoms are built; None for every one that can be
    repair: frozenset[int] | None = None
    #: waters kept whatever :attr:`Options.waters` says, and waters dropped whatever it says
    keep_waters: frozenset[int] = frozenset()
    drop_waters: frozenset[int] = frozenset()
    #: residue → which alternate location to keep (the highest occupancy by default)
    altlocs: Mapping[int, str] = field(default_factory=dict)
    #: a site's key (:func:`site_key`) → the state that site must take, by its label
    states: Mapping[str, str] = field(default_factory=dict)

    def pinned(self, key: str, site: sites.Site) -> int | None:
        """Which state the site *key* is held to, if any (None: it chooses for itself)."""
        label = self.states.get(key)
        if label is None:
            return None
        for i, state in enumerate(site.states):
            if state.label == label:
                return i
        return None


def site_key(s: Structure, site: sites.Site) -> str:
    """A name for *site* that survives a second preparation of the same structure: its kind,
    its residue, and the atoms it hangs from and moves (a phosphate has three rotors, one
    residue)."""
    anchor = str(s.name[site.anchors[0]]) if len(site.anchors) else ""
    moving = str(s.name[site.moving[0]]) if len(site.moving) else ""
    return f"{site.kind}:{site.residue}:{anchor}/{moving}"


@dataclass(frozen=True)
class Options:
    ph: float = 7.0
    waters: str = ALL
    repair: bool = True
    optimise: bool = True
    weights: Weights | None = None
    titrate: str = sites.OFF  # OFF, POCKET or ALL (see :func:`.sites.find_sites`)
    ligand_states: bool = True  # each ligand's protomers and tautomers (:mod:`.ligand_sites`)
    orient_waters: bool = True  # hydrogens on the waters near ligands (:mod:`.waters`)


@dataclass(frozen=True)
class Prepared:
    structure: Structure
    states: dict[int, str]  # residue → HID, HIE, HIP, CYM…
    report: tuple[Finding, ...]
    choices: tuple[solve.Choice, ...]
    pka: dict[int, float] = field(default_factory=dict)  # residue → pKa (with titration)
    #: residue → protonated share; a ligand's: the share of the form it takes
    population: dict[int, float] = field(default_factory=dict)
    #: with ``detail``: every state's energy at each site of :attr:`choices`, in its order, so
    #: a review can show what the chosen state was chosen against (:mod:`.review`)
    energies: tuple[np.ndarray, ...] = ()
    #: each choice's :func:`site_key`, in the order of :attr:`choices`
    keys: tuple[str, ...] = ()

    @property
    def net_charge(self) -> int:
        return int(self.structure.charge.sum())


def prepare(
    s: Structure,
    options: Options | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
    decisions: Decisions | None = None,
    detail: bool = False,
) -> Prepared:
    """*s* prepared. With *decisions*, the chemist's choices overrule what it would do on its
    own, case by case (:class:`Decisions`); with *detail*, every state's energy is kept as
    well, for a review of the choices (:mod:`.review`)."""
    options = options or Options()
    chosen = decisions or Decisions()
    library = Library(ccd)
    report = inspect(s, library, cancel)
    if options.repair:
        s = repair.build_missing(repair.one_location(s, dict(chosen.altlocs) or None),
                                 library, cancel, only=chosen.repair)  # fmt: skip
    check(cancel)
    s = keep_waters(s, options.waters, chosen)
    protonated = hydrogens.add_hydrogens(s, library, ccd, cancel)
    kept = frozenset(n.residue for n in protonated.notes if n.text == hydrogens.KEPT)
    report += [Finding(NOTE, s.residue_label(n.residue), n.text, warn=True, residue=n.residue)
               for n in protonated.notes if n.text != hydrogens.KEPT]  # fmt: skip
    choices: list[solve.Choice] = []
    energies: tuple[np.ndarray, ...] = ()
    keys: tuple[str, ...] = ()
    pka: dict[int, float] = {}
    share: dict[int, float] = {}
    if options.optimise:
        weights = options.weights or Weights()
        found = _sites(protonated, options, ccd, kept, report, cancel)
        check(cancel)
        table = energy.tables(protonated.structure, found, options.ph)
        check(cancel)
        keys = tuple(site_key(protonated.structure, site) for site in found)
        held = _held(found, keys, chosen)
        picked = solve.solve(table, weights, allowed=held or None)
        if any(site.titrating for site in found):
            done = populations.titration(table, weights, cancel=cancel)
            at = done.at(options.ph)
            states_at = [np.array([np.interp(options.ph, done.grid, column) for column in pop.T])
                         for pop in done.populations]  # fmt: skip
            picked = _majority(table, weights, picked, at, states_at)
            for i, (site, k, fraction) in enumerate(zip(found, done.pka, at, strict=True)):
                if site.kind == sites.LIGAND:
                    form = site.states[picked[i]].form
                    share[site.residue] = float(sum(
                        w for st, w in zip(site.states, states_at[i], strict=True)
                        if st.form == form))  # fmt: skip
                elif site.titrating:
                    pka[site.residue], share[site.residue] = k, fraction
        choices = solve.choices(table, weights, picked)
        if detail:
            energies = tuple(solve.energies_of(table, weights, picked))
        protonated = solve.apply(protonated, found, picked)
        report += _site_findings(protonated.structure, choices, share)
        report += _rotor_finding(choices)
        report += _water_finding(choices)
    return Prepared(protonated.structure, dict(protonated.states), tuple(report), tuple(choices),
                    pka, share, energies, keys)  # fmt: skip


def _sites(
    protonated: hydrogens.Protonated,
    options: Options,
    ccd: CcdCache | None,
    kept: frozenset[int],
    report: list[Finding],
    cancel: CancelToken | None,
) -> list[sites.Site]:
    """The protein's sites and (with :attr:`Options.ligand_states`) the ligands'."""
    s = protonated.structure
    found: list[sites.Site] = []
    if options.ligand_states:
        found, notes = ligand_sites.find(protonated, options.ph, ccd, kept, cancel)
        report += [Finding(NOTE, s.residue_label(n.residue), n.text, warn=True,
                           residue=n.residue) for n in notes]  # fmt: skip
    taken = np.zeros(s.n_atoms, dtype=bool)
    for site in found:
        taken[site.moving] = True
    if options.orient_waters:
        found += waters.find(protonated)
    pocket = _pocket(s) if options.titrate == sites.POCKET else None
    return sites.find_sites(protonated, options.titrate, pocket=pocket, taken=taken) + found


def _pocket(s: Structure) -> np.ndarray:
    """The residues with an atom within :data:`POCKET_REACH` of a ligand's (the ones that
    titrate with :attr:`Options.titrate` ``POCKET``)."""
    near = np.zeros(s.n_residues, dtype=bool)
    ligand = s.atom_kind == LIGAND
    if not ligand.any():
        return near
    x = s.xyz()
    distance = cKDTree(x[ligand]).query(x, distance_upper_bound=POCKET_REACH)[0]
    near[s.residue_of[np.isfinite(distance) & ~ligand]] = True
    return near


def _majority(
    table: energy.Tables,
    weights: Weights,
    picked: list[int],
    share: list[float],
    populations_at: list[np.ndarray],
) -> list[int]:
    """*picked*, with every titrating site in its most populated form (the lowest energy can
    differ, as a form with many rotor states gains from their number). A protein group's form
    is its proton count; a ligand's, its protomer or tautomer."""
    allowed = {}
    for i, (n, fraction) in enumerate(zip(table.protons, share, strict=True)):
        site = table.sites[i]
        if site.kind == sites.LIGAND:
            weight: dict[str, float] = {}
            for st, w in zip(site.states, populations_at[i], strict=True):
                weight[st.form] = weight.get(st.form, 0.0) + float(w)
            best = max(weight, key=lambda form: weight[form])
            allowed[i] = np.array([st.form == best for st in site.states])
        elif len(n) and n.max() > n.min():
            allowed[i] = n == (n.max() if fraction >= 0.5 else n.min())
    if all(mask[picked[i]] for i, mask in allowed.items()):
        return picked
    return solve.solve(table, weights, allowed=allowed)


def _held(
    found: list[sites.Site], keys: tuple[str, ...], chosen: Decisions
) -> dict[int, np.ndarray]:
    """Each pinned site held to its one state (a mask for ``solve.solve``), so the others are
    solved around it."""
    held: dict[int, np.ndarray] = {}
    for i, (site, key) in enumerate(zip(found, keys, strict=True)):
        state = chosen.pinned(key, site)
        if state is not None:
            mask = np.zeros(len(site.states), dtype=bool)
            mask[state] = True
            held[i] = mask
    return held


def keep_waters(s: Structure, which: str, chosen: Decisions | None = None) -> Structure:
    """*s* with all its waters, none, or those within :data:`WATER_REACH` of a ligand; the
    chemist's own keeps and drops (:class:`Decisions`) overrule that either way."""
    chosen = chosen or Decisions()
    water = s.atom_kind == WATER
    if which == ALL and not chosen.drop_waters:
        return s
    if not water.any():
        return s
    drop = water.copy()
    if which == ALL:
        drop[:] = False
    elif which == NEAR_LIGANDS:
        ligand = np.flatnonzero(s.atom_kind == LIGAND)
        if len(ligand):
            x = s.xyz()
            near = cKDTree(x[ligand]).query(x[water], distance_upper_bound=WATER_REACH)[0]
            keep = np.flatnonzero(water)[np.isfinite(near)]
            residues = np.unique(s.residue_of[keep])
            drop &= ~np.isin(s.residue_of, residues)
    if chosen.keep_waters:
        drop &= ~np.isin(s.residue_of, list(chosen.keep_waters))
    if chosen.drop_waters:
        drop |= water & np.isin(s.residue_of, list(chosen.drop_waters))
    return s.subset(~drop)


def warns(c: solve.Choice, fraction: float | None = None) -> bool:
    """Whether a choice is worth the chemist's eye (⚠), the one rule the report and a review
    (:mod:`.review`) both use. *fraction* is a titrating site's protonated share (a ligand's:
    the share of the form it takes).

    - A water never: its free hydrogen nearly always has a second orientation as good.
    - A titrating site, or a ligand, when both its forms are present at this pH.
    - A rotor not when it is free to spin (no hydrogen bond) or its next best makes the same
      bonds (:data:`BONDED`, :data:`DEGENERATE`): either answer is the same answer.
    - Anything else when the next best state is within :data:`AMBIGUOUS`.
    """
    kind = c.site.kind
    if kind == sites.WATER:
        return False
    if kind == sites.LIGAND:
        return fraction is not None and fraction < 1 - MIXED
    if fraction is not None:
        return MIXED < fraction < 1 - MIXED
    if kind == sites.ROTOR and (c.energy > BONDED or c.gap <= DEGENERATE):
        return False
    return c.gap < AMBIGUOUS


def site_fraction(c: solve.Choice, share: Mapping[int, float]) -> float | None:
    """A site's protonated share, if it titrated (its rotors do not)."""
    return share.get(c.site.residue) if c.site.titrating else None


def _site_findings(
    s: Structure, choices: list[solve.Choice], share: dict[int, float] | None = None
) -> list[Finding]:
    """A row for each flip, each His, each uncommon protonation state, and each site whose
    choice is close (or, with titration, whose forms are mixed)."""
    share = share or {}
    found = []
    for c in choices:
        state = c.site.states[c.state]
        fraction = site_fraction(c, share)
        if c.site.kind == sites.LIGAND:
            found += _ligand_finding(s, c, fraction)
            continue
        if c.site.kind == sites.WATER:
            continue  # one row for them all (_water_finding)
        if fraction is not None and c.site.kind != sites.HIS:  # titration decides, not the gap
            mixed = warns(c, fraction)
            if state.form not in UNCOMMON and not mixed:
                continue
            what = UNCOMMON.get(state.form) or state.form.title()
            text = f"{what}: {state.form}, {fraction:.0%} protonated"
            found.append(Finding(SITE, s.residue_label(c.site.residue), text, warn=mixed,
                                 residue=c.site.residue))  # fmt: skip
            continue
        if c.site.kind == sites.ROTOR and (c.energy > BONDED or c.gap <= DEGENERATE):
            continue  # free to spin, or the same bonds either way (_rotor_finding counts them)
        if c.site.kind in (sites.ROTOR, sites.ACID) and c.gap >= AMBIGUOUS:
            continue
        if c.site.kind == sites.AMIDE and not state.flipped and c.gap >= AMBIGUOUS:
            continue
        close = warns(c, fraction)
        what = {sites.AMIDE: "Amide", sites.HIS: "His", sites.ROTOR: "Rotor",
                sites.ACID: "Acid"}[c.site.kind]  # fmt: skip
        text = f"{what}: {state.label}" + (f" (next best within {c.gap:.1f})" if close else "")
        if fraction is not None:
            text = f"{what}: {state.label}, {fraction:.0%} protonated"
        found.append(Finding(SITE, s.residue_label(c.site.residue), text, warn=close,
                             residue=c.site.residue))  # fmt: skip
    return found


def _ligand_finding(s: Structure, c: solve.Choice, fraction: float | None) -> list[Finding]:
    """A row for a ligand not as drawn, or whose forms are mixed."""
    state = c.site.states[c.state]
    mixed = warns(c, fraction)
    if state.form == ligand_sites.AS_DRAWN and not mixed:
        return []
    text = f"Ligand: {state.form}"
    if state.form != ligand_sites.AS_DRAWN:
        text += f" (charge on these atoms {int(sum(state.formal)):+d})"
    if fraction is not None:
        text += f", {fraction:.0%} of it"
    return [Finding(SITE, s.residue_label(c.site.residue), text, warn=mixed,
                    residue=c.site.residue)]  # fmt: skip


def _rotor_finding(choices: list[solve.Choice]) -> list[Finding]:
    """One row for the hydroxyl and amine hydrogens with nothing to bind to. No ⚠: they are
    free to spin, so the orientation they were given means nothing (as a water's does not)."""
    free = sum(1 for c in choices
               if c.site.kind == sites.ROTOR and c.energy > BONDED)  # fmt: skip
    if not free:
        return []
    text = f"Hydroxyl and amine hydrogens with no hydrogen bond: {free}"
    return [Finding(SITE, "Rotors", text)]


def _water_finding(choices: list[solve.Choice]) -> list[Finding]:
    """One row for the waters near ligands given hydrogens. No ⚠: a water's free hydrogen
    nearly always has a second orientation as good, so a near tie says nothing."""
    oriented = sum(1 for c in choices if c.site.kind == sites.WATER)
    if not oriented:
        return []
    return [Finding(SITE, "Waters", f"Waters near ligands oriented: {oriented}")]
