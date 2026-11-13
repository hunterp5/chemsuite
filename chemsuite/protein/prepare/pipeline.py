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

Every titrating group (:attr:`Options.titrate`) takes its protonation state at the pH from
:mod:`.populations`: the structure carries each group's most populated form (the best state
within it), and the result each group's pKa and protonated share. A group between
:data:`MIXED` and 1 − :data:`MIXED` is marked ⚠: both forms are present at this pH.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

from ..structure.components import CcdCache
from ..structure.model import LIGAND, WATER, Structure
from . import energy, hydrogens, populations, repair, sites, solve
from .energy import Weights
from .inspect import Finding, inspect
from .templates import Library

ALL, NEAR_LIGANDS, NONE = "all", "near ligands", "none"
WATER_REACH = 5.0  # Å from a ligand atom, for NEAR_LIGANDS
AMBIGUOUS = 1.0  # a site whose next best state is closer than this is marked ⚠
MIXED = 0.1  # a titrating group with more than this share of each form is marked ⚠
#: the forms a report names (the others are the textbook ones at pH 7)
UNCOMMON = {"ASH": "Asp", "GLH": "Glu", "LYN": "Lys", "TYM": "Tyr", "CYM": "Cys",
            "NH2": "N-terminus", "HIP": "His"}  # fmt: skip
SITE = "site"
NOTE = "note"


@dataclass(frozen=True)
class Options:
    ph: float = 7.0
    waters: str = ALL
    repair: bool = True
    optimise: bool = True
    weights: Weights | None = None
    titrate: str = sites.OFF  # OFF, POCKET or ALL (see :func:`.sites.find_sites`)


@dataclass(frozen=True)
class Prepared:
    structure: Structure
    states: dict[int, str]  # residue → HID, HIE, HIP, CYM…
    report: tuple[Finding, ...]
    choices: tuple[solve.Choice, ...]
    pka: dict[int, float] = field(default_factory=dict)  # residue → pKa (with titration)
    population: dict[int, float] = field(default_factory=dict)  # residue → protonated share

    @property
    def net_charge(self) -> int:
        return int(self.structure.charge.sum())


def prepare(
    s: Structure,
    options: Options | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
) -> Prepared:
    options = options or Options()
    library = Library(ccd)
    report = inspect(s, library, cancel)
    if options.repair:
        s = repair.build_missing(repair.one_location(s), library, cancel)
    check(cancel)
    s = keep_waters(s, options.waters)
    protonated = hydrogens.add_hydrogens(s, library, ccd, cancel)
    report += [Finding(NOTE, s.residue_label(n.residue), n.text, warn=True, residue=n.residue)
               for n in protonated.notes if n.text != "Kept the hydrogens it has"]  # fmt: skip
    choices: list[solve.Choice] = []
    pka: dict[int, float] = {}
    share: dict[int, float] = {}
    if options.optimise:
        weights = options.weights or Weights()
        found = sites.find_sites(protonated, options.titrate)
        check(cancel)
        table = energy.tables(protonated.structure, found, options.ph)
        check(cancel)
        picked = solve.solve(table, weights)
        if any(site.titrating for site in found):
            done = populations.titration(table, weights, cancel=cancel)
            at = done.at(options.ph)
            picked = _majority(table, weights, picked, at)
            for site, k, fraction in zip(found, done.pka, at, strict=True):
                if site.titrating:
                    pka[site.residue], share[site.residue] = k, fraction
        choices = solve.choices(table, weights, picked)
        protonated = solve.apply(protonated, found, picked)
        report += _site_findings(protonated.structure, choices, share)
    return Prepared(protonated.structure, dict(protonated.states), tuple(report), tuple(choices),
                    pka, share)  # fmt: skip


def _majority(
    table: energy.Tables, weights: Weights, picked: list[int], share: list[float]
) -> list[int]:
    """*picked*, with every titrating site in its most populated form (the lowest energy can
    differ, as a form with many rotor states gains from their number)."""
    allowed = {}
    for i, (n, fraction) in enumerate(zip(table.protons, share, strict=True)):
        if len(n) and n.max() > n.min():
            allowed[i] = n == (n.max() if fraction >= 0.5 else n.min())
    if all(mask[picked[i]] for i, mask in allowed.items()):
        return picked
    return solve.solve(table, weights, allowed=allowed)


def keep_waters(s: Structure, which: str) -> Structure:
    """*s* with all its waters, none, or those within :data:`WATER_REACH` of a ligand."""
    water = s.atom_kind == WATER
    if which == ALL or not water.any():
        return s
    drop = water.copy()
    if which == NEAR_LIGANDS:
        ligand = np.flatnonzero(s.atom_kind == LIGAND)
        if len(ligand):
            x = s.xyz()
            near = cKDTree(x[ligand]).query(x[water], distance_upper_bound=WATER_REACH)[0]
            keep = np.flatnonzero(water)[np.isfinite(near)]
            residues = np.unique(s.residue_of[keep])
            drop &= ~np.isin(s.residue_of, residues)
    return s.subset(~drop)


def _site_findings(
    s: Structure, choices: list[solve.Choice], share: dict[int, float] | None = None
) -> list[Finding]:
    """A row for each flip, each His, each uncommon protonation state, and each site whose
    choice is close (or, with titration, whose forms are mixed)."""
    share = share or {}
    found = []
    for c in choices:
        state = c.site.states[c.state]
        fraction = share.get(c.site.residue)
        if fraction is not None and c.site.kind != sites.HIS:  # titration decides, not the gap
            mixed = MIXED < fraction < 1 - MIXED
            if state.form not in UNCOMMON and not mixed:
                continue
            what = UNCOMMON.get(state.form) or state.form.title()
            text = f"{what}: {state.form}, {fraction:.0%} protonated"
            found.append(Finding(SITE, s.residue_label(c.site.residue), text, warn=mixed,
                                 residue=c.site.residue))  # fmt: skip
            continue
        if c.site.kind in (sites.ROTOR, sites.ACID) and c.gap >= AMBIGUOUS:
            continue
        if c.site.kind == sites.AMIDE and not state.flipped and c.gap >= AMBIGUOUS:
            continue
        close = c.gap < AMBIGUOUS
        what = {sites.AMIDE: "Amide", sites.HIS: "His", sites.ROTOR: "Rotor",
                sites.ACID: "Acid"}[c.site.kind]  # fmt: skip
        text = f"{what}: {state.label}" + (f" (next best within {c.gap:.1f})" if close else "")
        if fraction is not None:
            close = MIXED < fraction < 1 - MIXED
            text = f"{what}: {state.label}, {fraction:.0%} protonated"
        found.append(Finding(SITE, s.residue_label(c.site.residue), text, warn=close,
                             residue=c.site.residue))  # fmt: skip
    return found
