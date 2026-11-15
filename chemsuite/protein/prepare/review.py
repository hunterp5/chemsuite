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
"""Everything a preparation would find and decide, laid out to be reviewed and overruled.

:func:`survey` prepares a structure and reports what it did as a list of **issues** (what the
file is missing or has twice: alternate locations, chain breaks, missing atoms, links,
disulfides, metal sites) and a list of **choices** (each flip, His tautomer, rotor, titratable
group, ligand form and water orientation, with every state it could take and what that state
would cost). Two kinds of choice are not solved for, only chosen: which **alternate location**
of a residue to keep (by occupancy), and whether each **crystal water** is kept.
Nothing is committed: the survey hands back the structure it would produce, and the chemist
decides what of it to keep.

Each line carries a stable ``id``, so a panel can remember what was ticked and what was
overruled across a re-survey, and turn the answers into
:class:`~chemsuite.protein.prepare.pipeline.Decisions` (:func:`decisions_from`). Applying them
re-runs the preparation rather than editing the proposal: a pinned site changes what its
neighbours should do, and only a fresh solve keeps the network honest.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from chemsuite.core.jobs import CancelToken

from ..structure.components import CcdCache
from ..structure.model import WATER, Structure
from . import inspect as inspection
from . import pipeline, sites
from .pipeline import Decisions, Options, Prepared

#: the issues a preparation can act on; the rest are reported and left alone (alternate
#: locations are acted on residue by residue, as :data:`ALTLOC` choices)
FIXABLE = (inspection.MISSING,)
#: the choices that are chosen, not solved for: a residue's alternate location, and whether a
#: crystal water is kept (its states :data:`KEPT` and :data:`DROPPED`)
ALTLOC, KEEP = "altloc", "keep"
KEPT, DROPPED = "Kept", "Dropped"
#: how a choice's kind reads in a review, in the order they are listed
KINDS = {
    ALTLOC: "Alternate location",
    sites.AMIDE: "Amide flip",
    sites.HIS: "His state",
    sites.ACID: "Protonation",
    sites.ROTOR: "Rotor",
    sites.LIGAND: "Ligand form",
    KEEP: "Crystal water",
    sites.WATER: "Water orientation",
}


@dataclass(frozen=True)
class Issue:
    """Something the file is missing or carries twice, and whether preparation can act on it."""

    id: str
    kind: str  # inspect.ALTLOCS, BREAK, GAP, MISSING, UNKNOWN, LINK, DISULFIDE, METAL
    label: str
    text: str
    warn: bool
    residue: int
    atoms: tuple[str, ...] = ()

    @property
    def fixable(self) -> bool:
        """Whether ticking this line changes anything (a chain break cannot be built)."""
        return self.kind in FIXABLE and self.residue >= 0


@dataclass(frozen=True)
class State:
    """One state a site could take, and what it would cost there."""

    label: str
    form: str
    energy: float  # kcal/mol-like, the other sites as they were chosen
    chosen: bool
    cost: float = 0.0  # how much worse than the site's best state (0 for the best)
    note: str = ""  # for a choice not solved for, what tells its states apart (occupancy)


@dataclass(frozen=True)
class Choice:
    """One decision the preparation made, with the states it chose between."""

    id: str
    kind: str  # sites.AMIDE, HIS, ACID, ROTOR, LIGAND, WATER
    label: str  # the residue, as the viewer names it
    residue: int
    states: tuple[State, ...]
    chosen: int
    gap: float  # to the next best state
    pka: float | None = None  # with titration
    population: float | None = None
    #: worth a look (⚠), by the same rule as the preparation report (``pipeline.warns``)
    warn: bool = False

    @property
    def what(self) -> str:
        return KINDS.get(self.kind, self.kind.title())

    @property
    def scored(self) -> bool:
        """Whether its states have energies (a site solved for), not just names."""
        return self.kind not in (ALTLOC, KEEP)

    @property
    def state(self) -> State | None:
        return self.states[self.chosen] if 0 <= self.chosen < len(self.states) else None


@dataclass(frozen=True)
class Survey:
    """What a preparation found and would do, and the structure it would produce."""

    protein: str
    issues: tuple[Issue, ...]
    choices: tuple[Choice, ...]
    prepared: Prepared
    options: Options
    #: the decisions this survey was run with (a re-survey after an override carries them)
    decisions: Decisions = field(default_factory=Decisions)

    @property
    def warnings(self) -> int:
        return sum(1 for i in self.issues if i.warn) + sum(1 for c in self.choices if c.warn)

    @property
    def net_charge(self) -> int:
        return self.prepared.net_charge

    def choice(self, id: str) -> Choice | None:
        return next((c for c in self.choices if c.id == id), None)


def survey(
    protein: str,
    s: Structure,
    options: Options | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
    decisions: Decisions | None = None,
) -> Survey:
    """Prepare *s* without committing it, and report everything found and decided."""
    options = options or Options()
    chosen = decisions or Decisions()
    prepared = pipeline.prepare(s, options, ccd, cancel, chosen, detail=True)
    issues = tuple(
        Issue(_issue_id(f), f.kind, f.label, f.text, f.warn, f.residue, f.atoms)
        for f in prepared.report
        if f.kind != pipeline.SITE and f.kind != pipeline.NOTE
    )
    altlocs = _altloc_choices(s, chosen) if options.repair else []
    choices = (*altlocs, *_choices(prepared), *_water_choices(s, prepared))
    return Survey(protein, issues, choices, prepared, options, chosen)


def _altloc_choices(s: Structure, decisions: Decisions) -> list[Choice]:
    """A line for each residue with alternate locations: its letters with their occupancies,
    the one kept chosen (the highest occupancy, unless the chemist picked another)."""
    out = []
    has = s.altloc != ""
    for r, default in sorted(s.default_altlocs.items()):
        label = s.residue_label(r)
        atoms = np.flatnonzero(has & (s.residue_of == r))
        letters = sorted({str(x) for x in s.altloc[atoms]})
        picked = decisions.altlocs.get(label, default)
        picked = picked if picked in letters else default
        states = []
        for x in letters:
            occupancy = float(s.occupancy[atoms][s.altloc[atoms] == x].mean())
            states.append(State(x, "", 0.0, x == picked, 0.0, f"occupancy {occupancy:.2f}"))
        out.append(Choice(f"{ALTLOC}:{label}", ALTLOC, label, r, tuple(states),
                          letters.index(picked), float("inf")))  # fmt: skip
    return out


def _water_choices(s: Structure, prepared: Prepared) -> list[Choice]:
    """A line for each crystal water: kept or dropped, as the preparation left it."""
    after = prepared.structure
    kept = {after.residue_label(int(r))
            for r in np.unique(after.residue_of[after.atom_kind == WATER])}  # fmt: skip
    out = []
    for r in np.unique(s.residue_of[s.atom_kind == WATER]):
        label = s.residue_label(int(r))
        on = label in kept
        states = (State(KEPT, "", 0.0, on), State(DROPPED, "", 0.0, not on))
        out.append(Choice(f"{KEEP}:{label}", KEEP, label, int(r), states, 0 if on else 1,
                          float("inf")))  # fmt: skip
    return out


def _issue_id(finding: inspection.Finding) -> str:
    """A name for the issue that still means it when the document reopens, or after an edit
    elsewhere in the structure: its kind, its residue's label and its atoms."""
    atoms = ",".join(finding.atoms)
    return f"{finding.kind}:{finding.label}:{atoms}"


def _choices(prepared: Prepared) -> tuple[Choice, ...]:
    """Every decision, with each state's energy where the energies were kept."""
    energies = prepared.energies
    keys = prepared.keys
    out = []
    for i, choice in enumerate(prepared.choices):
        site = choice.site
        found = energies[i] if i < len(energies) else None
        best = float(np.min(found)) if found is not None and len(found) else 0.0
        states = tuple(
            State(
                state.label,
                state.form,
                float(found[k]) if found is not None and k < len(found) else 0.0,
                k == choice.state,
                (float(found[k]) - best) if found is not None and k < len(found) else 0.0,
            )
            for k, state in enumerate(site.states)
        )
        out.append(
            Choice(
                keys[i] if i < len(keys) else pipeline.site_key(prepared.structure, site),
                site.kind,
                prepared.structure.residue_label(site.residue),
                site.residue,
                states,
                choice.state,
                choice.gap,
                prepared.pka.get(site.residue),
                prepared.population.get(site.residue),
                pipeline.warns(choice, pipeline.site_fraction(choice, prepared.population)),
            )
        )
    return tuple(out)


def decisions_from(
    found: Survey,
    *,
    skip_issues: Mapping[str, bool] | None = None,
    overrides: Mapping[str, str] | None = None,
) -> Decisions:
    """The chemist's answers as :class:`Decisions`.

    *found* must be the survey run **without** decisions (the proposals): the answers are read
    against it, and nothing it was itself run with is carried over, so a re-survey that held a
    choice never turns that choice into a silent default. *skip_issues* marks issue ids to
    leave alone (untick a repair); *overrides* holds choice ids against the state label each
    must take (a site's state, a residue's alternate location, a water kept or dropped). What
    is not named is left to the preparation.
    """
    skipped = {id for id, skip in (skip_issues or {}).items() if skip}
    repair: set[str] | None = None
    if any(i.kind == inspection.MISSING and i.id in skipped for i in found.issues):
        repair = {
            i.label
            for i in found.issues
            if i.kind == inspection.MISSING and i.residue >= 0 and i.id not in skipped
        }
    states: dict[str, str] = {}
    altlocs: dict[str, str] = {}
    keep: set[str] = set()
    drop: set[str] = set()
    for id, label in (overrides or {}).items():
        choice = found.choice(id)
        if choice is None or all(state.label != label for state in choice.states):
            continue
        if choice.kind == ALTLOC:
            altlocs[choice.label] = label
        elif choice.kind == KEEP:
            (keep if label == KEPT else drop).add(choice.label)
        else:
            states[choice.id] = label
    return Decisions(
        repair=frozenset(repair) if repair is not None else None,
        keep_waters=frozenset(keep),
        drop_waters=frozenset(drop),
        altlocs=altlocs,
        states=states,
    )
