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
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

from ..structure.components import CcdCache
from ..structure.model import LIGAND, WATER, Structure
from . import hydrogens, repair, solve
from .energy import Weights
from .inspect import Finding, inspect
from .templates import Library

ALL, NEAR_LIGANDS, NONE = "all", "near ligands", "none"
WATER_REACH = 5.0  # Å from a ligand atom, for NEAR_LIGANDS
AMBIGUOUS = 1.0  # a site whose next best state is closer than this is marked ⚠
SITE = "site"
NOTE = "note"


@dataclass(frozen=True)
class Options:
    ph: float = 7.0
    waters: str = ALL
    repair: bool = True
    optimise: bool = True
    weights: Weights | None = None


@dataclass(frozen=True)
class Prepared:
    structure: Structure
    states: dict[int, str]  # residue → HID, HIE, HIP, CYM…
    report: tuple[Finding, ...]
    choices: tuple[solve.Choice, ...]

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
    if options.optimise:
        protonated, choices = solve.optimise(protonated, options.ph, options.weights, cancel)
        report += _site_findings(protonated.structure, choices)
    return Prepared(protonated.structure, dict(protonated.states), tuple(report), tuple(choices))


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


def _site_findings(s: Structure, choices: list[solve.Choice]) -> list[Finding]:
    """A row for each flip, each His, and each site whose choice is close."""
    found = []
    for c in choices:
        if c.site.kind == "rotor" and c.gap >= AMBIGUOUS:
            continue
        state = c.site.states[c.state]
        if c.site.kind == "amide" and not state.flipped and c.gap >= AMBIGUOUS:
            continue
        close = c.gap < AMBIGUOUS
        what = {"amide": "Amide", "his": "His", "rotor": "Rotor"}[c.site.kind]
        text = f"{what}: {state.label}" + (f" (next best within {c.gap:.1f})" if close else "")
        found.append(Finding(SITE, s.residue_label(c.site.residue), text, warn=close,
                             residue=c.site.residue))  # fmt: skip
    return found
