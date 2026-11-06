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

"""Implicit hydrogens and valence problems, by RDKit's model (SKETCHER_PLAN §4.1).

Chosen in P0 (docs/sketcher/valence-review.md): RDKit's model got 59 of 61 reference species
right, Ketcher's CTfile model 57. A bad atom is marked with the reason, not refused, as Ketcher
does with ``badConn``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem, rdBase

from chemsuite.sketch.convert import SketchError, atom_order, to_rdkit
from chemsuite.sketch.model import Sketch


@dataclass(frozen=True)
class Report:
    hydrogens: dict[int, int] = field(default_factory=dict)  # atom id -> implicit H
    problems: dict[int, str] = field(default_factory=dict)  # atom id -> what is wrong

    @property
    def valid(self) -> bool:
        return not self.problems


def check(sketch: Sketch) -> Report:
    """Hydrogens to show on each atom, and the atoms RDKit rejects with the reason."""
    if sketch.is_empty():
        return Report()
    ids = atom_order(sketch)
    try:
        mol = to_rdkit(sketch)
    except SketchError:
        return _diagnose(sketch, ids)
    return Report({ids[a.GetIdx()]: a.GetTotalNumHs() for a in mol.GetAtoms()})


def _diagnose(sketch: Sketch, ids: list[int]) -> Report:
    raw = to_rdkit(sketch, sanitize=False)
    with rdBase.BlockLogs():
        found = Chem.DetectChemistryProblems(raw)
    problems: dict[int, str] = {}
    for problem in found:
        if problem.GetType() == "AtomValenceException":
            index = problem.GetAtomIdx()
            problems[ids[index]] = _too_many_bonds(raw.GetAtomWithIdx(index))
        elif problem.GetType() == "KekulizeException":
            for index in problem.GetAtomIndices():
                problems[ids[index]] = (
                    "These aromatic bonds cannot be made alternating single and double bonds."
                )
        elif hasattr(problem, "GetAtomIdx"):
            problems[ids[problem.GetAtomIdx()]] = problem.Message()
    if not problems:  # rejected for a reason not tied to an atom: blame every atom
        problems = dict.fromkeys(ids, "RDKit cannot read this drawing.")
    raw.UpdatePropertyCache(strict=False)
    hydrogens = {ids[a.GetIdx()]: a.GetTotalNumHs() for a in raw.GetAtoms()}
    for aid in problems:
        hydrogens[aid] = 0
    return Report(hydrogens, problems)


def _too_many_bonds(atom: Chem.Atom) -> str:
    bonds = sum(b.GetValenceContrib(atom) for b in atom.GetBonds())
    label = atom.GetSymbol()
    if atom.GetFormalCharge():
        sign = "+" if atom.GetFormalCharge() > 0 else "−"
        label += f"{abs(atom.GetFormalCharge()) if abs(atom.GetFormalCharge()) > 1 else ''}{sign}"
    return f"Too many bonds: {label} with {bonds:g} bond orders is not an allowed valence."
