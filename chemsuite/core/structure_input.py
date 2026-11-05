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

"""Reading one structure typed or pasted by the user, with a reason when it is not valid.

Bulk loading uses :func:`chemsuite.core.mols.parse_structure`; this is the slower, chattier
path behind the structure editor: it says *why* text is not a structure and offers the
formula and weight for a preview.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Literal

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors

from chemsuite.core.mols import looks_like_molblock, mol_to_blob, smiles_of
from chemsuite.core.table.columns import StructureColumn

InputKind = Literal["blank", "smiles", "molfile", "inchi"]


@dataclass(frozen=True)
class StructureInput:
    text: str
    kind: InputKind
    mol: Chem.Mol | None
    error: str = ""

    @property
    def is_blank(self) -> bool:
        return self.kind == "blank"

    @property
    def is_valid(self) -> bool:
        """True for a structure or for blank text (which clears a cell)."""
        return self.is_blank or self.mol is not None

    @cached_property
    def display(self) -> str:
        """The cell text: SMILES or InChI as typed, canonical SMILES for a molfile."""
        if self.is_blank:
            return ""
        if self.kind == "molfile":
            return smiles_of(self.mol) if self.mol is not None else self.text
        return self.text.strip()

    @cached_property
    def formula(self) -> str:
        return rdMolDescriptors.CalcMolFormula(self.mol) if self.mol is not None else ""

    @cached_property
    def mol_weight(self) -> float | None:
        return Descriptors.MolWt(self.mol) if self.mol is not None else None

    def canonical_smiles(self) -> str:
        return smiles_of(self.mol) if self.mol is not None else ""

    def as_column(self) -> StructureColumn:
        """A one-cell structure column, already parsed (a molfile keeps its coordinates)."""
        if self.is_blank:
            return StructureColumn.blanks(1)
        blob = mol_to_blob(self.mol) if self.mol is not None else b""
        return StructureColumn.from_parsed([blob], [self.display])


def read_structure(text: str) -> StructureInput:
    """Parse SMILES, a molfile block or an InChI, explaining failures."""
    molfile = looks_like_molblock(text)
    clean = text.rstrip() if molfile else text.strip()
    if not clean:
        return StructureInput(text, "blank", None)
    kind: InputKind = "molfile" if molfile else "inchi" if clean.startswith("InChI=") else "smiles"
    with rdBase.BlockLogs():
        if kind == "inchi":
            mol = Chem.MolFromInchi(clean)
            return StructureInput(clean, kind, mol, "" if mol else "Not a valid InChI.")
        if kind == "molfile":
            mol = Chem.MolFromMolBlock(clean)
            raw = None if mol else Chem.MolFromMolBlock(clean, sanitize=False)
        else:
            mol = Chem.MolFromSmiles(clean)
            raw = None if mol else Chem.MolFromSmiles(clean, sanitize=False)
    if mol is not None:
        return StructureInput(clean, kind, mol)
    return StructureInput(clean, kind, None, _explain(kind, raw))


def _explain(kind: InputKind, raw: Chem.Mol | None) -> str:
    if raw is None:
        return (
            "Not a valid molfile." if kind == "molfile" else "Not valid SMILES (check the syntax)."
        )
    with rdBase.BlockLogs():
        problems = Chem.DetectChemistryProblems(raw)
    if not problems:
        return "The structure could not be read."
    problem = problems[0]
    kind_name = problem.GetType()
    if kind_name == "AtomValenceException":
        atom = raw.GetAtomWithIdx(problem.GetAtomIdx())
        return f"Atom {atom.GetIdx() + 1} ({atom.GetSymbol()}) has too many bonds."
    if kind_name == "KekulizeException":
        atoms = ", ".join(str(i + 1) for i in problem.GetAtomIndices())
        return (
            f"Aromatic ring cannot be drawn with alternating bonds (atoms {atoms}). "
            "Check lowercase atoms and [nH]."
        )
    message = problem.Message()
    return message[:1].upper() + message[1:] + ("" if message.endswith(".") else ".")


def molblock_of(mol: Chem.Mol) -> str:
    """A V2000 molfile (2D coordinates computed when the molecule has none)."""
    if mol.GetNumConformers() == 0:
        mol = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)
