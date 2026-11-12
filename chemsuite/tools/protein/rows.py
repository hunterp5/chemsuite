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

"""Ligands as table rows: the columns that link a row to a protein, and the rows themselves.

A ligand row names its protein (``Protein``, the asset id), its residue (``Ligand ID``, as
:meth:`~chemsuite.protein.structure.model.Structure.residue_label` writes it) and holds the
ligand in its pose (``Ligand``, a molfile), with where its bond orders came from.
"""

from __future__ import annotations

from rdkit import Chem

from chemsuite.core.jobs import CancelToken
from chemsuite.protein.chem import ligands
from chemsuite.protein.structure.components import CcdCache
from chemsuite.protein.structure.model import Structure

PROTEIN_COLUMN = "Protein"
LIGAND_ID_COLUMN = "Ligand ID"
LIGAND_COLUMN = "Ligand"
BOND_ORDERS_COLUMN = "Bond Orders"
#: the columns that say what a row shows in the workbench: editing one follows the row again
ROW_COLUMNS = (PROTEIN_COLUMN, LIGAND_ID_COLUMN, LIGAND_COLUMN)
SOURCE_TITLES = {
    ligands.SMILES: "SMILES",
    ligands.FILE: "File",
    ligands.CCD: "CCD (downloaded)",
    ligands.TEMPLATE: "CCD template",
    ligands.GEOMETRY: "Guessed from geometry",
}


def ligand_rows(
    structure: Structure,
    protein: str,
    residues: list[int] | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
) -> dict[str, list[str]]:
    """One row per ligand (or per residue of ``residues``): the protein's asset id, the residue,
    a molfile in the ligand's pose and where its bond orders came from. A ligand RDKit cannot
    read is left out, and so are cofactors when the protein has another ligand."""
    rows: dict[str, list[str]] = {PROTEIN_COLUMN: [], LIGAND_ID_COLUMN: [], LIGAND_COLUMN: [],
                                  BOND_ORDERS_COLUMN: []}  # fmt: skip
    if residues is None:
        residues = ligands.row_residues(structure)
    for residue in residues:
        try:
            ligand = ligands.to_mol(structure, residue, ccd=ccd, cancel=cancel)
            block = Chem.MolToMolBlock(ligand.mol)
        except (ValueError, RuntimeError):
            continue
        rows[PROTEIN_COLUMN].append(protein)
        rows[LIGAND_ID_COLUMN].append(structure.residue_label(residue))
        rows[LIGAND_COLUMN].append(block)
        rows[BOND_ORDERS_COLUMN].append(SOURCE_TITLES[ligand.source])
    return rows


def guessed(rows: dict[str, list[str]]) -> list[str]:
    """The ligand IDs of ``rows`` whose bond orders are a guess."""
    pairs = zip(rows[LIGAND_ID_COLUMN], rows[BOND_ORDERS_COLUMN], strict=True)
    return [label for label, source in pairs if source == SOURCE_TITLES[ligands.GEOMETRY]]


def mol_of(given: bytes | str | None) -> Chem.Mol | None:
    """A molecule from :meth:`ToolContext.structure_inputs` (a blob, unparsed text or None)."""
    from chemsuite.core.mols import mol_from_blob, mol_from_text

    if isinstance(given, bytes):
        return mol_from_blob(given) if given else None
    return mol_from_text(given) if given else None
