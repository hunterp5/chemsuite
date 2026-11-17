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
"""The document's proteins' ligands, in their bound poses, for tools that take a reference
(Chemistry ▸ 3D ▸ Superpose)."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from chemsuite.app.document.tools import ToolContext
from chemsuite.macromolecule.chem import ligands
from chemsuite.macromolecule.structure import store

from .tool import KIND, ccd


@dataclass(frozen=True)
class LigandRef:
    protein: str  # the protein asset's id
    protein_name: str
    label: str  # the residue, as ``A/ATP 355``

    @property
    def title(self) -> str:
        return f"{self.protein_name} ▸ {self.label}"


def ligand_refs(ctx: ToolContext) -> list[LigandRef]:
    """Every ligand of every protein in the document, in document order."""
    found = []
    for asset in ctx.table.assets.values():
        if asset.kind != KIND:
            continue
        structure = store.as_structure(asset.payload)
        for residue in ligands.ligand_residues(structure).tolist():
            found.append(LigandRef(asset.id, asset.name, structure.residue_label(residue)))
    return found


def ligand_mol(ctx: ToolContext, ref: LigandRef) -> Chem.Mol | None:
    """The ligand's heavy atoms in their pose (bond orders from the file, the dictionary
    already on disk, or a template); None if the protein or the ligand has gone."""
    asset = ctx.table.assets.get(ref.protein)
    if asset is None or asset.kind != KIND:
        return None
    structure = store.as_structure(asset.payload)
    residue = ligands.residue_named(structure, ref.label)
    if residue < 0:
        return None
    return ligands.to_mol(structure, residue, ccd=ccd(ctx, offline=True)).mol
