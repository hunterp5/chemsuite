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
"""A prepared receptor and its bound ligand, made once for the docking tests."""

from __future__ import annotations

from dataclasses import dataclass

import pytest
from rdkit import Chem

from chemsuite.macromolecule.chem import ligands
from chemsuite.macromolecule.dock import box as box_tools
from chemsuite.macromolecule.dock import pdbqt
from chemsuite.macromolecule.dock.box import Box
from chemsuite.macromolecule.prepare import pipeline
from chemsuite.macromolecule.structure.model import Structure
from tests.macromolecule.conftest import load


@dataclass(frozen=True)
class Pocket:
    """The EGFR pocket of 5GNK, prepared, with its covalent inhibitor as a molecule."""

    structure: Structure
    receptor: pdbqt.Receptor
    mol: Chem.Mol  # 80U where the crystal put it
    residue: int
    box: Box


@pytest.fixture(scope="session")
def pocket() -> Pocket:
    prepared = pipeline.prepare(load("5gnk_pocket.cif.gz")).structure
    residue = next(
        int(r) for r in ligands.ligand_residues(prepared) if str(prepared.res_name[r]) == "80U"
    )
    return Pocket(
        prepared,
        pdbqt.receptor(prepared),
        ligands.to_mol(prepared, residue).mol,
        residue,
        box_tools.from_residue(prepared, residue),
    )
