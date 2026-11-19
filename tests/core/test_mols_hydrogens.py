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

"""Hydrogens as the table keeps them: on 3D structures, not on 2D ones."""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.core.mols import as_read, has_3d, with_hydrogens, without_hydrogens


def ensemble(smiles: str, count: int = 3) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, count, randomSeed=7)
    for conformer in mol.GetConformers():
        conformer.SetDoubleProp("Energy", 1.0 + conformer.GetId())
    return mol


def hydrogens(mol: Chem.Mol) -> int:
    return sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())


def test_a_3d_structure_keeps_its_hydrogens_and_a_2d_one_loses_them() -> None:
    three_d = ensemble("CCO", 1)
    assert has_3d(three_d) and as_read(three_d) is three_d
    flat = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.Compute2DCoords(flat)
    assert not has_3d(flat)
    read = as_read(flat)
    assert hydrogens(read) == 0 and Chem.MolToSmiles(read) == "CCO"
    assert as_read(None) is None


def test_hydrogens_are_placed_on_every_conformer_with_its_properties() -> None:
    full = ensemble("CC(=O)Nc1ccccc1")
    heavy = Chem.RemoveHs(full)
    back = with_hydrogens(heavy)
    assert back.GetNumAtoms() == full.GetNumAtoms() and back.GetNumConformers() == 3
    for conformer in back.GetConformers():
        positions = conformer.GetPositions()
        assert np.isfinite(positions).all() and not np.allclose(positions[-1], 0)
        assert conformer.GetDoubleProp("Energy") == 1.0 + conformer.GetId()
    one = with_hydrogens(heavy, first_conformer=True)
    assert one.GetNumConformers() == 1
    assert np.allclose(one.GetConformer().GetPositions()[: heavy.GetNumAtoms()],
                       heavy.GetConformer().GetPositions())  # fmt: skip


def test_a_complete_structure_is_left_as_it_is() -> None:
    full = ensemble("CCO")
    assert with_hydrogens(full) is full
    bare = Chem.MolFromSmiles("CCO")
    assert without_hydrogens(bare) is bare


def test_removing_hydrogens_keeps_isotopes() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("[2H]C([2H])Cl"))
    assert Chem.MolToSmiles(without_hydrogens(mol)) == "[2H]C([2H])Cl"
