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

"""The salts-and-solvents table. Each entry was checked against PubChem by InChIKey
(scripts/review/salts_review.py, docs/sketcher/salts-review.md); these tests keep the table
sound offline and pin the corrections made to Ketcher's records."""

from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import rdMolTransforms as transforms

from chemsuite.sketch import convert, salts
from chemsuite.sketch.edits import History

BY_NAME = {e.name: e for e in salts.TABLE}


def test_every_entry_is_a_valid_distinct_structure():
    assert len(salts.TABLE) == 135  # one per Ketcher record
    keys = {}
    for entry in salts.TABLE:
        mol = Chem.MolFromSmiles(entry.smiles)
        assert mol is not None, entry.name
        key = Chem.MolToInchiKey(mol)
        assert key not in keys, f"{entry.name} duplicates {keys.get(key)}"
        keys[key] = entry.name
        assert entry.group in salts.GROUPS


@pytest.mark.parametrize("entry", salts.TABLE, ids=lambda e: e.name)
def test_each_entry_draws_and_reads_back(entry):
    sketch = salts.sketch_of(entry)
    assert convert.to_smiles(sketch) == Chem.CanonSmiles(entry.smiles)
    assert convert.bond_length_spread(sketch) < 0.05


def test_dbu_is_the_amidine_not_ketchers_hydrazine():
    mol = Chem.MolFromSmiles(BY_NAME["1,8-diazabicyclo[5.4.0]undec-7-ene"].smiles)
    assert rdMolDescriptors.CalcMolFormula(mol) == "C9H16N2"
    assert not mol.HasSubstructMatch(Chem.MolFromSmarts("[#7]-[#7]"))
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[NX2])[NX3]"))  # amidine


def test_guanidinium_carries_the_charge_on_the_cn3_centre():
    mol = Chem.MolFromSmiles(BY_NAME["guanidine hydrochloride"].smiles)
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3][CX3](=[NX3+])[NX3]"))
    assert not mol.HasSubstructMatch(Chem.MolFromSmarts("[NX4+]"))


def test_cis_decalin_is_cis():
    mol = Chem.AddHs(Chem.MolFromSmiles(BY_NAME["cis-decalin"].smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    heads = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C" and a.GetDegree() == 4
             and sum(n.GetSymbol() == "C" for n in a.GetNeighbors()) == 3]  # fmt: skip
    a, b = heads
    ha = next(n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetSymbol() == "H")
    hb = next(n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetSymbol() == "H")
    # cis: the bridgehead hydrogens are gauche or eclipsed, never anti (about 180 degrees)
    assert abs(transforms.GetDihedralDeg(mol.GetConformer(), ha, a, b, hb)) < 90


def test_place_puts_the_entry_to_the_right_of_the_drawing():
    history = History(convert.from_text("c1ccccc1C(=O)O"))
    right = history.sketch.bounds()[2]
    change = history.change("Add")
    added = salts.place(change, BY_NAME["sodium chloride"])
    history.commit(change)
    assert len(added) == 2
    assert min(history.sketch.atoms[a].x for a in added) > right
    assert convert.to_smiles(history.sketch) == Chem.CanonSmiles("O=C(O)c1ccccc1.[Na+].[Cl-]")
