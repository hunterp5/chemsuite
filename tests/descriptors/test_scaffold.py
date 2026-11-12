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
"""The Bemis-Murcko scaffold, and its stand-in where RDKit's fails (2026.9.2 on macOS)."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.descriptors.scaffold import murcko, pruned

#: molecules and RDKit's own scaffolds (GetScaffoldForMol, 2026.9.2 on Windows)
CASES = [
    ("CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1"),
    ("CS(=O)(=O)c1ccc(C)cc1", "c1ccccc1"),
    ("Cn1ccc2ccccc21", "c1ccc2[nH]ccc2c1"),
    ("CCO", ""),
    ("CN(C)CCC=C1c2ccccc2CCc2ccccc21", "C=C1c2ccccc2CCc2ccccc21"),
    ("C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O",
     "c1ccc2c(c1)CC[C@@H]1[C@@H]2CCC2CCC[C@H]21"),
    ("CC(C)OC(=O)C(C(=O)OC(C)C)=C1SC=CS1", "C=C1SC=CS1"),
    ("CN(C)c1ccc(C(=C2C=CC(=[N+](C)C)C=C2)c2ccccc2)cc1",
     "[N+]=C1C=CC(=C(c2ccccc2)c2ccccc2)C=C1"),
]  # fmt: skip


@pytest.mark.parametrize(("smiles", "scaffold"), CASES)
def test_the_stand_in_gives_rdkits_scaffold(smiles: str, scaffold: str) -> None:
    mol = Chem.MolFromSmiles(smiles)
    assert Chem.MolToSmiles(pruned(mol)) == scaffold
    assert Chem.MolToSmiles(murcko(mol)) == scaffold
