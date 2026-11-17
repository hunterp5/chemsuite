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

from __future__ import annotations

import io

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from chemsuite.sketch import convert, files, reaction


def test_a_structure_saves_as_a_molfile_sd_file_or_smiles() -> None:
    sketch = convert.from_text("CC(=O)Oc1ccccc1C(=O)O")
    assert list(files.file_types(sketch)) == [".mol", ".sdf", ".smi"]
    aspirin = Chem.CanonSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mol = Chem.MolFromMolBlock(files.file_text(sketch, ".mol"))
    assert Chem.MolToSmiles(mol) == aspirin
    sdf = files.file_text(sketch, ".SDF")
    assert sdf.endswith("M  END\n$$$$\n")
    read = next(Chem.ForwardSDMolSupplier(io.BytesIO(sdf.encode())))
    assert Chem.MolToSmiles(read) == aspirin
    assert files.file_text(sketch, ".smi") == aspirin + "\n"
    with pytest.raises(ValueError):
        files.file_text(sketch, ".rxn")


def test_a_reaction_saves_as_an_rxn_file_or_reaction_smiles() -> None:
    sketch = reaction.from_text("CCO.CC(=O)O>>CC(=O)OCC")
    assert list(files.file_types(sketch)) == [".rxn", ".smi"]
    rxn = rdChemReactions.ReactionFromRxnBlock(files.file_text(sketch, ".rxn"))
    assert rxn.GetNumReactantTemplates() == 2 and rxn.GetNumProductTemplates() == 1
    assert ">>" in files.file_text(sketch, ".smi")
    with pytest.raises(ValueError):
        files.file_text(sketch, ".mol")
