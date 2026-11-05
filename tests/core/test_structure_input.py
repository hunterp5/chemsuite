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

import pytest
from rdkit import Chem

from chemsuite.core.structure_input import molblock_of, read_structure


@pytest.mark.parametrize(
    ("text", "kind", "valid"),
    [
        ("CCO", "smiles", True),
        (" c1ccccc1 \n", "smiles", True),
        ("InChI=1S/CH4/h1H4", "inchi", True),
        ("", "blank", True),
        ("   ", "blank", True),
        ("C1CC", "smiles", False),
        ("InChI=1S/nonsense", "inchi", False),
    ],
)
def test_kinds_and_validity(text: str, kind: str, valid: bool) -> None:
    parsed = read_structure(text)
    assert parsed.kind == kind
    assert parsed.is_valid is valid
    assert bool(parsed.error) is not valid


def test_errors_explain_the_problem() -> None:
    assert "syntax" in read_structure("C1CC(").error
    assert read_structure("CN(C)(C)(C)C").error == "Atom 2 (N) has too many bonds."
    assert "Aromatic ring" in read_structure("c1cccc1").error


def test_molfile_keeps_coordinates_and_displays_smiles() -> None:
    block = molblock_of(Chem.MolFromSmiles("OC(=O)c1ccccc1"))
    parsed = read_structure(block)
    assert parsed.kind == "molfile" and parsed.display == "O=C(O)c1ccccc1"
    column = parsed.as_column()
    assert column.display(0) == "O=C(O)c1ccccc1"
    assert column.mol(0).GetNumConformers() == 1


def test_smiles_kept_as_typed_with_formula_and_weight() -> None:
    parsed = read_structure("  OCC ")
    assert parsed.display == "OCC"
    assert parsed.canonical_smiles() == "CCO"
    assert parsed.formula == "C2H6O"
    assert parsed.mol_weight == pytest.approx(46.07, abs=0.01)
    assert read_structure("").as_column().display(0) == ""
