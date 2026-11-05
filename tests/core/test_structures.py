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

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.core.mols import mol_from_text, parse_structure, parse_structures
from chemsuite.core.table import TEXT, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import (
    STRUCTURE,
    NumberColumn,
    StructureColumn,
    TextColumn,
    concat_columns,
)

from .conftest import DRUGS, molblock


def test_text_parsing_accepts_smiles_inchi_and_molfile() -> None:
    assert mol_from_text("c1ccccc1").GetNumAtoms() == 6
    inchi = Chem.MolToInchi(Chem.MolFromSmiles("CCO"))
    assert Chem.MolToSmiles(mol_from_text(inchi)) == "CCO"
    assert mol_from_text(molblock("CCO")).GetNumAtoms() == 3
    assert mol_from_text("C1CC") is None
    assert mol_from_text("not a molecule") is None
    assert mol_from_text("") is None


def test_molfile_with_empty_title_line_parses() -> None:
    block = molblock("CCO", title="")
    assert block.startswith("\n")
    assert mol_from_text(block) is not None


def test_parse_structure_keeps_typed_smiles_and_flags_invalid() -> None:
    blob, text = parse_structure(" CCO ")
    assert blob and text == "CCO"
    blob, text = parse_structure(molblock("OCC"))
    assert blob and text == "CCO"  # molfiles display as canonical SMILES
    assert parse_structure("C1CC") == (b"", "C1CC")
    assert parse_structures(["CCO", "", "xx"])[1:] == [b"", b""]


def test_blob_keeps_coordinates_and_stereo_but_not_properties() -> None:
    mol = Chem.MolFromMolBlock(molblock(DRUGS["alanine"]))
    mol.SetProp("pIC50", "7")
    column = StructureColumn.from_parsed([parse_structure(Chem.MolToMolBlock(mol))[0]], ["x"])
    back = column.mol(0)
    assert back.GetNumConformers() == 1
    assert Chem.MolToSmiles(back) == "C[C@H](N)C(=O)O"
    assert not back.HasProp("pIC50")


def test_cell_states() -> None:
    column = StructureColumn.from_texts(["CCO", "", "C1CC"])
    assert column.pending_mask().tolist() == [True, False, True]
    assert column.blank_mask().tolist() == [False, True, False]
    assert column.mol(0).GetNumAtoms() == 3  # parsed on first use
    assert column.mol(2) is None
    assert column.mol(1) is None
    assert column.pending_mask().tolist() == [False, False, False]
    assert column.invalid_mask().tolist() == [False, False, True]
    assert column.displays(np.arange(3)) == ["CCO", "", "C1CC"]


def test_molecules_are_cached() -> None:
    column = StructureColumn.from_texts(["c1ccccc1"])
    assert column.mol(0) is column.mol(0)


def test_typing_smiles_into_a_structure_cell_is_undoable() -> None:
    table = Table()
    table.add_rows({"Structure": StructureColumn.from_texts(["CCO", "CCN"])})
    table.set_cells("Structure", [1], ["c1ccccc1"])
    column = table.column("Structure")
    assert table.kind("Structure") == STRUCTURE
    assert column.mol(0).GetNumAtoms() == 6
    table.undo()
    assert table.cell("Structure", 1) == "CCO"
    assert table.column("Structure").mol(0).GetNumAtoms() == 3


def test_structure_columns_cannot_be_sorted_but_can_be_filtered() -> None:
    from chemsuite.core.table import TextFilter

    table = Table()
    table.add_rows({"S": StructureColumn.from_texts(["CCO", "c1ccccc1"])})
    with pytest.raises(ValueError):
        table.set_sort([("S", True)])
    table.set_filter("f", TextFilter("S", "c1"))
    assert table.visible_oids().tolist() == [2]


def test_fill_structures_is_not_an_edit_but_notifies() -> None:
    table = Table()
    table.add_rows({"S": StructureColumn.from_texts(["CCO", "xx"])}, undoable=False)
    seen: list[object] = []
    table.subscribe(seen.append)
    table.fill_structures("S", np.array([1, 2]), parse_structures(["CCO", "xx"]))
    assert not table.history.can_undo
    assert isinstance(seen[-1], ev.Changed) and not seen[-1].structural
    assert table.column("S").invalid_mask().tolist() == [False, True]


def test_results_can_create_structure_columns() -> None:
    table = Table()
    table.add_rows({"ID": ["a", "b"]})
    names = table.apply_results({"Core": ([2], StructureColumn.from_texts(["c1ccccc1"]))})
    assert names == {"Core": "Core"}
    assert table.kind("Core") == STRUCTURE
    assert table.column("Core").displays(np.arange(2)) == ["", "c1ccccc1"]


def test_concat_reconciles_kinds() -> None:
    numbers = NumberColumn.from_floats(np.array([1.0, 2.0]))
    blank_text = TextColumn.blanks(2)
    text = TextColumn.from_texts(["x", "y"])
    assert concat_columns([numbers, blank_text]).kind == "number"
    mixed = concat_columns([numbers, text])
    assert mixed.kind == TEXT
    assert mixed.displays(np.arange(4)) == ["1", "2", "x", "y"]
    structures = concat_columns([StructureColumn.from_texts(["C"]), TextColumn.from_texts(["N"])])
    assert structures.kind == STRUCTURE and structures.displays(np.arange(2)) == ["C", "N"]
