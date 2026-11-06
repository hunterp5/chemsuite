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

"""Reaction columns (REACTIONS_PLAN R3): parsing, the table, documents and files."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.core import document
from chemsuite.core.io import rdf
from chemsuite.core.io.load import load_table, parse_pending_structures, save_table
from chemsuite.core.reactions import (
    parse_reaction,
    reaction_from_text,
    reaction_smiles,
    rxn_block,
)
from chemsuite.core.table import Table
from chemsuite.core.table.columns import (
    REACTION,
    STRUCTURE,
    TEXT,
    ReactionColumn,
    StructureColumn,
    concat_columns,
)

ESTER = "CC(=O)O.OCC>[H+]>CC(=O)OCC"
AMIDE = "[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>CN(C)C=O>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]"
SUZUKI = "Brc1ccccc1.OB(O)c1ccccc1>[Pd]>c1ccc(-c2ccccc2)cc1"


# -- parsing ----------------------------------------------------------------------------------


@pytest.mark.parametrize("text", [ESTER, AMIDE, SUZUKI, "C>>", "[#6:1]=[O:2]>>[#6:1][O:2]"])
def test_valid_reactions_parse_and_keep_their_text(text: str) -> None:
    blob, shown = parse_reaction(text)
    assert blob and shown == text


@pytest.mark.parametrize("text", ["CCO", "x>>y", "C(C)(C)(C)(C)C>>C", ">>", "", "C>C"])
def test_invalid_reactions_do_not_parse(text: str) -> None:
    assert parse_reaction(text)[0] == b""


def test_an_rxn_block_is_shown_as_reaction_smiles() -> None:
    block = rxn_block(reaction_from_text(AMIDE))
    assert block.startswith("$RXN")
    blob, shown = parse_reaction(block)
    assert blob and shown == reaction_smiles(reaction_from_text(AMIDE))


def test_reaction_column_parses_on_demand_and_flags_invalid_cells() -> None:
    column = ReactionColumn.from_texts([ESTER, "", "not a reaction"])
    assert column.kind == REACTION and not column.sortable
    assert column.pending_mask().tolist() == [True, False, True]
    rxn = column.reaction(0)
    assert rxn.GetNumReactantTemplates() == 2 and rxn.GetNumAgentTemplates() == 1
    assert column.reaction(1) is None and column.reaction(2) is None
    assert column.invalid_mask().tolist() == [False, False, True]
    assert column.reaction(0) is rxn  # cached


def test_structures_and_reactions_are_separate_kinds() -> None:
    reactions = ReactionColumn.from_texts([ESTER])
    structures = StructureColumn.from_texts(["CCO"])
    assert not isinstance(reactions, StructureColumn)
    with pytest.raises(TypeError):
        structures.put(np.array([0]), reactions)
    assert concat_columns([reactions, ReactionColumn.blanks(2)]).kind == REACTION
    assert concat_columns([reactions, structures]).kind == TEXT  # no cell is lost


# -- the table ----------------------------------------------------------------------------------


def test_convert_a_text_column_to_reactions_and_edit_it() -> None:
    table = Table()
    table.add_rows({"rxn": [ESTER, SUZUKI]}, undoable=False)
    table.convert_column("rxn", REACTION)
    assert table.kind("rxn") == REACTION
    assert parse_pending_structures(table) == 0
    table.set_cells("rxn", [int(table.oids[0])], ["C=C>>CC"])  # typed text stays a reaction
    assert table.kind("rxn") == REACTION
    assert table.column("rxn").reaction(0).GetNumProductTemplates() == 1
    table.undo()
    assert table.cell("rxn", int(table.oids[0])) == ESTER


def test_reaction_columns_save_in_documents(tmp_path) -> None:
    table = Table()
    table.add_rows({"ID": ["a", "b", "c"], "rxn": [ESTER, "", "bad>>"]}, undoable=False)
    table.convert_column("rxn", REACTION)
    parse_pending_structures(table)
    path = tmp_path / "rxn.csd"
    document.save_document(table, path)
    loaded, _ = document.load_document(path)
    column = loaded.column("rxn")
    assert isinstance(column, ReactionColumn)
    assert column.displays(np.arange(3)) == [ESTER, "", "bad>>"]
    assert column.invalid_mask().tolist() == [False, False, True]
    assert not column.pending_mask().any()  # parse results are saved too


# -- table files ---------------------------------------------------------------------------------


def test_csv_reaction_smiles_column_is_detected(tmp_path) -> None:
    path = tmp_path / "uspto.csv"
    path.write_text(
        f"id,prod_smiles,rxn_smiles\nr1,CCOC(C)=O,{ESTER}\nr2,c1ccc(-c2ccccc2)cc1,{SUZUKI}\n"
    )
    table = load_table(path)
    # rxn_smiles is named like a structure column, but holds reactions
    assert table.kind("rxn_smiles") == REACTION
    assert table.kind("prod_smiles") == STRUCTURE
    assert parse_pending_structures(table) == 0


def test_a_column_named_as_structures_on_purpose_stays_structures(tmp_path) -> None:
    path = tmp_path / "r.csv"
    path.write_text(f"rxn_smiles\n{ESTER}\n")
    assert load_table(path, structure_column="rxn_smiles").kind("rxn_smiles") == STRUCTURE


def test_a_column_with_some_non_reactions_stays_text(tmp_path) -> None:
    path = tmp_path / "r.csv"
    path.write_text(f"note\n{ESTER}\nsee above\n")
    assert load_table(path).kind("note") == TEXT


def test_reaction_columns_round_trip_through_csv_and_xlsx(tmp_path) -> None:
    table = Table()
    table.add_rows({"ID": ["a", "b"], "Reaction": [AMIDE, SUZUKI]}, undoable=False)
    table.convert_column("Reaction", REACTION)
    for name in ("r.csv", "r.xlsx"):
        save_table(table, tmp_path / name)
        back = load_table(tmp_path / name)
        assert back.kind("Reaction") == REACTION
        assert back.column("Reaction").displays(np.arange(2)) == [AMIDE, SUZUKI]


# -- RXN and RD files ------------------------------------------------------------------------------


def _reaction_table() -> Table:
    table = Table()
    table.add_rows(
        {
            "Reaction": [AMIDE, ESTER, SUZUKI],
            "Yield": ["85", "", "92.5"],
            "Notes": ["two lines\nof notes", "x" * 120, "plain"],
        },
        undoable=False,
    )
    table.convert_column("Reaction", REACTION)
    return table


def test_rd_file_round_trip(tmp_path) -> None:
    table = _reaction_table()
    path = tmp_path / "set.rdf"
    save_table(table, path)
    text = path.read_text()
    assert text.startswith("$RDFILE 1\n$DATM")
    assert text.count("$RFMT\n$RXN") == 3
    assert "$DTYPE Yield\n$DATUM 85\n" in text
    back = load_table(path)
    assert back.kind("Reaction") == REACTION
    assert parse_pending_structures(back) == 0
    for i in range(3):
        before = reaction_smiles(table.column("Reaction").reaction(i))
        assert back.column("Reaction").display(i) == before
    assert back.column("Yield").displays(np.arange(3)) == ["85", "", "92.5"]
    assert back.column("Notes").displays(np.arange(3)) == table.column("Notes").displays(
        np.arange(3)
    )


def test_rd_file_keeps_atom_maps_and_agents(tmp_path) -> None:
    table = _reaction_table()
    save_table(table, tmp_path / "set.rdf")
    rxn = load_table(tmp_path / "set.rdf").column("Reaction").reaction(0)
    assert rxn.GetNumAgentTemplates() == 1
    maps = sorted(a.GetAtomMapNum() for a in rxn.GetReactantTemplate(0).GetAtoms())
    assert maps == [1, 2, 3, 4]


def test_rd_file_data_wrapped_at_80_characters_is_joined(tmp_path) -> None:
    block = rxn_block(reaction_from_text(ESTER))
    long = "A" * 100
    path = tmp_path / "wrapped.rdf"
    path.write_text(
        "$RDFILE 1\n$DATM    01/01/26 00:00\n"
        f"$RFMT $RIREG 7\n{block}\n"
        f"$DTYPE RXN:VARIATION(1):LONG\n$DATUM {long[:73]}\n{long[73:]}\n"
        "$DTYPE Two\n$DATUM first\nsecond\n"
    )
    table = load_table(path)
    assert table.cell("RXN:VARIATION(1):LONG", int(table.oids[0])) == long
    assert table.cell("Two", int(table.oids[0])) == "first\nsecond"
    assert table.cell("ID", int(table.oids[0])) == "7"


def test_unreadable_reaction_record_is_kept_with_its_data(tmp_path) -> None:
    path = tmp_path / "bad.rdf"
    path.write_text(
        "$RDFILE 1\n$RFMT\n$RXN\n\n\n\n  1  0\n$MOL\nbroken\n$DTYPE Name\n$DATUM kept\n"
    )
    table = load_table(path)
    assert len(table) == 1
    assert table.column("Reaction").invalid_mask().tolist() == [True]
    assert table.cell("Name", int(table.oids[0])) == "kept"


def test_rxn_file_open_and_save(tmp_path) -> None:
    path = tmp_path / "one.rxn"
    path.write_text(rxn_block(reaction_from_text(AMIDE)) + "\n")
    table = load_table(path)
    assert len(table) == 1 and table.kind("Reaction") == REACTION
    assert table.column("Reaction").reaction(0).GetNumProductTemplates() == 1
    save_table(table, tmp_path / "again.rxn")
    assert (tmp_path / "again.rxn").read_text().startswith("$RXN")
    with pytest.raises(ValueError, match="one reaction"):
        save_table(_reaction_table(), tmp_path / "many.rxn")


def test_rd_file_needs_a_reaction_column(tmp_path) -> None:
    table = Table()
    table.add_rows({"ID": ["a"]}, undoable=False)
    with pytest.raises(ValueError, match="no reaction column"):
        rdf.write_rdf(table, tmp_path / "x.rdf")
