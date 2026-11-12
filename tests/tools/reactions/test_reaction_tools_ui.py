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

"""Reactions: Reaction Search, Split Reaction, Enumerate Products, and the
reaction filter card (REACTIONS_PLAN R4)."""

from __future__ import annotations

import pytest

from tests.gui import gui

pytestmark = gui

REACTIONS = [
    "CC(=O)O.NCc1ccccc1>CN(C)C=O>CC(=O)NCc1ccccc1",
    "Brc1ccccc1.OB(O)c1ccccc1>[Pd]>c1ccc(-c2ccccc2)cc1",
    "CC(C)(C)OC(=O)N1CCNCC1>>C1CNCCN1",
]
AMIDE = "[C:1](=[O:2])[OH].[N;H2,H1;!$(NC=O):3]>>[C:1](=[O:2])[N:3]"


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["R1", "R2", "R3"],
            "Rxn": REACTIONS,
            "Acid": ["CC(=O)O", "OC(=O)c1ccccc1", "CCO"],
            "Amine": ["NC", "NCC", "NC"],
        },
        undoable=False,
    )
    table.convert_column("Rxn", "reaction")
    table.convert_column("Acid", "structure")
    table.convert_column("Amine", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_reaction_tools_need_a_reaction_column(chem) -> None:
    actions = chem._tool_actions
    assert actions["reaction-search"].isEnabled() and actions["split-reaction"].isEnabled()
    chem.workspace.table.convert_column("Rxn", "text")
    chem._update_actions()
    assert not actions["reaction-search"].isEnabled()
    assert actions["enumerate-products"].isEnabled()  # it needs structures, not reactions


def test_reaction_search_selects_matches(chem, qtbot) -> None:
    dialog = chem.open_tool("reaction-search")
    assert "Enter a query" in dialog.status.text() and not dialog.primary.isEnabled()
    dialog.query.set_text("C(=O)O.N>>C(=O)N")
    assert dialog.primary.isEnabled() and dialog.primary.text() == "Search 3 Rows"
    dialog.run()
    first = int(chem.workspace.table.oids[0])
    qtbot.waitUntil(lambda: chem.workspace.selection.oids == {first}, timeout=20_000)


def test_reaction_search_marks_a_column_by_role(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("reaction-search")
    dialog.role.setCurrentIndex(dialog.role.findData("agent"))
    dialog.query.set_text("C>>C")
    assert "choose “Reaction contains”" in dialog.status.text()
    dialog.query.set_text("[Pd]")
    dialog.output.set_value("column")
    dialog.run()
    qtbot.waitUntil(lambda: "Reaction Search" in table.column_names, timeout=20_000)
    assert table.column("Reaction Search").displays(table.slots_of(table.oids)) == [
        "No",
        "Yes",
        "No",
    ]


def test_split_reaction_by_role(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("split-reaction")
    assert dialog.primary.text() == "Split 3 Reactions"
    dialog.run()
    qtbot.waitUntil(lambda: "Rxn Products" in table.column_names, timeout=20_000)
    assert table.kind("Rxn Reactants") == "structure"
    r1, r2, r3 = (int(o) for o in table.oids)
    assert table.cell("Rxn Reactants", r1) == "CC(=O)O.NCc1ccccc1"
    assert table.cell("Rxn Agents", r2) == "[Pd]"
    assert table.cell("Rxn Agents", r3) == ""
    assert table.cell("Rxn Products", r3) == "C1CNCCN1"
    table.undo()
    assert "Rxn Products" not in table.column_names


def test_split_reaction_one_column_per_component(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("split-reaction")
    dialog.layout_choice.set_value("each")
    dialog.run()
    qtbot.waitUntil(lambda: "Rxn Reactant 2" in table.column_names, timeout=20_000)
    r1, _r2, r3 = (int(o) for o in table.oids)
    assert table.cell("Rxn Reactant 2", r1) == "NCc1ccccc1"
    assert table.cell("Rxn Reactant 2", r3) == ""


def test_enumerate_row_by_row(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("enumerate-products")
    dialog.template.set_text("CC>>CO")
    assert "Map the atoms" in dialog.status.text()
    dialog.template.set_text(AMIDE)
    dialog.sources[0].setCurrentText("Acid")
    dialog.sources[1].setCurrentText("Amine")
    assert dialog._blocks_section.is_row_visible(dialog.sources[1])
    assert not dialog._blocks_section.is_row_visible(dialog.sources[2])
    dialog.mode.set_value("rows")
    assert dialog.primary.text() == "Enumerate 3 Rows"
    dialog.run()
    qtbot.waitUntil(lambda: "Product" in table.column_names, timeout=20_000)
    r1, r2, r3 = (int(o) for o in table.oids)
    assert table.cell("Product", r1) == "CNC(C)=O"
    assert table.cell("Product", r2) == "CCNC(=O)c1ccccc1"
    assert table.cell("Product", r3) == ""  # ethanol is not an acid
    assert table.cell("Products Found", r3) == "0"


def test_enumerate_every_combination_opens_a_library(chem, qtbot) -> None:
    dialog = chem.open_tool("enumerate-products")
    dialog.template.set_text(AMIDE)
    dialog.sources[0].setCurrentText("Acid")
    dialog.sources[1].setCurrentText("Amine")
    # 3 distinct acids × 2 distinct amines
    assert dialog.primary.text() == "Enumerate 6 Combinations"
    opened = []
    chem.open_table = lambda table, message="": opened.append((table, message))
    dialog.run()
    qtbot.waitUntil(lambda: bool(opened), timeout=20_000)
    library, message = opened[0]
    assert "4 products from 6 combinations" in message
    assert library.column_names == ["Reactant 1", "Reactant 2", "Product"]
    assert library.cell("Product", int(library.oids[0])) == "CNC(C)=O"


def test_enumerate_from_a_building_block_file(chem, qtbot, tmp_path) -> None:
    path = tmp_path / "amines.smi"
    path.write_text("NC methylamine\nNCC ethylamine\nC1CCNCC1 piperidine\n")
    dialog = chem.open_tool("enumerate-products")
    dialog.template.set_text(AMIDE)
    dialog.sources[0].setCurrentText("Acid")
    dialog.set_file(1, path)
    assert dialog.primary.text() == "Enumerate"  # the count is known once the file is read
    dialog.mode.set_value("rows")
    assert "choose columns" in dialog.status.text()
    dialog.mode.set_value("all")
    opened = []
    chem.open_table = lambda table, message="": opened.append((table, message))
    dialog.run()
    qtbot.waitUntil(lambda: bool(opened), timeout=20_000)
    assert "6 products from 9 combinations" in opened[0][1]


def test_reaction_filter_card(chem, qtbot) -> None:
    table = chem.workspace.table
    card = chem.filters.add_filter("Rxn")
    assert card.kind == "reaction"
    card.role.setCurrentIndex(card.role.findData("product"))
    card.query.set_text("c1ccccc1-c1ccccc1")
    card.search_now()
    second = int(table.oids[1])
    qtbot.waitUntil(lambda: card.hits == frozenset({second}), timeout=20_000)
    assert card.body_state() == {"query": "c1ccccc1-c1ccccc1", "role": "product"}
