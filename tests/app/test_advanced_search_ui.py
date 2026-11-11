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
from PySide6.QtCore import Qt

from .conftest import gui

pytestmark = gui


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["C0", "C1", "C2", "C3", "C4"],
            "SMILES": ["c1ccccc1O", "CCO", "Cc1ccccc1", "c1ccncc1", "CCCC"],
            "pIC50": ["7.5", "6.0", "8.2", "5.1", ""],
            "Series": ["A", "B", "A", "B", "A"],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.convert_column("pIC50", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def shown(window) -> list[str]:
    table = window.workspace.table
    return [table.cell("ID", int(o)) for o in table.oids[table.view()]]


def fill(row, column: str, op: str, value: str = "", value2: str = "") -> None:
    row.column.setCurrentIndex(row.column.findData(column) if column else 0)
    row.op.setCurrentIndex(row.op.findData(op))
    if op in ("substructure", "exact"):
        row.structure.set_text(value)
    else:
        row.value.setText(value)
        row.value2.setText(value2)


def test_show_only_with_all_and_any(chem) -> None:
    panel = chem.show_advanced_search()
    assert chem.search_dock.isFloating()  # its own window, never beside Filters
    fill(panel.rows[0], "pIC50", "ge", "7")
    fill(panel.add_row(), "Series", "equals", "A")
    assert panel.run("show")
    assert shown(chem) == ["C0", "C2"]
    assert panel.clear_button.isEnabled()
    panel.combine.setCurrentIndex(panel.combine.findData("any"))
    panel.run("show")
    assert shown(chem) == ["C0", "C2", "C4"]
    panel.rows[1].negate.setChecked(True)  # pIC50 ≥ 7 or not Series A
    panel.run("show")
    assert shown(chem) == ["C0", "C1", "C2", "C3"]
    # Text and number conditions follow edits without searching again.
    chem.workspace.table.set_cells("pIC50", [int(chem.workspace.table.oids[4])], ["9"])
    assert shown(chem) == ["C0", "C1", "C2", "C3", "C4"]
    chem.show_all_rows()
    assert len(shown(chem)) == 5 and panel.showing is None
    assert len(panel.recent) == 3


def test_invalid_rows_block_the_search(chem) -> None:
    panel = chem.show_advanced_search()
    fill(panel.rows[0], "pIC50", "between", "6", "lots")
    assert not panel.show_button.isEnabled()
    assert "not a number" in panel.rows[0].error.text()
    assert not panel.run("show")
    fill(panel.rows[0], "pIC50", "between", "6", "8")
    assert panel.show_button.isEnabled() and not panel.rows[0].error.isVisible()
    # Operators follow the column's kind.
    row = panel.rows[0]
    row.column.setCurrentIndex(row.column.findData("Series"))
    assert row.op.findData("between") < 0 and row.op.findData("regex") >= 0


def test_select_and_add_to_list(chem) -> None:
    panel = chem.show_advanced_search()
    fill(panel.rows[0], "", "contains", "b")  # any column
    assert panel.run("select")
    table = chem.workspace.table
    assert sorted(chem.workspace.selection.oids) == [int(table.oids[1]), int(table.oids[3])]
    assert not panel.list_button.isEnabled()  # no list name yet
    panel.list_name.setCurrentText("Series B")
    assert panel.run("list")
    assert chem.workspace.lists.oids("Series B").tolist() == [
        int(table.oids[1]),
        int(table.oids[3]),
    ]


def test_structure_condition(chem, qtbot) -> None:
    panel = chem.show_advanced_search()
    fill(panel.rows[0], "SMILES", "substructure", "c1ccccc1")
    fill(panel.add_row(), "pIC50", "gt", "8")
    assert panel.run("show")
    qtbot.waitUntil(lambda: panel.job is None and panel.showing is not None, timeout=20_000)
    assert shown(chem) == ["C2"]
    # Editing a structure searches again.
    table = chem.workspace.table
    table.set_cells("pIC50", [int(table.oids[0])], ["9"])
    assert shown(chem) == ["C0", "C2"]
    table.set_cells("SMILES", [int(table.oids[0])], ["CCN"])
    qtbot.waitUntil(lambda: shown(chem) == ["C2"], timeout=20_000)


def test_removing_a_column_stops_the_search(chem) -> None:
    panel = chem.show_advanced_search()
    fill(panel.rows[0], "Series", "equals", "A")
    panel.run("show")
    chem.workspace.table.remove_columns(["Series"])
    assert panel.showing is None and len(shown(chem)) == 5
    assert "changed" in panel.status.text()
    assert "missing" in panel.rows[0].column.currentText()


def test_saved_with_the_document(chem, tmp_path, qtbot) -> None:
    panel = chem.show_advanced_search()
    fill(panel.rows[0], "pIC50", "lt", "7")
    panel.run("show")
    target = tmp_path / "search.csd"
    assert chem._save_to(target)
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    chem._install(Workspace(Table()))
    assert not chem.advanced_search.recent  # a new document starts with no searches
    chem.open_path(target)
    qtbot.waitUntil(lambda: chem.workspace.path == target, timeout=10_000)
    restored = chem.advanced_search
    assert restored.current_search().criteria[0].value == "7"
    assert restored.recent and restored.recent[0] == restored.current_search()
    qtbot.waitUntil(lambda: shown(chem) == ["C1", "C3"], timeout=5000)


def test_advanced_search_is_its_own_window(chem, qtbot) -> None:
    from chemsuite.app.document import dock_layout

    chem.show()
    qtbot.waitExposed(chem)
    chem.filter_column("Series")
    dock = chem.search_dock
    chem.show_advanced_search()
    assert dock.isFloating() and dock.isVisible()
    assert dock.allowedAreas() == Qt.DockWidgetArea.NoDockWidgetArea
    header = dock.titleBarWidget()
    assert not header.float_button.isVisibleTo(dock) and not header.arrange_button.isVisibleTo(dock)
    assert dock not in dock_layout.docked(chem) and dock not in chem.tabifiedDockWidgets(
        chem.filter_dock
    )
    moved = dock.pos()
    dock.close()
    chem.show_advanced_search()  # reopening keeps where it was put
    assert dock.isFloating() and dock.pos() == moved


def test_nested_groups(chem, qtbot) -> None:
    panel = chem.show_advanced_search()
    # Series A and not (pIC50 ≥ 8 or SMILES contains phenol)
    fill(panel.rows[0], "Series", "equals", "A")
    group = panel.add_group()
    group.negate.setChecked(True)
    group.combine.setCurrentIndex(group.combine.findData("any"))
    fill(group.rows()[0], "pIC50", "ge", "8")
    fill(group.add_row(), "SMILES", "substructure", "c1ccccc1O")
    search = panel.current_search()
    assert search.describe() == (
        "Series equals A and not (pIC50 ≥ 8 or SMILES contains substructure c1ccccc1O)"
    )
    assert len(panel.rows) == 3
    assert panel.run("show")
    qtbot.waitUntil(lambda: shown(chem) == ["C4"], timeout=20_000)
    # an empty group blocks the search until it has a condition
    empty = panel.add_group()
    panel.root.remove(panel.rows[0])
    empty.remove(empty.items[0])
    assert not panel.show_button.isEnabled()
    panel.root.remove(empty)
    assert panel.show_button.isEnabled()


def test_not_in_one_column_and_not_in_another(chem) -> None:
    panel = chem.show_advanced_search()
    fill(panel.rows[0], "Series", "equals", "A")
    panel.rows[0].negate.setChecked(True)
    fill(panel.add_row(), "ID", "equals", "C3")
    panel.rows[1].negate.setChecked(True)
    assert panel.run("show")
    assert shown(chem) == ["C1"]


def test_groups_are_saved_with_the_document(chem, qtbot, tmp_path) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    panel = chem.show_advanced_search()
    fill(panel.rows[0], "Series", "equals", "B")
    group = panel.add_group()
    group.negate.setChecked(True)
    fill(group.rows()[0], "ID", "equals", "C3")
    saved = panel.current_search()
    target = tmp_path / "groups.csd"
    assert chem._save_to(target)
    chem._install(Workspace(Table()))
    chem.open_path(target)
    qtbot.waitUntil(lambda: chem.advanced_search.current_search() == saved, timeout=10_000)
