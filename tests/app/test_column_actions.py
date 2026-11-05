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
from PySide6.QtCore import QItemSelection, QItemSelectionModel, Qt
from PySide6.QtWidgets import QInputDialog, QMessageBox

from .conftest import gui

pytestmark = gui


def answer_text(monkeypatch, text: str, ok: bool = True) -> None:
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: (text, ok))


def answer_question(monkeypatch, button) -> None:
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: button)


def menu_titles(menu) -> list[str]:
    return [a.text() for a in menu.actions() if not a.isSeparator()]


def trigger(menu, title: str) -> None:
    for action in menu.actions():
        if action.text() == title:
            action.trigger()
            return
        if action.menu() is not None:
            for sub in action.menu().actions():
                if sub.text() == title:
                    sub.trigger()
                    return
    raise AssertionError(f"no menu item {title!r} in {menu_titles(menu)}")


def test_menu_offers_then_sort_only_with_an_existing_sort(window) -> None:
    table = window.workspace.table
    menu = window.columns.build_menu("pIC50", [], 1)
    assert "Then Sort Ascending" not in menu_titles(menu)
    assert not next(a for a in menu.actions() if a.text() == "Clear Sort").isEnabled()
    trigger(menu, "Sort Descending")
    assert table.sort_keys == (("pIC50", False),)
    trigger(window.columns.build_menu("ID", [], 0), "Then Sort Ascending")
    assert table.sort_keys == (("pIC50", False), ("ID", True))
    header = window.view.horizontalHeader()
    assert header.sortIndicatorSection() == 1  # primary key shown, secondary kept
    assert header.sortIndicatorOrder() == Qt.SortOrder.DescendingOrder


def test_rename_insert_delete_are_undoable(window, monkeypatch) -> None:
    table = window.workspace.table
    answer_text(monkeypatch, "Potency")
    trigger(window.columns.build_menu("pIC50", [], 1), "Rename…")
    assert table.column_names == ["ID", "Potency"]
    answer_text(monkeypatch, "Notes")
    trigger(window.columns.build_menu("ID", [], 0), "Insert Column Right…")
    assert table.column_names == ["ID", "Notes", "Potency"]
    answer_text(monkeypatch, "ID")  # duplicate name: refused with a warning
    warnings: list[str] = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: warnings.append(a[2]))
    trigger(window.columns.build_menu("ID", [], 0), "Insert Column Left…")
    assert warnings and table.column_names == ["ID", "Notes", "Potency"]
    table.undo()
    table.undo()
    assert table.column_names == ["ID", "pIC50"]


def test_delete_uses_fully_selected_columns(window) -> None:
    table = window.workspace.table
    table.add_column("Notes")
    view = window.view
    model = view.table_model
    selection = QItemSelection()
    for col in (0, 2):
        selection.select(model.index(0, col), model.index(model.rowCount() - 1, col))
    view.selectionModel().select(selection, QItemSelectionModel.SelectionFlag.ClearAndSelect)
    assert view.selected_column_names() == ["ID", "Notes"]
    menu = window.columns.build_menu("Notes", view.selected_column_names(), 2)
    trigger(menu, "Delete 2 Columns")
    assert table.column_names == ["pIC50"]
    table.undo()
    assert table.column_names == ["ID", "pIC50", "Notes"]
    # right-clicking a column outside the selection deletes just that column
    menu = window.columns.build_menu("pIC50", ["ID", "Notes"], 1)
    assert "Delete Column" in menu_titles(menu)


@pytest.mark.parametrize(
    ("button", "converted"),
    [
        (QMessageBox.StandardButton.Yes, True),
        (QMessageBox.StandardButton.Cancel, False),
    ],
)
def test_convert_to_number_asks_before_clearing(window, monkeypatch, button, converted) -> None:
    table = window.workspace.table
    table.add_column("Mixed", ["1.5", "n/a", "3"])
    answer_question(monkeypatch, button)
    trigger(window.columns.build_menu("Mixed", [], 2), "Number")
    assert table.kind("Mixed") == ("number" if converted else "text")
    if converted:
        assert [table.cell("Mixed", o) for o in (1, 2, 3)] == ["1.5", "", "3"]


def test_convert_to_structure_parses_in_background(qtbot, window) -> None:
    table = window.workspace.table
    table.add_column("SMILES", ["CCO", "C1CC", "c1ccccc1"])
    default_height = window.view.verticalHeader().defaultSectionSize()
    trigger(window.columns.build_menu("SMILES", [], 2), "Structure")
    column = table.column("SMILES")
    qtbot.waitUntil(lambda: not column.pending_mask().any(), timeout=20_000)
    assert column.invalid_mask().tolist() == [False, True, False]
    assert window.view.verticalHeader().defaultSectionSize() > default_height
    assert table.history.undo_label == "Convert SMILES to structure"


def test_pasted_smiles_in_a_structure_column_get_parsed(qtbot, window) -> None:
    from chemsuite.core.table.columns import StructureColumn

    table = window.workspace.table
    table.add_column("Structure", StructureColumn.from_texts(["", "", ""]))
    table.set_cells("Structure", [2], ["c1ccncc1"])
    column = table.column("Structure")
    qtbot.waitUntil(lambda: not column.pending_mask().any(), timeout=20_000)
    assert table.column("Structure").mol(1) is not None


def test_widths_survive_resets_and_header_drag_moves_columns(qtbot, window) -> None:
    table = window.workspace.table
    view = window.view
    view.setColumnWidth(0, 237)
    table.add_rows({"ID": ["D"], "pIC50": ["8"]})  # model reset
    assert view.columnWidth(0) == 237
    view.horizontalHeader().moveSection(0, 1)  # drag "ID" to the right
    qtbot.waitUntil(lambda: table.column_names == ["pIC50", "ID"], timeout=2000)
    assert view.horizontalHeader().visualIndex(0) == 0  # view order follows the table
    assert view.columnWidth(1) == 237  # the width moved with the column
    table.undo()
    assert table.column_names == ["ID", "pIC50"]


def test_clicking_an_unsortable_header_keeps_the_indicator(qtbot, window) -> None:
    from chemsuite.core.table.columns import StructureColumn

    table = window.workspace.table
    table.add_column("Structure", StructureColumn.from_texts(["", "", ""]), position=0)
    table.set_sort([("ID", True)])
    header = window.view.horizontalHeader()
    assert header.sortIndicatorSection() == 1
    header.setSortIndicator(0, Qt.SortOrder.AscendingOrder)  # what a click does
    qtbot.waitUntil(lambda: header.sortIndicatorSection() == 1, timeout=2000)
    assert table.sort_keys == (("ID", True),)
