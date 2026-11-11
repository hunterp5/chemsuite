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

from .conftest import gui

pytestmark = gui


def ids_shown(window) -> list[str]:
    table = window.workspace.table
    return [table.cell("ID", int(o)) for o in table.visible_oids()]


def test_search_bar_filters_rows(window) -> None:
    bar = window.search
    bar.edit.setText("b")
    bar.apply_now()
    assert ids_shown(window) == ["B"]
    bar.edit.setText("[")
    bar.regex.setChecked(True)  # invalid regex: no filter, red box
    assert len(window.workspace.table.view()) == 3 and "Invalid" in bar.edit.toolTip()
    bar.edit.setText("^[AC]$")
    bar.apply_now()
    assert ids_shown(window) == ["A", "C"]
    bar.edit.clear()
    bar.apply_now()
    assert "search" not in window.workspace.table.filters
    assert window.workspace.table.history.is_clean


def test_find_next_walks_and_wraps(window) -> None:
    table = window.workspace.table
    table.add_column("Notes", ["x5", "", "5x"])
    dialog = window.show_find_replace()
    dialog.find_edit.setText("5")
    view = window.view
    seen = []
    for _ in range(4):
        assert dialog.find_next()
        current = view.currentIndex()
        seen.append((current.row(), view.table_model.column_name(current.column())))
    assert seen == [(0, "pIC50"), (0, "Notes"), (2, "Notes"), (0, "pIC50")]  # wraps
    assert dialog.status.text() == "1 of 3"
    assert dialog.find_next(backwards=True) and dialog.status.text() == "3 of 3"


def test_replace_and_replace_all(window) -> None:
    table = window.workspace.table
    dialog = window.show_find_replace()
    dialog.find_edit.setText("6")
    dialog.replace_edit.setText("6.5")
    dialog.replace_one()  # first press goes to the match
    assert table.cell("pIC50", 2) == "6"
    dialog.replace_one()  # second press replaces it
    assert table.cell("pIC50", 2) == "6.5"
    dialog.find_edit.setText("A")
    dialog.mode.setCurrentIndex(dialog.mode.findData("whole"))
    dialog.replace_edit.setText("Z")
    assert dialog.replace_everything() == 1
    assert table.cell("ID", 1) == "Z" and dialog.status.text() == "Replaced in 1 cell(s)"
    table.undo()
    assert table.cell("ID", 1) == "A"


def test_selected_columns_scope(window) -> None:
    table = window.workspace.table
    dialog = window.show_find_replace()
    view = window.view
    view.setCurrentIndex(view.table_model.index(0, 0))  # select a cell in ID only
    dialog.scope.setCurrentIndex(dialog.scope.findData("selected"))
    dialog.find_edit.setText("5")  # only in pIC50
    assert not dialog.find_next() and dialog.status.text() == "No matches"
    table.set_cells("ID", [3], ["C5"])
    assert dialog.find_next() and table.cell("ID", 3) == "C5"


def test_search_bar_starts_hidden_and_closing_clears_it(window, qtbot) -> None:
    from PySide6.QtCore import Qt

    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    bar = window.search
    action = window.command_action("view.search_bar")
    assert not bar.is_open and not action.isChecked()
    window.command_action("edit.search").trigger()  # Search Rows opens and focuses it
    assert bar.is_open and action.isChecked()
    bar.edit.setText("b")
    bar.apply_now()
    assert ids_shown(window) == ["B"]
    qtbot.keyClick(bar.edit, Qt.Key.Key_Escape)  # first Esc clears the text
    assert bar.is_open and bar.edit.text() == ""
    bar.edit.setText("b")
    bar.apply_now()
    action.trigger()  # View ▸ Search Bar off: hides and clears the search
    assert not bar.is_open and "search" not in window.workspace.table.filters
    assert ids_shown(window) == ["A", "B", "C"]
    action.trigger()
    assert bar.is_open
    qtbot.keyClick(bar.edit, Qt.Key.Key_Escape)  # Esc in an empty box closes the bar
    assert not bar.is_open and not action.isChecked()
    # A newly loaded document starts with the bar hidden too.
    window.search.focus()
    table = Table()
    table.add_rows({"ID": ["X"]}, undoable=False)
    window._install(Workspace(table))
    assert not window.search.is_open and not action.isChecked()


def test_search_bar_searches_one_column(window) -> None:
    table = window.workspace.table
    table.add_column("Notes", ["", "a note", ""])
    bar = window.search
    assert bar.column.itemText(0) == "All Columns" and bar.column_name() is None
    bar.edit.setText("a")
    bar.apply_now()
    assert ids_shown(window) == ["A", "B"]  # "A" in ID, "a note" in Notes
    bar.set_column("Notes")
    assert ids_shown(window) == ["B"] and "Notes" in bar.edit.placeholderText()
    table.rename_column("Notes", "Remarks")
    assert bar.column_name() == "Remarks" and ids_shown(window) == ["B"]
    table.remove_column("Remarks")  # the column goes: back to every column
    assert bar.column_name() is None and ids_shown(window) == ["A"]


def test_search_bar_substructure_search(window, qtbot) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["A", "B", "C"], "SMILES": ["c1ccccc1O", "CCO", "Oc1ccncc1"]}, undoable=False
    )
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    bar = window.search
    bar.set_column("SMILES")
    assert bar.structure_mode() and not bar.options_button.isEnabled()
    bar.edit.setText("c1ccccc1")
    bar.apply_now()
    qtbot.waitUntil(lambda: ids_shown(window) == ["A"], timeout=20_000)
    bar.edit.setText("[OH]a")  # SMARTS: an aromatic hydroxyl
    bar.apply_now()
    qtbot.waitUntil(lambda: ids_shown(window) == ["A", "C"], timeout=20_000)
    bar.edit.setText("C(")
    bar.apply_now()
    assert bar.edit.toolTip()  # not a valid query: says why, keeps the last result
    bar.set_column(None)  # back to text: the structure filter goes
    bar.edit.setText("B")
    bar.apply_now()
    assert ids_shown(window) == ["B"]


def test_search_bar_field_syntax_count_and_history(window, qtbot) -> None:
    from chemsuite.app.document import session

    bar = window.search
    bar.open()
    bar.edit.setText("pIC50>5.5 -ID:C")
    bar.apply_now()
    assert ids_shown(window) == ["B"]
    assert bar.count.isVisibleTo(bar) and bar.count.text() == "1 row"
    bar.edit.setText("pIC50>=6")
    bar.apply_now()
    assert ids_shown(window) == ["B", "C"] and bar.count.text() == "2 rows"
    table = window.workspace.table
    table.set_cells("pIC50", [int(table.oids[0])], ["9"])  # an edit updates the count
    qtbot.waitUntil(lambda: bar.count.text() == "3 rows", timeout=3000)
    bar.whole.setChecked(True)  # Whole Cell: the text as typed, no field syntax
    assert ids_shown(window) == []
    bar.whole.setChecked(False)

    bar.remember("pIC50>=6")
    bar.remember("ID:A")
    bar.remember("pIC50>=6")  # moves back to the top
    assert bar.history == ["pIC50>=6", "ID:A"]
    bar._fill_history()
    bar.history_menu.actions()[1].trigger()
    assert bar.edit.text() == "ID:A"

    extras = session.capture(window)
    bar.history = []
    session.restore(window, extras)
    assert window.search.history == ["pIC50>=6", "ID:A"]
    bar.dismiss()
    assert not bar.count.isVisibleTo(bar)


def test_find_replace_in_selected_cells_and_chosen_columns(window) -> None:
    from PySide6.QtCore import QItemSelectionModel

    table = window.workspace.table
    table.add_column("Notes", ["A-5", "B-6", "C-7"])
    view = window.view
    model = view.table_model
    notes = model.column_index("Notes")
    selection = view.selectionModel()
    flags = QItemSelectionModel.SelectionFlag
    selection.select(model.index(0, notes), flags.ClearAndSelect)
    selection.select(model.index(2, notes), flags.Select)
    dialog = window.show_find_replace()
    dialog.scope.setCurrentIndex(dialog.scope.findData("cells"))
    dialog.find_edit.setText(r"(\w)-(\d)")
    dialog.mode.setCurrentIndex(dialog.mode.findData("regex"))
    dialog.replace_edit.setText("$2$1")
    assert dialog.find_next()  # moves the cursor, so the selection changes...
    dialog.case.setChecked(True)  # ...but editing the query keeps the scope read at the start
    assert dialog.replace_everything() == 2
    assert [table.cell("Notes", int(o)) for o in table.oids] == ["5A", "B-6", "7C"]

    dialog.replace_edit.setText("$3")  # a group the expression does not have
    dialog.scope.setCurrentIndex(dialog.scope.findData("all"))
    assert dialog.replace_everything() == 0
    assert "invalid group reference 3" in dialog.status.text()
    assert table.cell("Notes", int(table.oids[1])) == "B-6"

    dialog.find_edit.setText("6")
    dialog.mode.setCurrentIndex(dialog.mode.findData("contains"))
    dialog.replace_edit.setText("9")
    dialog.set_chosen(["Notes"])
    assert dialog.columns_button.isVisibleTo(dialog)
    assert dialog.replace_everything() == 1  # Notes only; pIC50's 6 is left alone
    assert table.cell("Notes", int(table.oids[1])) == "B-9"
    assert table.cell("pIC50", int(table.oids[1])) == "6"
    dialog._fill_columns_menu()
    checked = [a.text() for a in dialog.columns_menu.actions() if a.isChecked()]
    assert checked == ["Notes"]


def test_search_options_are_glyph_toggles(window) -> None:
    from PySide6.QtWidgets import QToolButton

    bar = window.search
    buttons = bar.options_button.findChildren(QToolButton)
    assert [b.defaultAction() for b in buttons] == [bar.case, bar.whole, bar.regex]
    assert all(not b.icon().isNull() and b.isCheckable() for b in buttons)
    assert not bar.close_button.icon().isNull()
    buttons[2].click()  # Regular Expression
    assert bar.regex.isChecked() and bar.query().mode == "regex"
    buttons[1].click()  # Whole Cell turns the regular expression off
    assert bar.whole.isChecked() and not bar.regex.isChecked()
