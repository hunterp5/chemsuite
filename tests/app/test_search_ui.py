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
