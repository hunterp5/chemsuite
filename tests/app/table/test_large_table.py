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
"""Large tables: painting and resetting never ask the model about every row."""

from __future__ import annotations

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from tests.gui import gui

pytestmark = gui

ROWS = 4000
#: calls to the table model counted here (PySide looks an override up once per model, so the
#: counting methods are put on the class before the model is made)
CALLS = {"flags": 0, "row labels": 0}


@pytest.fixture
def big(qtbot, window, monkeypatch):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.app.table.table_model import TableModel
    from chemsuite.core.table import Table

    flags, header_data = TableModel.flags, TableModel.headerData

    def counted_flags(self, index):
        CALLS["flags"] += 1
        return flags(self, index)

    def counted_header_data(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Vertical:
            CALLS["row labels"] += 1
        return header_data(self, section, orientation, role)

    monkeypatch.setattr(TableModel, "flags", counted_flags)
    monkeypatch.setattr(TableModel, "headerData", counted_header_data)
    table = Table()
    table.add_rows(
        {"ID": [f"C{i}" for i in range(ROWS)], "pIC50": [str(i % 97) for i in range(ROWS)]},
        undoable=False,
    )
    table.convert_column("pIC50", "number")
    window._install(Workspace(table))
    window.show()
    qtbot.waitExposed(window)
    return window


@pytest.mark.parametrize("share", [1, 2])
def test_column_header_paint_does_not_ask_every_row(big, share) -> None:
    """With all or the first half of the rows selected, painting the column header (on every
    hover) must not check each row's flags, as Qt's isColumnSelected does for highlighting."""
    view, table = big.view, big.workspace.table
    big.workspace.selection.set(table.oids[: ROWS // share].tolist(), "test")
    QApplication.processEvents()  # the table's own repaint asks about its visible cells
    CALLS["flags"] = 0
    view.horizontalHeader().viewport().repaint()
    assert CALLS["flags"] < 100
    assert view.horizontalHeader().property("highlighted_sections") == [0, 1]


def test_a_reset_does_not_measure_every_row_number(big) -> None:
    """Qt sizes the row header after every model reset by reading hundreds of row labels;
    its width is fixed, so it should not."""
    QApplication.processEvents()
    CALLS["row labels"] = 0
    big.search.edit.setText("C1")
    big.search.apply_now()
    QApplication.processEvents()  # the layout after a reset is deferred
    assert CALLS["row labels"] < 500  # those painted; measuring alone reads over 800 more


@pytest.mark.parametrize("how", ["select all", "column"])
def test_header_paint_after_select_all_or_a_column_does_not_ask_every_row(qtbot, big, how):
    """Select All (one range over everything) or a column picked on its header: Qt asks the
    header's own selection model whether each shown column is selected, row by row, so it must
    not hold those columns while nothing is being Ctrl+clicked."""
    from tests.app.table.test_table_interactions import click_header

    view = big.view
    if how == "select all":
        view.selectAll()
    else:
        click_header(qtbot, view, 1)
    QApplication.processEvents()
    CALLS["flags"] = 0
    view.horizontalHeader().viewport().repaint()
    assert CALLS["flags"] < 100
    assert len(view.horizontalHeader().selectionModel().selection()) == 0


def test_selected_column_names_do_not_ask_every_row(big) -> None:
    view = big.view
    view.selectAll()
    CALLS["flags"] = 0
    assert view.selected_column_names() == ["ID", "pIC50"]
    assert CALLS["flags"] < 100


def test_a_scattered_selection_is_replaced_without_comparing_ranges(big, monkeypatch) -> None:
    """Swapping one selection of scattered rows for another with ClearAndSelect makes Qt
    compare the old ranges with the new pairwise (26 s at 50,000 rows): clear, then select."""
    from PySide6.QtCore import QItemSelectionModel

    view, table, selection = big.view, big.workspace.table, big.workspace.selection
    model = view.selectionModel()
    flags_used = []
    select = model.select
    monkeypatch.setattr(
        model, "select", lambda what, flags: (flags_used.append(flags), select(what, flags))[1]
    )
    selection.set(table.oids[::2].tolist(), "test")
    selection.set(table.oids[1::2].tolist(), "test")
    clear_and_select = QItemSelectionModel.SelectionFlag.ClearAndSelect
    assert flags_used and all(f != clear_and_select for f in flags_used)
    assert len(view.selected_oids()) == ROWS // 2
    assert view.selected_oids().tolist() == table.oids[1::2].tolist()
