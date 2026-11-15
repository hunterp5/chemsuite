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
"""A selection made in a plot brings its rows into view in the table."""

from __future__ import annotations

import pytest
from PySide6.QtWidgets import QApplication

from tests.gui import gui

pytestmark = gui

ROWS = 500


@pytest.fixture
def long(qtbot, window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": [f"C{i}" for i in range(ROWS)]}, undoable=False)
    window._install(Workspace(table))
    window.show()
    qtbot.waitExposed(window)
    return window


def _visible(view) -> tuple[int, int]:
    rect = view.viewport().rect()
    return view.rowAt(rect.top()), view.rowAt(rect.bottom())


def test_a_plot_selection_scrolls_to_its_row(long) -> None:
    view, table = long.view, long.workspace.table
    long.workspace.selection.set([int(table.oids[400])], "plot:1")
    QApplication.processEvents()
    top, bottom = _visible(view)
    assert top <= 400 <= bottom


def test_rows_already_in_view_leave_the_table_still(long) -> None:
    view, table = long.view, long.workspace.table
    long.workspace.selection.set([int(table.oids[1]), int(table.oids[400])], "plot:1")
    QApplication.processEvents()
    assert view.verticalScrollBar().value() == 0


def test_only_new_rows_count(long) -> None:
    """Shift-adding a far row to a selection on screen brings the far row into view."""
    view, table, selection = long.view, long.workspace.table, long.workspace.selection
    selection.set([int(table.oids[1])], "plot:1")
    selection.add([int(table.oids[450])], "plot:1")
    QApplication.processEvents()
    top, bottom = _visible(view)
    assert top <= 450 <= bottom


def test_other_selections_do_not_scroll(long) -> None:
    view, table = long.view, long.workspace.table
    long.workspace.selection.set([int(table.oids[400])], "similarity")
    QApplication.processEvents()
    assert view.verticalScrollBar().value() == 0
