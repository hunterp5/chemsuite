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

"""Structure row height: dragging a row divider, Ctrl+wheel and the Layout menu's zoom, and
the height saved with the document."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QApplication

from chemsuite.core.settings import DEFAULT_ROW_HEIGHT, MAX_ROW_HEIGHT
from tests.gui import gui

pytestmark = gui


def _install(window, columns: dict[str, list[str]], structure: str | None = "SMILES"):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(columns, undoable=False)
    if structure:
        table.convert_column(structure, "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


@pytest.fixture
def structures(window):
    return _install(window, {"ID": ["A", "B", "C"], "SMILES": ["CCO", "c1ccccc1", "CC(=O)O"]})


def _divider(view, row: int) -> QPoint:
    header = view.verticalHeader()
    y = header.sectionViewportPosition(row) + header.sectionSize(row)
    return QPoint(header.width() // 2, y)


def _drag(qtbot, view, row: int, dy: int) -> None:
    viewport = view.verticalHeader().viewport()
    start = _divider(view, row)
    qtbot.mousePress(viewport, Qt.MouseButton.LeftButton, pos=start)
    qtbot.mouseMove(viewport, start + QPoint(0, dy))
    qtbot.mouseRelease(viewport, Qt.MouseButton.LeftButton, pos=start + QPoint(0, dy))


def _wheel(view, notches: int, modifiers=Qt.KeyboardModifier.ControlModifier) -> None:
    pos = QPointF(view.viewport().rect().center())
    event = QWheelEvent(
        pos,
        QPointF(view.viewport().mapToGlobal(pos.toPoint())),
        QPoint(),
        QPoint(0, 120 * notches),
        Qt.MouseButton.NoButton,
        modifiers,
        Qt.ScrollPhase.NoScrollPhase,
        False,
    )
    QApplication.sendEvent(view.viewport(), event)


def test_dragging_a_divider_sizes_every_row(structures, qtbot) -> None:
    view = structures.view
    assert view.rowHeight(0) == DEFAULT_ROW_HEIGHT
    _drag(qtbot, view, 0, 50)
    assert view.row_sizer.height == DEFAULT_ROW_HEIGHT + 50
    assert {view.rowHeight(r) for r in range(3)} == {DEFAULT_ROW_HEIGHT + 50}
    assert view.selectionModel().selectedRows() == []  # the drag selected nothing
    # the second divider moves twice as far per pixel of height, so it stays under the pointer
    _drag(qtbot, view, 1, 40)
    assert view.row_sizer.height == DEFAULT_ROW_HEIGHT + 70
    _drag(qtbot, view, 0, 10_000)
    assert view.row_sizer.height == MAX_ROW_HEIGHT


def test_double_click_on_a_divider_puts_the_settings_height_back(structures, qtbot) -> None:
    view = structures.view
    view.row_sizer.set_height(250)
    qtbot.mouseDClick(
        view.verticalHeader().viewport(), Qt.MouseButton.LeftButton, pos=_divider(view, 0)
    )
    assert view.rowHeight(0) == DEFAULT_ROW_HEIGHT
    assert view.row_sizer.saved() is None


def test_ctrl_wheel_zooms_rows_and_drawing_columns_together(structures) -> None:
    view = structures.view
    id_width, drawing_width = view.columnWidth(0), view.columnWidth(1)
    scroll = view.verticalScrollBar().value()
    _wheel(view, 2)
    height = view.row_sizer.height
    assert height == round(DEFAULT_ROW_HEIGHT * 1.15**2)
    assert view.rowHeight(2) == height
    assert view.columnWidth(1) == round(drawing_width * height / DEFAULT_ROW_HEIGHT)
    assert view.columnWidth(0) == id_width  # text columns keep their width
    _wheel(view, -2)
    assert view.row_sizer.height == DEFAULT_ROW_HEIGHT
    assert view.columnWidth(1) == drawing_width  # back exactly
    _wheel(view, 3, Qt.KeyboardModifier.NoModifier)  # without Ctrl the table scrolls
    assert view.row_sizer.height == DEFAULT_ROW_HEIGHT
    assert view.verticalScrollBar().value() >= scroll


def test_layout_menu_commands_zoom(structures) -> None:
    view = structures.view
    commands = structures._commands
    commands["view.larger_structures"].trigger()
    assert view.row_sizer.height == round(DEFAULT_ROW_HEIGHT * 1.15)
    commands["view.smaller_structures"].trigger()
    commands["view.smaller_structures"].trigger()
    assert view.row_sizer.height == round(DEFAULT_ROW_HEIGHT * 1.15**-1)
    commands["view.reset_structure_size"].trigger()
    assert view.row_sizer.height == DEFAULT_ROW_HEIGHT and view.row_sizer.saved() is None


def test_text_only_tables_do_not_resize(window, qtbot) -> None:
    view = _install(window, {"ID": ["A", "B"], "Note": ["x", "y"]}, structure=None).view
    line = view.rowHeight(0)
    _drag(qtbot, view, 0, 50)
    _wheel(view, 2)
    assert view.rowHeight(0) == line and view.row_sizer.saved() is None


def test_the_height_is_saved_with_the_document(structures) -> None:
    from chemsuite.app.document import dock_layout

    view = structures.view
    assert "row_height" not in dock_layout.capture(structures)  # the Settings height
    view.row_sizer.set_height(260)
    saved = dock_layout.capture(structures)
    assert saved["row_height"] == 260
    view.row_sizer.reset()
    dock_layout.restore(structures, saved)
    assert view.rowHeight(0) == 260
    dock_layout.restore(structures, {"row_height": "tall"})  # damaged: ignored
    assert view.rowHeight(0) == 260


def test_a_new_settings_height_replaces_the_documents(structures) -> None:
    from chemsuite.app.shell.settings_hub import settings_hub

    view = structures.view
    view.row_sizer.set_height(260)
    settings_hub().update(lambda s: setattr(s.appearance, "status_bar", False))
    assert view.rowHeight(0) == 260  # another setting leaves it alone
    settings_hub().update(lambda s: setattr(s.appearance, "structure_row_height", 150))
    assert view.rowHeight(0) == 150 and view.row_sizer.saved() is None


def test_frozen_rows_follow_and_their_dividers_drag(structures, qtbot) -> None:
    view = structures.view
    view.set_frozen_rows(1)
    frozen = view._frozen_rows
    view.row_sizer.set_height(200)
    assert frozen.rowHeight(0) == 200 and frozen.frozen_height() == 200
    _drag(qtbot, frozen, 0, 20)
    assert view.rowHeight(1) == 220 and frozen.rowHeight(0) == 220
