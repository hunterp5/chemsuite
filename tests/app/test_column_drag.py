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

"""Dragging a column header: edge scrolling, and the view staying put after the drop."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QEvent, QPoint, QPointF, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication

from chemsuite.app.table.table_view import EDGE_MAX_SPEED, EDGE_MIN_SPEED, edge_speed

from .conftest import gui


def send(widget, kind, point: QPoint, buttons) -> None:
    button = Qt.MouseButton.NoButton if kind == QEvent.Type.MouseMove else Qt.MouseButton.LeftButton
    p = QPointF(point)
    QApplication.sendEvent(
        widget,
        QMouseEvent(kind, p, widget.mapToGlobal(p), button, buttons, Qt.KeyboardModifier(0)),
    )


def test_edge_speed_grows_toward_and_past_the_edge():
    width = 1000
    assert edge_speed(500, width) == 0
    assert edge_speed(100, width) == 0
    near, at, beyond = edge_speed(990, width), edge_speed(1000, width), edge_speed(1150, width)
    assert 0 < near < at < beyond <= EDGE_MAX_SPEED
    assert edge_speed(10, width) == -near
    assert edge_speed(5000, width) == EDGE_MAX_SPEED
    assert edge_speed(950, width) == pytest.approx(EDGE_MIN_SPEED, rel=0.05)
    assert edge_speed(5, 0) == 0


@gui
def test_dragging_a_column_to_the_edge_scrolls_and_drops_beyond_the_view(window, qtbot):
    table = window.workspace.table
    for i in range(30):
        table.add_column(f"Column number {i}", ["1", "2", "3"])
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    header = view.horizontalHeader()
    viewport = header.viewport()
    bar = view.horizontalScrollBar()
    assert bar.maximum() > 0
    start = QPoint(header.sectionViewportPosition(0) + 10, 5)
    edge = QPoint(viewport.width() + 40, 5)  # past the right edge: scrolls fast
    send(viewport, QEvent.Type.MouseButtonPress, start, Qt.MouseButton.LeftButton)
    for x in range(start.x(), edge.x(), 25):
        send(viewport, QEvent.Type.MouseMove, QPoint(x, 5), Qt.MouseButton.LeftButton)
    send(viewport, QEvent.Type.MouseMove, edge, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(lambda: bar.value() > 400, timeout=3000)
    edge.setX(viewport.width() - 2)
    send(viewport, QEvent.Type.MouseMove, edge, Qt.MouseButton.LeftButton)
    send(viewport, QEvent.Type.MouseButtonRelease, edge, Qt.MouseButton.NoButton)
    qtbot.waitUntil(lambda: table.column_names[0] != "ID", timeout=2000)
    qtbot.wait(20)
    position = table.column_names.index("ID")
    assert position > 5  # dropped where the table had scrolled to, far right of the start
    scrolled = bar.value()
    assert scrolled > 400  # the view did not jump back to the first column
    col = view.table_model.column_index("ID")
    left = header.sectionPosition(col)
    assert scrolled <= left < scrolled + view.viewport().width()  # the moved column is in view
    table.undo()
    assert table.column_names[0] == "ID"


@gui
def test_a_moved_column_stays_selected(window, qtbot):
    table = window.workspace.table
    view = window.view
    window.show()
    qtbot.waitExposed(window)
    header = view.horizontalHeader()
    point = QPoint(header.sectionViewportPosition(0) + 5, 5)
    send(header.viewport(), QEvent.Type.MouseButtonPress, point, Qt.MouseButton.LeftButton)
    send(header.viewport(), QEvent.Type.MouseButtonRelease, point, Qt.MouseButton.NoButton)
    assert view.whole_columns() == ["ID"]
    table.move_column("ID", 1)
    assert table.column_names[:2] == ["pIC50", "ID"]
    assert view.whole_columns() == ["ID"]
