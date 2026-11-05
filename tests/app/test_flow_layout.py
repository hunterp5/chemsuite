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

#: A panel docked beside the table must fit in about a quarter of a laptop screen.
MAX_PANEL_MIN_WIDTH = 300


def test_flow_layout_wraps_instead_of_widening(qtbot) -> None:
    from PySide6.QtWidgets import QPushButton, QWidget

    from chemsuite.app.flow_layout import FlowLayout

    host = QWidget()
    qtbot.addWidget(host)
    flow = FlowLayout(host)
    buttons = [QPushButton(f"Button number {i}") for i in range(6)]
    for button in buttons:
        flow.addWidget(button)
    widest = max(b.minimumSizeHint().width() for b in buttons)
    assert flow.minimumSize().width() == widest
    one_line = flow.heightForWidth(flow.sizeHint().width())
    assert flow.heightForWidth(widest) > 3 * one_line  # one button per line
    buttons[0].hide()
    assert flow.sizeHint().width() < sum(b.sizeHint().width() for b in buttons)


def test_plot_panels_can_be_docked_narrow(window) -> None:
    table = window.workspace.table
    table.add_column("A rather long column name for a measurement", ["1", "2", "3"])
    table.convert_column("A rather long column name for a measurement", "number")
    for kind in ("scatter", "histogram", "box"):
        panel = window.open_tool(kind)
        assert panel.minimumSizeHint().width() <= MAX_PANEL_MIN_WIDTH, kind


def test_statistics_panel_can_be_docked_narrow(window) -> None:
    panel = window.open_tool("statistics")
    for i in range(panel.tabs.count()):
        width = panel.tabs.widget(i).minimumSizeHint().width()
        assert width <= MAX_PANEL_MIN_WIDTH, panel.tabs.tabText(i)
