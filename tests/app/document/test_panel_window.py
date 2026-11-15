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
"""A tool window holding panel docks: closing one leaves the others' headers sound."""

from __future__ import annotations

from PySide6.QtCore import QCoreApplication, QEvent
from PySide6.QtWidgets import QWidget

from chemsuite.app.document.panel_dock import HeaderTabs, PanelHeader
from chemsuite.app.document.panel_window import PanelWindow
from chemsuite.app.widgets.panel import ToolPanel
from tests.gui import gui

pytestmark = gui


def test_a_closed_tab_leaves_its_neighbours_header_at_once(qtbot) -> None:
    """The headers learn of a close on a zero-delay timer; the closed dock can be deleted, and
    a neighbour's header resized, before that timer fires (seen on macOS). Its tabs must not
    still hold the deleted dock then."""
    window = PanelWindow(QWidget(), "Viewer", name=lambda: "Viewer")
    qtbot.addWidget(window)
    window.resize(900, 600)
    window.show()
    first = window.add_panel(ToolPanel(), "Report")
    second = window.add_panel(ToolPanel(), "Inspect")
    qtbot.waitUntil(lambda: second in first.titleBarWidget().tabs.docks, timeout=5_000)
    second.close()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    tabs = first.titleBarWidget().tabs
    assert isinstance(first.titleBarWidget(), PanelHeader) and isinstance(tabs, HeaderTabs)
    assert second not in tabs.docks
    tabs.resize(tabs.width() + 10, tabs.height())  # lays the tabs out again: no dead dock
    assert tabs.titles() == [] or tabs.titles() == ["Report"]
