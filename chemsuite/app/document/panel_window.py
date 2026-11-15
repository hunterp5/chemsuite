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
"""A tool window whose own results dock beside it.

A workspace with a window of its own — the Protein Viewer — computes things there: what a
preparation decided, the poses a docking run produced. Those belong beside the structure they
describe, not beside the table in another window, so :class:`PanelWindow` is a
:class:`~chemsuite.app.shell.windows.ToolWindow` that can hold panel docks of its own, with
the same headers, tabs and Arrange menu as the panes beside the table.

:class:`WindowPanels` is its :class:`~chemsuite.app.document.panel_dock.DockHost`. Its Arrange
menu is the plain Qt one (beside, below, tabbed, which edge, float): the table's
``dock_layout`` knows about Filters, Search and saved arrangements, which a tool window has
none of.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QDockWidget, QMenu, QMenuBar, QWidget

from ..shell.windows import ToolWindow
from ..widgets.panel import ToolPanel
from .panel_dock import HIDE_DOCK_TABS, DockHost, PanelDock

#: a first panel gets this much of the window, or half of it, whichever is larger
PANEL_MIN_WIDTH = 320
#: which edge a panel can be moved to, as the Arrange menu offers them
AREAS = (
    (Qt.DockWidgetArea.RightDockWidgetArea, "Right"),
    (Qt.DockWidgetArea.LeftDockWidgetArea, "Left"),
    (Qt.DockWidgetArea.BottomDockWidgetArea, "Bottom"),
    (Qt.DockWidgetArea.TopDockWidgetArea, "Top"),
)


class WindowPanels(DockHost):
    """The panel docks of one tool window."""

    def __init__(self, window: PanelWindow) -> None:
        super().__init__(window)
        self.docks: list[PanelDock] = []

    def in_order(self) -> list[QDockWidget]:
        return list(self.docks)

    def fill_arrange(self, menu: QMenu, dock: QDockWidget) -> None:
        """Beside, below, tabbed with another panel, which edge, or floating."""
        window = self.window
        others = [d for d in self.docks if d is not dock and d.isVisible() and not d.isFloating()]
        for title, orientation in (
            ("Beside", Qt.Orientation.Horizontal),
            ("Below", Qt.Orientation.Vertical),
        ):
            sub = menu.addMenu(title)
            sub.setEnabled(bool(others))
            for other in others:
                sub.addAction(
                    other.windowTitle(),
                    lambda o=other, r=orientation: window.splitDockWidget(o, dock, r),
                )
        tabs = menu.addMenu("Tab With")
        tabs.setEnabled(bool(others))
        group = set(window.tabifiedDockWidgets(dock))
        for other in others:
            action = tabs.addAction(
                other.windowTitle(), lambda o=other: window.tabifyDockWidget(o, dock)
            )
            action.setCheckable(True)
            action.setChecked(other in group)
        menu.addSeparator()
        here = None if dock.isFloating() else window.dockWidgetArea(dock)
        for area, title in AREAS:
            action = menu.addAction(title, lambda a=area: self._move_to(dock, a))
            action.setCheckable(True)
            action.setChecked(area == here)
        menu.addSeparator()
        floating = menu.addAction("Float", lambda: dock.setFloating(not dock.isFloating()))
        floating.setCheckable(True)
        floating.setChecked(dock.isFloating())

    def _move_to(self, dock: QDockWidget, area: Qt.DockWidgetArea) -> None:
        dock.setFloating(False)
        self.window.addDockWidget(area, dock, Qt.Orientation.Horizontal)
        dock.show()
        dock.raise_()

    def forget(self, dock: PanelDock) -> None:
        if dock in self.docks:
            self.docks.remove(dock)


class PanelWindow(ToolWindow):
    """A tool window that can hold panel docks beside its own widget."""

    def __init__(
        self,
        widget: QWidget,
        title: str,
        *,
        name: Callable[[], str],
        menu_bar: QMenuBar | None = None,
        on_close: Callable[[], None] | None = None,
    ) -> None:
        super().__init__(widget, title, name=name, menu_bar=menu_bar, on_close=on_close)
        self.setStyleSheet(HIDE_DOCK_TABS)  # the tabs are in the panes' own headers
        self.panels = WindowPanels(self)

    @property
    def docks(self) -> list[PanelDock]:
        return self.panels.docks

    def add_panel(self, widget: QWidget, title: str) -> PanelDock:
        """Show *widget* in a dock beside this window's own. Later panels join the last one
        as tabs, as they do beside the table."""
        previous = [d for d in self.docks if d.isVisible() and not d.isFloating()]
        dock = PanelDock(title, self.panels)
        self.panels.docks.append(dock)
        dock.destroyed.connect(lambda _o=None, d=dock: self.panels.forget(d))
        dock.setObjectName(f"Panel{id(dock)}")
        dock.setWidget(widget)
        if isinstance(widget, ToolPanel):
            widget.titleChanged.connect(dock.setWindowTitle)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock, Qt.Orientation.Horizontal)
        if previous:
            self.tabifyDockWidget(previous[-1], dock)
        else:
            self.resizeDocks(
                [dock], [max(PANEL_MIN_WIDTH, self.width() // 3)], Qt.Orientation.Horizontal
            )
        dock.show()
        dock.raise_()
        self.panels.schedule_header_refresh()
        return dock

    def find_panel(self, tool_id: str) -> QWidget | None:
        """The open panel of *tool_id* in this window, brought to the front."""
        for dock in self.docks:
            widget = dock.widget()
            if isinstance(widget, ToolPanel) and widget.tool_id == tool_id:
                self.bring_forward()
                dock.show()
                dock.raise_()
                return widget
        return None

    def panel_widgets(self) -> list[QWidget]:
        """The widgets of this window's docks (not its own central widget)."""
        return [dock.widget() for dock in self.docks]

    def closeEvent(self, event: QCloseEvent) -> None:
        for dock in list(self.docks):
            dock.close()
        super().closeEvent(event)
