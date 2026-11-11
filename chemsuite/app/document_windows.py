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
"""Tool windows, and saving and reopening tool panels with the document.

A tool shows its work either in a dock beside the table (``ToolContext.add_panel``) or, for a
workspace of its own such as the Protein Viewer, in a :class:`ToolWindow`
(``ToolContext.add_window``): a top-level window with its own menu bar, titled after the
document, that cannot be docked. Either way it belongs to the document: it closes with it, and
when its widget is a ``ToolPanel`` with a ``tool_id`` it is saved with it and reopened by the
tool's ``restore`` (a window also gets its size and place back). :class:`DocumentPanels` owns
both lists.
"""

from __future__ import annotations

import base64
import logging
from collections.abc import Callable, Iterable
from typing import TYPE_CHECKING

import shiboken6
from PySide6.QtCore import QByteArray, Qt
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QDockWidget, QMainWindow, QMenu, QMenuBar, QWidget

from .panel_dock import DockHost, PanelDock
from .widgets.panel import ToolPanel

if TYPE_CHECKING:
    from .main_window import MainWindow
    from .tools import Tool

log = logging.getLogger(__name__)

GEOMETRY = "geometry"  # a saved window's size and place (base64 of saveGeometry)


class ToolWindow(QMainWindow):
    """A top-level window around a tool's widget, closed with its document.

    *name* gives the document's name (the title reads "Protein Viewer — Kinases") and
    *on_close* lets the widget go of the document before the window is deleted.
    """

    def __init__(
        self,
        widget: QWidget,
        title: str,
        *,
        name: Callable[[], str],
        menu_bar: QMenuBar | None = None,
        on_close: Callable[[], None] | None = None,
    ) -> None:
        super().__init__(None)  # its own window: not kept above the table, its own taskbar entry
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.panel = widget
        self._title = title
        self._name = name
        self._on_close = on_close
        self.setCentralWidget(widget)
        if menu_bar is not None:
            self.setMenuBar(menu_bar)
            menu_bar.show()
        self.retitle()

    def retitle(self) -> None:
        name = self._name()
        self.setWindowTitle(f"{self._title} — {name}" if name else self._title)

    def bring_forward(self) -> None:
        if self.isMinimized():
            self.showNormal()
        self.show()
        self.raise_()
        self.activateWindow()

    def closeEvent(self, event: QCloseEvent) -> None:
        if self._on_close is not None:
            self._on_close()
            self._on_close = None
        super().closeEvent(event)


class DocumentPanels(DockHost):
    """The tool panels of one document window: docks beside the table and windows of their
    own. The only owner of those lists; the window, its tools and its saved session ask it."""

    window: MainWindow

    def __init__(self, window: MainWindow) -> None:
        super().__init__(window)
        self.docks: list[PanelDock] = []
        self.windows: list[ToolWindow] = []
        #: dock_layout is putting Filters back in place (moving it fires signals that would
        #: call back into it)
        self.pinning = False

    def widgets(self) -> list[QWidget]:
        """Every open panel, docked or in a window."""
        return [dock.widget() for dock in self.docks] + [shown.panel for shown in self.windows]

    # -- docks ---------------------------------------------------------------------------------
    def new_dock(self, title: str) -> PanelDock:
        dock = PanelDock(title, self)
        self.docks.append(dock)
        dock.destroyed.connect(lambda _o=None, d=dock: self._forget_dock(d))
        return dock

    def _forget_dock(self, dock: PanelDock) -> None:
        if dock in self.docks:
            self.docks.remove(dock)

    def in_order(self) -> list[QDockWidget]:
        window = self.window
        return [window.filter_dock, window.search_dock, *self.docks]

    def fill_arrange(self, menu: QMenu, dock: QDockWidget) -> None:
        from .dock_layout import fill_arrange_menu

        fill_arrange_menu(menu, self.window, dock)

    # -- windows -------------------------------------------------------------------------------
    def add_window(
        self,
        widget: QWidget,
        title: str,
        menu_bar: QMenuBar | None,
        on_close: Callable[[], None] | None,
    ) -> ToolWindow:
        """Show *widget* in a new tool window of the document."""
        window = self.window
        shown = ToolWindow(widget, title, name=lambda: window.document_name, menu_bar=menu_bar,
                           on_close=on_close)  # fmt: skip
        window.windowTitleChanged.connect(shown.retitle)
        self.windows.append(shown)
        shown.destroyed.connect(lambda _o=None, w=shown: self._forget_window(w))
        shown.show()
        return shown

    def _forget_window(self, shown: ToolWindow) -> None:
        if shown in self.windows:
            self.windows.remove(shown)

    def find(self, tool_id: str) -> QWidget | None:
        """The open panel of *tool_id*, brought to the front, if there is one."""
        for dock in self.docks:
            widget = dock.widget()
            if isinstance(widget, ToolPanel) and widget.tool_id == tool_id:
                dock.show()
                dock.raise_()
                return widget
        for shown in self.windows:
            if isinstance(shown.panel, ToolPanel) and shown.panel.tool_id == tool_id:
                shown.bring_forward()
                return shown.panel
        return None

    def close_all(self) -> None:
        """Close every panel (the document is closing or being replaced)."""
        for dock in list(self.docks):
            dock.close()
        self.docks.clear()
        for shown in list(self.windows):
            if shiboken6.isValid(shown):
                shown.close()
        self.windows.clear()

    # -- saving and reopening ------------------------------------------------------------------
    def states(self) -> list[dict]:
        """Saved form of the open panels and tool windows that can be restored
        (``Tool.restore``).

        The saved docks are named ``Panel0``… in this order, which the saved arrangement
        refers to; a window's entry keeps its size and place.
        """
        from . import dock_layout

        states, docks = [], []
        for dock in self.docks:
            entry = _entry(dock.widget())
            if entry is not None:
                states.append(entry)
                docks.append(dock)
        dock_layout.name_panels(self.window, docks)
        for shown in self.windows:
            entry = _entry(shown.panel)
            if entry is not None:
                geometry = bytes(shown.saveGeometry().data())
                entry[GEOMETRY] = base64.b64encode(geometry).decode("ascii")
                states.append(entry)
        return states

    def restore(self, saved: list, tools: Iterable[Tool]) -> None:
        from .tools import ToolContext

        by_id = {}
        for tool in tools:
            by_id.update(dict.fromkeys((*tool.aliases, tool.id), tool))
        for i, entry in enumerate(saved):
            tool = by_id.get(entry.get("tool")) if isinstance(entry, dict) else None
            if tool is None or tool.restore is None or tool.missing():
                log.warning("cannot restore saved panel %r", entry)
                continue
            docks, windows = list(self.docks), list(self.windows)
            try:
                tool.restore(ToolContext(self.window), dict(entry.get("state") or {}))
            except (KeyError, ValueError, TypeError):  # a damaged or outdated entry
                log.exception("failed to restore panel %r", entry)
                continue
            added = [d for d in self.docks if d not in docks]
            if added:
                added[-1].setObjectName(f"Panel{i}")  # as named when saved
            opened = [w for w in self.windows if w not in windows]
            if opened and isinstance(entry.get(GEOMETRY), str):
                _place(opened[-1], entry[GEOMETRY])


def _entry(widget: QWidget) -> dict | None:
    """The saved form of a tool's widget, if it can be restored."""
    if isinstance(widget, ToolPanel) and widget.tool_id:
        return {"tool": widget.tool_id, "state": widget.state()}
    return None


def _place(shown: ToolWindow, text: str) -> None:
    try:
        geometry = base64.b64decode(text.encode("ascii"), validate=True)
    except ValueError:
        return
    shown.restoreGeometry(QByteArray(geometry))
