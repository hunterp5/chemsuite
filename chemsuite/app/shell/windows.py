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

"""Top-level windows other than a document's main window."""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QMainWindow, QMenuBar, QWidget


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
