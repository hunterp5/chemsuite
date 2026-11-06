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
"""Docks for tool panels, with a header that is the top of the panel's card.

The header holds, left to right: the panel's title, the panel's own toolbar (a plot's pan, box,
lasso…), then settings, undock and close. A panel opts into each part by what it has:

* ``toolbar``: a ``QToolBar`` moved into the header;
* ``options``: a popup with ``popup(anchor)``, opened by the settings button;
* ``set_custom_title(text | None)``: makes the title editable (double-click it; an empty title
  goes back to the automatic one);
* ``attach_header(bool)``: told when its top edge now meets the header, to square its corners.

Mouse presses on the header's empty space and on the title (when not editing) fall through to
the dock, so dragging the header still moves the panel.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QEvent, QSize, Qt
from PySide6.QtGui import QFocusEvent, QKeyEvent, QMouseEvent
from PySide6.QtWidgets import (
    QDockWidget,
    QFrame,
    QHBoxLayout,
    QLineEdit,
    QMenu,
    QSizePolicy,
    QToolButton,
    QWidget,
)

from .glyphs import glyph

#: panes sit flush against each other, as in native macOS split views
INSET = 0


class TitleEdit(QLineEdit):
    """A title that looks like a label until double-clicked, then edits in place.

    Enter keeps the new title, Escape or clicking away with nothing typed puts it back.
    """

    def __init__(self, on_rename: Callable[[str], None] | None = None) -> None:
        super().__init__()
        self.on_rename = on_rename
        self._before = ""
        self.setFrame(False)
        self.setReadOnly(True)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setMinimumWidth(40)  # a long title is cut off, not a wide minimum
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.setCursorPosition(0)
        self.editingFinished.connect(self._finish)
        self._style()

    @property
    def editable(self) -> bool:
        return self.on_rename is not None

    @property
    def editing(self) -> bool:
        return not self.isReadOnly()

    def show_title(self, text: str) -> None:
        if not self.editing:
            self.setText(text)
            self.setCursorPosition(0)  # show the start of a long title
            self.setToolTip(text + ("\nDouble-click to rename" if self.editable else ""))

    def start_editing(self) -> None:
        if not self.editable or self.editing:
            return
        self._before = self.text()
        self.setReadOnly(False)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._style()
        self.setFocus(Qt.FocusReason.MouseFocusReason)
        self.selectAll()

    def _stop(self) -> None:
        self.setReadOnly(True)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.deselect()
        self._style()

    def _finish(self) -> None:
        if not self.editing:
            return
        text = self.text().strip()
        self._stop()
        if self.on_rename is not None and text != self._before:
            self.on_rename(text)  # the panel answers with show_title (automatic if empty)
        else:
            self.show_title(self._before)

    def cancel(self) -> None:
        if self.editing:
            self.setText(self._before)
            self._stop()
            self.show_title(self._before)

    def _style(self) -> None:
        weight = "font-weight: 600;"
        if self.editing:
            self.setStyleSheet(
                f"QLineEdit {{ {weight} background: palette(base);"
                " border: 1px solid palette(highlight); border-radius: 3px; padding: 0 3px; }"
            )
        else:
            self.setStyleSheet(
                f"QLineEdit {{ {weight} background: transparent; border: none; padding: 0 4px; }}"
            )

    # -- events: read-only presses go to the dock (drag); a double-click edits ----------------
    def mousePressEvent(self, event: QMouseEvent) -> None:
        if self.editing:
            super().mousePressEvent(event)
        else:
            event.ignore()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        if self.editing:
            super().mouseMoveEvent(event)
        else:
            event.ignore()

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if self.editing:
            super().mouseReleaseEvent(event)
        else:
            event.ignore()

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:
        if self.editable and not self.editing:
            self.start_editing()
            event.accept()
        elif self.editing:
            super().mouseDoubleClickEvent(event)
        else:
            event.ignore()  # the dock floats or docks, as with Qt's own title bar

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if self.editing and event.key() == Qt.Key.Key_Escape:
            self.cancel()
            return
        super().keyPressEvent(event)

    def focusOutEvent(self, event: QFocusEvent) -> None:
        super().focusOutEvent(event)
        self._finish()


class PanelHeader(QWidget):
    """The dock's title bar: a flat strip over the pane, divided from it by a hairline."""

    def __init__(self, dock: QDockWidget) -> None:
        super().__init__(dock)
        self.dock = dock
        self.strip = QFrame()
        self.strip.setObjectName("PanelHeader")
        self.title = TitleEdit()
        self.title.show_title(dock.windowTitle())

        def button(tip: str) -> QToolButton:
            b = QToolButton()
            b.setAutoRaise(True)
            b.setToolTip(tip)
            b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            b.setIconSize(QSize(14, 14))
            return b

        self.settings = button("Settings")
        self.settings.setVisible(False)
        self.arrange_button = button("Arrange: beside or below another pane, tabs, float")
        self.arrange_menu = QMenu(self.arrange_button)
        self.arrange_button.setMenu(self.arrange_menu)
        self.arrange_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.arrange_button.setStyleSheet("QToolButton::menu-indicator { image: none; }")
        self.arrange_menu.aboutToShow.connect(self._fill_arrange)
        self.float_button = button("Undock (or dock again)")
        self.close_button = button("Close")
        self._glyphs = {
            self.settings: "gear",
            self.arrange_button: "arrange",
            self.float_button: "undock",
            self.close_button: "close",
        }
        self.float_button.clicked.connect(lambda: dock.setFloating(not dock.isFloating()))
        self.close_button.clicked.connect(dock.close)

        self.row = QHBoxLayout(self.strip)
        self.row.setContentsMargins(8, 3, 4, 3)
        self.strip.setFixedHeight(32)  # every pane's title bar, with or without a toolbar
        self.row.setSpacing(2)
        self.row.addWidget(self.title, 1)
        self.tools_at = self.row.count()  # where a panel's toolbar goes
        for b in (self.settings, self.arrange_button, self.float_button, self.close_button):
            self.row.addWidget(b)
        outer = QHBoxLayout(self)
        outer.setContentsMargins(INSET, INSET, INSET, 0)
        outer.addWidget(self.strip)
        self.toolbar: QWidget | None = None
        dock.windowTitleChanged.connect(self.title.show_title)
        self._restyle()

    def _fill_arrange(self) -> None:
        from .dock_layout import fill_arrange_menu

        window = self.dock.parentWidget()
        if window is not None and hasattr(window, "_panels"):
            fill_arrange_menu(self.arrange_menu, window, self.dock)

    def add_tools(self, toolbar: QWidget) -> None:
        self.toolbar = toolbar
        self.row.insertWidget(self.tools_at, toolbar)

    def _restyle(self) -> None:
        palette = self.palette()
        ink = palette.text().color()
        for b, name in self._glyphs.items():
            b.setIcon(glyph(name, ink))
        line = palette.mid().color().name()
        window = palette.window().color().name()
        self.strip.setStyleSheet(
            f"QFrame#PanelHeader {{ background: {window}; border: none;"
            f" border-bottom: 1px solid {line}; }}"
            " QToolButton { border: none; border-radius: 4px; padding: 2px; }"
            " QToolButton:hover { background: palette(midlight); }"
        )

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "strip"):
            self._restyle()


class PanelDock(QDockWidget):
    """A dock holding a tool's panel; closing it closes the panel for good."""

    def __init__(self, title: str, parent: QWidget | None = None) -> None:
        super().__init__(title, parent)
        self.title_bar = PanelHeader(self)
        self.setTitleBarWidget(self.title_bar)

    def setWidget(self, widget: QWidget) -> None:
        super().setWidget(widget)
        bar = self.title_bar
        options = getattr(widget, "options", None)
        bar.settings.setVisible(options is not None)
        if options is not None:
            bar.settings.clicked.connect(lambda: options.popup(bar.settings))
        toolbar = getattr(widget, "toolbar", None)
        if toolbar is not None:
            bar.add_tools(toolbar)
        if hasattr(widget, "set_custom_title"):
            bar.title.on_rename = widget.set_custom_title
            bar.title.show_title(self.windowTitle())
        if hasattr(widget, "attach_header"):
            widget.attach_header(True)

    def closeEvent(self, event) -> None:
        widget = self.widget()
        if widget is not None and hasattr(widget, "close_panel"):
            widget.close_panel()
        super().closeEvent(event)
        self.deleteLater()
