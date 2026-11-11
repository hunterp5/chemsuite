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
lasso…), then settings, undock and close. A panel (a ``ToolPanel``) opts into each part:

* ``toolbar``: a ``QToolBar`` moved into the header;
* ``options``: an ``OptionsPopup``, opened by the settings button;
* ``renamable`` with ``set_custom_title(text | None)``: makes the title editable (double-click
  it; an empty title goes back to the automatic one);
* ``attach_header(bool)``: told when its top edge now meets the header, to square its corners.

Where the headers' tabs come from and what Arrange offers is the window's :class:`DockHost`.

Mouse presses on the header's empty space and on the title (when not editing) fall through to
the dock, so dragging the header still moves the panel.

Panes tabbed together show their tabs in the header in place of the title (Qt's own tab bar
below the panes is hidden, see :data:`HIDE_DOCK_TABS`): the current pane's tab is raised, a
click on another tab brings that pane forward, and double-clicking the current tab renames it.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QEvent, QPointF, QRectF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QFocusEvent,
    QFont,
    QFontMetrics,
    QKeyEvent,
    QMouseEvent,
    QPainter,
    QPaintEvent,
)
from PySide6.QtWidgets import (
    QDockWidget,
    QFrame,
    QHBoxLayout,
    QLineEdit,
    QMainWindow,
    QMenu,
    QSizePolicy,
    QToolButton,
    QToolTip,
    QWidget,
)

from .glyphs import glyph
from .widgets.panel import ToolPanel

#: panes sit flush against each other, as in native macOS split views
INSET = 0
#: the height of every pane's title bar, with or without a toolbar or tabs
HEADER_HEIGHT = 26
#: Qt's tab bar under tabbed panes, hidden because the tabs are in the panes' headers (only
#: the main window's own tab bars match, not tab widgets inside panes)
HIDE_DOCK_TABS = (
    "QMainWindow > QTabBar { max-height: 0px; }"
    " QMainWindow > QTabBar::tab { height: 0px; min-height: 0px; max-height: 0px;"
    " padding: 0px; margin: 0px; border: none; }"
)


class TitleEdit(QLineEdit):
    """A title that looks like a label until double-clicked, then edits in place.

    Enter keeps the new title, Escape or clicking away with nothing typed puts it back.
    """

    #: editing ended, kept or cancelled
    stopped = Signal()

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
        self.stopped.emit()

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


class HeaderTabs(QWidget):
    """The tabs of panes sharing one space, drawn in the current pane's header.

    The current tab is a raised segment; the others are plain text that darkens on hover.
    Pressing the current tab falls through to the dock (dragging it moves the pane), a click on
    another tab raises that pane, and double-clicking the current tab asks to rename it.
    """

    rename_requested = Signal()
    PAD = 9  # text inset within a tab
    GAP = 2
    MIN_TAB = 44

    def __init__(self, dock: QDockWidget) -> None:
        super().__init__()
        self.dock = dock
        self.docks: list[QDockWidget] = []
        self._rects: list[QRectF] = []
        self._hover = -1
        #: double-clicking the current tab renames the pane
        self.renamable = False
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.setMinimumWidth(self.MIN_TAB)

    def set_docks(self, docks: list[QDockWidget]) -> None:
        self.docks = list(docks)
        self._hover = -1
        self._layout_tabs()
        self.update()

    def titles(self) -> list[str]:
        return [d.windowTitle() for d in self.docks]

    def _font(self, current: bool) -> QFont:
        font = self.font()
        if current:
            font.setWeight(font.Weight.DemiBold)
        return font

    def _layout_tabs(self) -> None:
        """Tabs take their text's width; when that does not fit, they share the width evenly."""
        self._rects = []
        if not self.docks:
            return
        height = self.height()
        widths = []
        for d in self.docks:
            metrics = self.fontMetrics() if d is not self.dock else _metrics(self._font(True))
            widths.append(metrics.horizontalAdvance(d.windowTitle()) + 2 * self.PAD + 1)
        room = self.width() - self.GAP * (len(widths) - 1)
        if sum(widths) > room:
            # Shrink the widest tabs first, as far as an even share
            share = max(self.MIN_TAB, room / len(widths))
            widths = _fit(widths, max(room, share * len(widths)), self.MIN_TAB)
        x = 0.0
        tab_h = min(height - 4, 20)
        y = (height - tab_h) / 2
        for w in widths:
            self._rects.append(QRectF(x, y, w, tab_h))
            x += w + self.GAP

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._layout_tabs()

    def tab_at(self, pos) -> int:
        for i, rect in enumerate(self._rects):
            if rect.contains(pos):
                return i
        return -1

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        palette = self.palette()
        ink = palette.text().color()
        for i, (dock, rect) in enumerate(zip(self.docks, self._rects, strict=False)):
            if rect.left() >= self.width():
                break
            current = dock is self.dock
            box = rect.adjusted(0.5, 0.5, -0.5, -0.5)
            if current:
                painter.setPen(QColor(palette.mid().color()))
                painter.setBrush(palette.base())
                painter.drawRoundedRect(box, 5, 5)
            elif i == self._hover:
                hover = QColor(ink)
                hover.setAlphaF(0.07)
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(hover)
                painter.drawRoundedRect(box, 5, 5)
            font = self._font(current)
            painter.setFont(font)
            colour = QColor(ink)
            if not current:
                colour.setAlphaF(0.85 if i == self._hover else 0.6)
            painter.setPen(colour)
            text_rect = rect.adjusted(self.PAD, 0, -self.PAD, 0)
            text = _metrics(font).elidedText(
                dock.windowTitle(), Qt.TextElideMode.ElideRight, int(text_rect.width())
            )
            painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, text)
        painter.end()

    # -- events -------------------------------------------------------------------------------
    def event(self, event: QEvent) -> bool:
        if event.type() == QEvent.Type.ToolTip:
            i = self.tab_at(QPointF(event.pos()))
            if i < 0:
                QToolTip.hideText()
                return True
            tip = self.docks[i].windowTitle()
            if self.docks[i] is not self.dock:
                tip += "\nClick to show"
            elif self.renamable:
                tip += "\nDouble-click to rename"
            QToolTip.showText(event.globalPos(), tip, self)
            return True
        return super().event(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        i = self.tab_at(event.position())
        if i != self._hover:
            self._hover = i
            self.update()
        event.ignore()  # a drag that began on the current tab moves the dock

    def leaveEvent(self, event) -> None:
        self._hover = -1
        self.update()
        super().leaveEvent(event)

    def mousePressEvent(self, event: QMouseEvent) -> None:
        i = self.tab_at(event.position())
        if (
            event.button() == Qt.MouseButton.LeftButton
            and i >= 0
            and self.docks[i] is not self.dock
        ):
            self.docks[i].raise_()
            event.accept()
            return
        event.ignore()  # the current tab or empty space: the dock handles it (drag)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        event.ignore()

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:
        i = self.tab_at(event.position())
        if i >= 0 and self.docks[i] is self.dock and self.renamable:
            self.rename_requested.emit()
            event.accept()
            return
        event.ignore()


def _metrics(font: QFont) -> QFontMetrics:
    return QFontMetrics(font)


def _fit(widths: list[float], room: float, floor: float) -> list[float]:
    """Cap the widest of *widths* so they sum to *room* (never below *floor*)."""
    if sum(widths) <= room:
        return widths
    lo, hi = floor, max(widths)
    for _ in range(30):  # the cap that makes them fit, by bisection
        cap = (lo + hi) / 2
        if sum(min(w, cap) for w in widths) > room:
            hi = cap
        else:
            lo = cap
    return [min(w, lo) for w in widths]


class DockHost:
    """What the panes' headers ask of the window they sit in: the panes in tab order, the
    Arrange menu, and a refresh once Qt has finished regrouping them. The document window's
    panels (``document_windows.DocumentPanels``) are one."""

    def __init__(self, window: QMainWindow) -> None:
        self.window = window
        self._refresh_pending = False

    def in_order(self) -> list[QDockWidget]:
        """Every pane in a fixed order: the order of their tabs."""
        return []

    def fill_arrange(self, menu: QMenu, dock: QDockWidget) -> None:
        """Fill a pane's Arrange menu."""

    def refresh_headers(self) -> None:
        """Bring every pane header's tabs up to date with how the panes are grouped."""
        self._refresh_pending = False
        for dock in self.window.findChildren(QDockWidget):
            header = dock.titleBarWidget()
            if isinstance(header, PanelHeader):
                header.update_tabs()

    def schedule_header_refresh(self) -> None:
        """Refresh the headers once control returns to the event loop (Qt regroups docks in
        several steps, and the grouping is only final afterwards)."""
        if self._refresh_pending:
            return
        self._refresh_pending = True
        QTimer.singleShot(0, self.window, self.refresh_headers)


class PanelHeader(QWidget):
    """The dock's title bar: a flat strip over the pane, divided from it by a hairline."""

    def __init__(self, dock: QDockWidget, host: DockHost) -> None:
        self._built = False  # a palette change can arrive while it is being built
        super().__init__(dock)
        self.dock = dock
        self.host = host
        self.strip = QFrame()
        self.strip.setObjectName("PanelHeader")
        self.title = TitleEdit()
        self.title.show_title(dock.windowTitle())

        def button(tip: str) -> QToolButton:
            b = QToolButton()
            b.setAutoRaise(True)
            b.setToolTip(tip)
            b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            b.setIconSize(QSize(13, 13))
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
        # a pinned pane (Filters) is neither floated nor arranged; an always-floating one
        # (Advanced Search) cannot be docked
        pinned = not dock.features() & QDockWidget.DockWidgetFeature.DockWidgetFloatable
        undockable = dock.allowedAreas() == Qt.DockWidgetArea.NoDockWidgetArea
        self.float_button.setVisible(not pinned and not undockable)
        self.arrange_button.setVisible(not pinned and not undockable)
        self.close_button.clicked.connect(dock.close)

        self.tabs = HeaderTabs(dock)
        self.tabs.setVisible(False)
        self.tabs.rename_requested.connect(self._rename_from_tab)
        self.title.stopped.connect(self.update_tabs)

        self.row = QHBoxLayout(self.strip)
        self.row.setContentsMargins(6, 2, 4, 2)
        self.strip.setFixedHeight(HEADER_HEIGHT)
        self.row.setSpacing(2)
        self.row.addWidget(self.title, 1)
        self.row.addWidget(self.tabs, 1)
        self.tools_at = self.row.count()  # where a panel's toolbar goes
        for b in (self.settings, self.arrange_button, self.float_button, self.close_button):
            self.row.addWidget(b)
        outer = QHBoxLayout(self)
        outer.setContentsMargins(INSET, INSET, INSET, 0)
        outer.addWidget(self.strip)
        self.toolbar: QWidget | None = None
        dock.windowTitleChanged.connect(self.title.show_title)
        for signal in (
            dock.windowTitleChanged,
            dock.visibilityChanged,
            dock.topLevelChanged,
            dock.dockLocationChanged,
        ):
            signal.connect(lambda *_: host.schedule_header_refresh())
        self._built = True
        self._restyle()

    # -- tabs -----------------------------------------------------------------------------------
    def tab_docks(self) -> list[QDockWidget]:
        """The panes tabbed with this one, in tab order (empty when it is not tabbed)."""
        if self.dock.isFloating():
            return []
        from .dock_layout import is_open

        group = {d for d in self.host.window.tabifiedDockWidgets(self.dock) if is_open(d)}
        if not group:
            return []
        group.add(self.dock)
        ordered = [d for d in self.host.in_order() if d in group]
        return ordered + [d for d in group if d not in ordered]

    def update_tabs(self) -> None:
        docks = self.tab_docks()
        tabbed = bool(docks) and not self.title.editing
        self.tabs.renamable = self.title.editable
        self.tabs.set_docks(docks)
        self.tabs.setVisible(tabbed)
        self.title.setVisible(not tabbed)

    def _rename_from_tab(self) -> None:
        self.tabs.setVisible(False)
        self.title.setVisible(True)
        self.title.start_editing()

    def _fill_arrange(self) -> None:
        self.host.fill_arrange(self.arrange_menu, self.dock)

    def add_tools(self, toolbar: QWidget) -> None:
        """Put a pane's own toolbar in the header (replacing one shown before, e.g. when the
        Filters pane is rebuilt for another document)."""
        if self.toolbar is not None and self.toolbar is not toolbar:
            self.row.removeWidget(self.toolbar)
            self.toolbar.setParent(None)
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
        if event.type() == QEvent.Type.PaletteChange and self._built:
            self._restyle()


class PanelDock(QDockWidget):
    """A dock holding a tool's panel; closing it closes the panel for good."""

    def __init__(self, title: str, host: DockHost) -> None:
        super().__init__(title, host.window)
        self.title_bar = PanelHeader(self, host)
        self.setTitleBarWidget(self.title_bar)

    def setWidget(self, widget: QWidget) -> None:
        super().setWidget(widget)
        if not isinstance(widget, ToolPanel):
            return
        bar = self.title_bar
        options = widget.options
        bar.settings.setVisible(options is not None)
        if options is not None:
            bar.settings.clicked.connect(lambda: options.popup(bar.settings))
        if widget.toolbar is not None:
            bar.add_tools(widget.toolbar)
        if widget.renamable:
            bar.title.on_rename = widget.set_custom_title
            bar.title.show_title(self.windowTitle())
        widget.attach_header(True)

    def closeEvent(self, event) -> None:
        widget = self.widget()
        if isinstance(widget, ToolPanel):
            widget.close_panel()
        super().closeEvent(event)
        self.deleteLater()
