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

"""Arranging panes around the table, and saving that arrangement with the document.

Panes (plots, the Structure Viewer) dock left of, right of or below the table, side by side,
stacked or as tabs. Two are exceptions: Filters is pinned at the right edge of the window at a
fixed width (``FILTER_WIDTH``), and every arrangement keeps it there; Advanced Search is always
its own floating window (``float_search``).

The other panes can be dragged by their header, and each header's **Arrange** menu places
the pane exactly:

* **Beside** / **Below** another pane: split that pane's space in two;
* **Tab With** another pane;
* **Left of / Right of / Below the Table**;
* **Float** in its own window.

The **Layout** menu (at the right of the menu bar; Window ▸ Layout on macOS) lays out every
docked pane at once: as tabs, side by side, stacked or in a grid; right of, left of or below
the table; and with the table or the panes given more of the window (:class:`Arrangement`).
The last choice is kept with the document.

Qt cannot split a pane that shares its space with tabs (``splitDockWidget`` tabs the new pane
instead), so splitting first takes the target's tab companions out and tabs them back onto the
target afterwards.

The arrangement (Qt's ``saveState``, which includes the splitter positions between the table
and the panes) is saved in the document with the column widths and the
search bar. Panels get stable object names (``Panel0``, ``Panel1``… in the order they are
saved), so the same arrangement applies when the document's panels are reopened.
"""

from __future__ import annotations

import base64
import math
from dataclasses import asdict, dataclass, replace
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QByteArray, QLineF, QRectF, Qt, QTimer
from PySide6.QtGui import QIcon, QPainter, QPalette, QPen, QPixmap
from PySide6.QtWidgets import QDockWidget, QMenu, QWidget

if TYPE_CHECKING:
    from .main_window import MainWindow

#: the pinned Filters pane's width
FILTER_WIDTH = 320
#: the Advanced Search window's first size
SEARCH_SIZE = (780, 560)
LAYOUT_KEY = "layout"
LAYOUT_VERSION = 1  # QMainWindow.saveState version tag
TABS, SIDE, STACKED, GRID = "tabs", "side", "stacked", "grid"
ARRANGEMENTS = (
    (TABS, "As Tabs"),
    (SIDE, "Side by Side"),
    (STACKED, "Stacked"),
    (GRID, "Grid"),
)
AREAS = (
    (Qt.DockWidgetArea.LeftDockWidgetArea, "Left of the Table"),
    (Qt.DockWidgetArea.RightDockWidgetArea, "Right of the Table"),
    (Qt.DockWidgetArea.BottomDockWidgetArea, "Below the Table"),
)
LEFT, RIGHT, BELOW = "left", "right", "below"
SIDES = (
    (RIGHT, "Right of the Table"),
    (LEFT, "Left of the Table"),
    (BELOW, "Below the Table"),
)
SIDE_AREAS = {
    LEFT: Qt.DockWidgetArea.LeftDockWidgetArea,
    RIGHT: Qt.DockWidgetArea.RightDockWidgetArea,
    BELOW: Qt.DockWidgetArea.BottomDockWidgetArea,
}
#: the panes' share of the window
SHARES = (
    (1 / 3, "Table Focus"),
    (1 / 2, "Balanced"),
    (2 / 3, "Pane Focus"),
)


@dataclass(frozen=True)
class Arrangement:
    """How the Layout menu lays out the panes."""

    how: str = TABS
    side: str = RIGHT
    share: float = 1 / 2

    @classmethod
    def from_saved(cls, saved: object) -> Arrangement:
        if not isinstance(saved, dict):
            return cls()
        how, side, share = saved.get("how"), saved.get("side"), saved.get("share")
        default = cls()
        ok_share = isinstance(share, int | float) and 0.1 <= share <= 0.9
        return cls(
            how if how in dict(ARRANGEMENTS) else default.how,
            side if side in SIDE_AREAS else default.side,
            float(share) if ok_share else default.share,
        )


def is_open(dock: QDockWidget) -> bool:
    """Open, even when behind another tab (Qt reports a background tab as not visible)."""
    return dock.toggleViewAction().isChecked()


def docked(window: MainWindow) -> list[QDockWidget]:
    """Visible, docked panes in a stable order (fixed panes, then panels as opened)."""
    fixed = [window.filter_dock, window.search_dock]
    return [
        d
        for d in [*fixed, *window.tool_panels.docks]
        if is_open(d) and not d.isFloating() and d.widget() is not None
    ]


def arrangeable(window: MainWindow) -> list[QDockWidget]:
    """Docked panes the Arrange menus may move (all but the pinned Filters pane)."""
    return [d for d in docked(window) if d is not window.filter_dock]


def pin_filters(window: MainWindow) -> None:
    """Put the Filters pane (and anything tabbed with it) back at the window's right edge."""
    dock = window.filter_dock
    panels = window.tool_panels
    if not is_open(dock) or panels.pinning:
        return
    panels.pinning = True  # re-adding the pane fires visibility signals that call back here
    try:
        _pin(window, dock)
    finally:
        panels.pinning = False


def _pin(window: MainWindow, dock: QDockWidget) -> None:
    dock.setFloating(False)
    companions = tab_group(window, dock)
    right = Qt.DockWidgetArea.RightDockWidgetArea
    current = next((d for d in [dock, *companions] if d.isVisible()), dock)
    for d in [dock, *companions]:
        _take(window, d)
    window.addDockWidget(right, dock, Qt.Orientation.Horizontal)
    _put(dock)
    for d in companions:
        window.tabifyDockWidget(dock, d)
        _put(d)
    current.raise_()


def widen(window: MainWindow, dock: QDockWidget, width: int) -> None:
    """Give a newly docked panel ``width`` pixels. Done now and again once the layout has
    settled: re-pinning Filters and showing the dock both share the space out again."""

    def resize() -> None:
        if dock.isVisible() and not dock.isFloating():
            window.resizeDocks([dock], [width], Qt.Orientation.Horizontal)

    resize()
    QTimer.singleShot(0, resize)


def float_search(window: MainWindow, *, show: bool = True) -> None:
    """Show Advanced Search as its own window (never docked). The first time, it gets a
    useful size over the right side of the main window; after that it keeps where you put it.
    With *show* False, only re-float it if it is open (after restoring a saved layout)."""
    dock = window.search_dock
    if not show and not is_open(dock):
        return
    first = not dock.property("placed")
    if not dock.isFloating():
        dock.setFloating(True)
    if first:
        frame = window.frameGeometry()
        width, height = SEARCH_SIZE
        width, height = min(width, frame.width()), min(height, frame.height())
        dock.resize(width, height)
        dock.move(frame.right() - width - 40, frame.top() + 80)
        dock.setProperty("placed", True)
    dock.show()
    dock.raise_()
    dock.activateWindow()


def tab_group(window: MainWindow, dock: QDockWidget) -> list[QDockWidget]:
    """The panes sharing *dock*'s space as tabs (not including *dock*)."""
    return [d for d in window.tabifiedDockWidgets(dock) if is_open(d)]


def _take(window: MainWindow, dock: QDockWidget) -> None:
    window.removeDockWidget(dock)  # also hides it


def _put(dock: QDockWidget) -> None:
    dock.show()


def split(
    window: MainWindow, target: QDockWidget, dock: QDockWidget, orientation: Qt.Orientation
) -> None:
    """Put *dock* beside (Horizontal) or below (Vertical) *target*, halving its space."""
    if dock is target:
        return
    dock.setFloating(False)
    _take(window, dock)
    if target.isFloating():
        target.setFloating(False)
    companions = tab_group(window, target)
    for d in companions:
        _take(window, d)
    window.splitDockWidget(target, dock, orientation)
    if dock in window.tabifiedDockWidgets(target):
        # Qt keeps a former tab group "tabbed" even with one pane left and tabs the new pane
        # into it; re-adding the target gives it a plain slot that splits.
        area = window.dockWidgetArea(target)
        _take(window, dock)
        _take(window, target)
        window.addDockWidget(area, target)
        _put(target)
        window.splitDockWidget(target, dock, orientation)
    _put(dock)
    for d in companions:
        window.tabifyDockWidget(target, d)
        _put(d)
    target.raise_()
    _even(window, [target, dock], orientation)
    dock.raise_()
    pin_filters(window)


def tab_with(window: MainWindow, target: QDockWidget, dock: QDockWidget) -> None:
    if dock is target:
        return
    dock.setFloating(False)
    _take(window, dock)
    window.tabifyDockWidget(target, dock)
    _put(dock)
    dock.raise_()
    pin_filters(window)


def move_to_area(window: MainWindow, dock: QDockWidget, area: Qt.DockWidgetArea) -> None:
    """Dock *dock* on one side of the table, beside whatever is already there."""
    dock.setFloating(False)
    _take(window, dock)
    there = [d for d in docked(window) if window.dockWidgetArea(d) == area and d is not dock]
    orientation = (
        Qt.Orientation.Horizontal
        if area == Qt.DockWidgetArea.BottomDockWidgetArea
        else Qt.Orientation.Vertical
    )
    window.addDockWidget(area, dock, orientation)
    _put(dock)
    if not there:  # a new side: give it a useful share of the window
        size = window.width() // 3 if orientation is Qt.Orientation.Vertical else None
        if size is not None:
            window.resizeDocks([dock], [max(320, size)], Qt.Orientation.Horizontal)
        else:
            window.resizeDocks([dock], [max(240, window.height() // 3)], Qt.Orientation.Vertical)
    dock.raise_()
    pin_filters(window)


def floating_panels(window: MainWindow) -> list[QDockWidget]:
    """Open panels in their own windows (Advanced Search is always one, so it is not)."""
    return [
        d
        for d in window.tool_panels.docks
        if is_open(d) and d.isFloating() and d.widget() is not None
    ]


def arrange(
    window: MainWindow,
    how: str,
    docks: list[QDockWidget] | None = None,
    side: str = RIGHT,
    share: float = 1 / 2,
) -> int:
    """Lay out every docked pane (or *docks*) on *side* of the table, the panes taking
    *share* of the window; returns how many moved."""
    docks = list(docks if docks is not None else arrangeable(window))
    if not docks:
        return 0
    area = SIDE_AREAS[side]
    for d in docks:
        d.setFloating(False)
        _take(window, d)
    first = docks[0]
    window.addDockWidget(area, first, Qt.Orientation.Horizontal)
    _put(first)
    rest = docks[1:]
    across, down = [first], [first]  # the panes along each direction, for sizing
    if how == TABS:
        for d in rest:
            window.tabifyDockWidget(first, d)
            _put(d)
        first.raise_()
    elif how in (SIDE, STACKED):
        orientation = Qt.Orientation.Horizontal if how == SIDE else Qt.Orientation.Vertical
        previous = first
        for d in rest:
            window.splitDockWidget(previous, d, orientation)
            _put(d)
            previous = d
        if how == SIDE:
            across = docks
        else:
            down = docks
    elif how == GRID:
        columns = math.ceil(math.sqrt(len(docks)))
        heads = docks[:columns]
        for previous, d in zip(heads, heads[1:], strict=False):
            window.splitDockWidget(previous, d, Qt.Orientation.Horizontal)
            _put(d)
        for k, d in enumerate(docks[columns:]):
            above = docks[k]  # the pane above in the same column
            window.splitDockWidget(above, d, Qt.Orientation.Vertical)
            _put(d)
        across, down = heads, docks[::columns]
    else:
        raise ValueError(f"unknown arrangement {how!r}")
    pin_filters(window)

    def size() -> None:  # now, and again once Qt has settled the new layout
        if side == BELOW:
            total = max(200, round(window.centralWidget().height() * share))
            window.resizeDocks(down, [total // len(down)] * len(down), Qt.Orientation.Vertical)
            _even(window, across, Qt.Orientation.Horizontal)
        else:
            room = window.width() - (FILTER_WIDTH if is_open(window.filter_dock) else 0)
            total = max(320, round(room * share))
            horizontal = Qt.Orientation.Horizontal
            window.resizeDocks(across, [total // len(across)] * len(across), horizontal)
            _even(window, down, Qt.Orientation.Vertical)

    size()
    QTimer.singleShot(0, size)
    return len(docks)


def apply_arrangement(window: MainWindow, choice: Arrangement, *, gather: bool = False) -> int:
    """Remember *choice* as the window's arrangement and lay the panes out by it; *gather*
    docks the floating panels too."""
    window.arrangement = choice
    docks = arrangeable(window) + (floating_panels(window) if gather else [])
    return arrange(window, choice.how, docks, choice.side, choice.share)


def _even(window: MainWindow, docks: list[QDockWidget], orientation: Qt.Orientation) -> None:
    """Share the space equally between *docks* along *orientation*."""
    if len(docks) < 2:
        return
    if orientation is Qt.Orientation.Horizontal:
        sizes = [d.width() for d in docks]
    else:
        sizes = [d.height() for d in docks]
    share = max(1, sum(sizes) // len(docks))
    window.resizeDocks(docks, [share] * len(docks), orientation)


# -- the Arrange menu ---------------------------------------------------------------------------
def fill_arrange_menu(menu: QMenu, window: MainWindow, dock: QDockWidget) -> None:
    """The Arrange menu of one pane's header (rebuilt each time it opens)."""
    menu.clear()
    others = arrangeable(window)
    others = [d for d in others if d is not dock]
    for title, orientation in (
        ("Beside", Qt.Orientation.Horizontal),
        ("Below", Qt.Orientation.Vertical),
    ):
        sub = menu.addMenu(title)
        sub.setEnabled(bool(others))
        for other in others:
            sub.addAction(
                other.windowTitle(),
                lambda o=other, r=orientation: split(window, o, dock, r),
            )
    tabs = menu.addMenu("Tab With")
    tabs.setEnabled(bool(others))
    group = set(tab_group(window, dock))
    for other in others:
        action = tabs.addAction(other.windowTitle(), lambda o=other: tab_with(window, o, dock))
        action.setCheckable(True)
        action.setChecked(other in group)
    menu.addSeparator()
    here = None if dock.isFloating() else window.dockWidgetArea(dock)
    for area, title in AREAS:
        action = menu.addAction(title, lambda a=area: move_to_area(window, dock, a))
        action.setCheckable(True)
        action.setChecked(area == here)
    menu.addSeparator()
    floating = menu.addAction("Float", lambda: dock.setFloating(not dock.isFloating()))
    floating.setCheckable(True)
    floating.setChecked(dock.isFloating())


def fill_layout_menu(menu: QMenu, window: MainWindow) -> None:
    """The Layout menu (rebuilt each time it opens): how, where and how large the panes are.
    Each choice keeps the other two and lays the panes out at once."""
    menu.clear()
    menu.setToolTipsVisible(True)
    current = window.arrangement
    has_panes = bool(arrangeable(window) or floating_panels(window))
    groups = (
        [(replace(current, how=how), title) for how, title in ARRANGEMENTS],
        [(replace(current, side=side), title) for side, title in SIDES],
        [(replace(current, share=share), title) for share, title in SHARES],
    )
    for k, group in enumerate(groups):
        if k:
            menu.addSeparator()
        for choice, title in group:
            action = menu.addAction(layout_icon(menu, choice), title)
            action.setCheckable(True)
            action.setChecked(choice == current)
            action.setEnabled(has_panes)
            action.triggered.connect(lambda _on=False, c=choice: apply_arrangement(window, c))
    menu.addSeparator()
    dock_all = menu.addAction(
        "Dock Floating Panes", lambda: apply_arrangement(window, current, gather=True)
    )
    dock_all.setEnabled(bool(floating_panels(window)))
    reset = menu.addAction(
        "Reset Layout", lambda: apply_arrangement(window, Arrangement(), gather=True)
    )
    reset.setToolTip("Every pane docked as tabs right of the table, half the window wide")
    reset.setEnabled(has_panes)


def layout_icon(widget: QWidget, choice: Arrangement) -> QIcon:
    """A small picture of *choice*: the window outlined, the panes filled."""
    width, height, ratio = 22, 16, widget.devicePixelRatioF()
    pixmap = QPixmap(round(width * ratio), round(height * ratio))
    pixmap.setDevicePixelRatio(ratio)
    pixmap.fill(Qt.GlobalColor.transparent)
    palette = widget.palette()
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    frame = QRectF(1.5, 1.5, width - 3, height - 3)
    if choice.side == BELOW:
        cut = frame.height() * choice.share
        panes = QRectF(frame.left(), frame.bottom() - cut, frame.width(), cut)
    else:
        cut = frame.width() * choice.share
        left = frame.left() if choice.side == LEFT else frame.right() - cut
        panes = QRectF(left, frame.top(), cut, frame.height())
    painter.fillRect(panes, palette.color(QPalette.ColorRole.Highlight))
    painter.setPen(QPen(palette.color(QPalette.ColorRole.HighlightedText), 1))
    for line in _pane_lines(panes, choice.how):
        painter.drawLine(line)
    painter.setPen(QPen(palette.color(QPalette.ColorRole.WindowText), 1))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRoundedRect(frame, 1.5, 1.5)
    painter.end()
    return QIcon(pixmap)


def _pane_lines(r: QRectF, how: str) -> list[QLineF]:
    """The lines between the panes drawn inside *r* for each arrangement."""
    mid_x, mid_y = r.center().x(), r.center().y()
    vertical = QLineF(mid_x, r.top(), mid_x, r.bottom())
    horizontal = QLineF(r.left(), mid_y, r.right(), mid_y)
    if how == SIDE:
        return [vertical]
    if how == STACKED:
        return [horizontal]
    if how == GRID:
        return [vertical, horizontal]
    tab = min(4.0, r.height() / 3)  # As Tabs: a tab strip along the top
    return [
        QLineF(r.left(), r.top() + tab, r.right(), r.top() + tab),
        QLineF(mid_x, r.top(), mid_x, r.top() + tab),
    ]


# -- saving ------------------------------------------------------------------------------------
def name_panels(window: MainWindow, restorable: list[QDockWidget]) -> None:
    """Name the panels that are saved ``Panel0``, ``Panel1``… in saving order, so a reopened
    document's panels match the saved arrangement."""
    for i, dock in enumerate(restorable):
        dock.setObjectName(f"Panel{i}")
    for dock in window.tool_panels.docks:
        if dock not in restorable:
            dock.setObjectName(f"Transient{id(dock)}")


def capture(window: MainWindow) -> dict[str, Any]:
    """The window's arrangement, column widths and search bar, as JSON-friendly data."""
    state = bytes(window.saveState(LAYOUT_VERSION).data())
    out: dict[str, Any] = {"docks": base64.b64encode(state).decode("ascii")}
    if window.arrangement != Arrangement():
        out["arrangement"] = asdict(window.arrangement)
    widths = window.view.column_widths()
    if widths:
        out["widths"] = widths
    if window.view.frozen_count:
        out["frozen"] = window.view.frozen_count
    if window.view.frozen_row_count:
        out["frozen_rows"] = window.view.frozen_row_count
    alignments = window.view.alignments()
    if alignments:
        out["alignments"] = alignments
    search = window.search
    if search.is_open or search.edit.text() or search.history:
        out["search"] = {
            "open": search.is_open,
            "text": search.edit.text(),
            "case": search.case.isChecked(),
            "whole": search.whole.isChecked(),
            "regex": search.regex.isChecked(),
            "column": search.column_name(),
            "history": list(search.history),
        }
    return out


def restore(window: MainWindow, saved: dict[str, Any]) -> None:
    """Apply what :func:`capture` saved; damaged parts are skipped."""
    if not isinstance(saved, dict):
        return
    window.arrangement = Arrangement.from_saved(saved.get("arrangement"))
    widths = saved.get("widths")
    if isinstance(widths, dict):
        window.view.set_column_widths(
            {str(k): int(v) for k, v in widths.items() if isinstance(v, int | float) and v > 0}
        )
    frozen = saved.get("frozen")
    if isinstance(frozen, int) and not isinstance(frozen, bool) and frozen > 0:
        window.view.set_frozen(frozen)
    frozen_rows = saved.get("frozen_rows")
    if isinstance(frozen_rows, int) and not isinstance(frozen_rows, bool) and frozen_rows > 0:
        window.view.set_frozen_rows(frozen_rows)
    alignments = saved.get("alignments")
    if isinstance(alignments, dict):
        table = window.workspace.table
        for column, core in alignments.items():
            if isinstance(column, str) and isinstance(core, str) and column in table:
                window.view.set_alignment(column, core)
    search = saved.get("search")
    if isinstance(search, dict):
        bar = window.search
        for key, action in (("case", bar.case), ("whole", bar.whole), ("regex", bar.regex)):
            action.setChecked(bool(search.get(key)))
        column = search.get("column")
        bar.set_column(column if isinstance(column, str) else None)
        history = search.get("history")
        if isinstance(history, list):
            bar.history = [str(item) for item in history if isinstance(item, str) and item][:10]
        text = str(search.get("text") or "")
        if search.get("open") or text:
            bar.open()
        bar.edit.setText(text)
        bar.apply_now()
    restore_docks(window, saved)


def restore_docks(window: MainWindow, saved: dict[str, Any]) -> None:
    """Apply the saved pane arrangement and splitter positions only."""
    docks = saved.get("docks") if isinstance(saved, dict) else None
    if isinstance(docks, str):
        try:
            data = base64.b64decode(docks.encode("ascii"), validate=True)
        except (ValueError, UnicodeError):
            return
        window.restoreState(QByteArray(data), LAYOUT_VERSION)
        pin_filters(window)  # documents saved before Filters was pinned may have moved it
        float_search(window, show=False)  # ... or docked Advanced Search
