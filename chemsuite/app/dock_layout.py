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

Panes (plots, the Browser) dock left of, right of or below the table, side by side,
stacked or as tabs. Two are exceptions: Filters is pinned at the right edge of the window at a
fixed width (``FILTER_WIDTH``), and every arrangement keeps it there; Advanced Search is always
its own floating window (``float_search``).

The other panes can be dragged by their header, and each header's **Arrange** menu places
the pane exactly:

* **Beside** / **Below** another pane: split that pane's space in two;
* **Tab With** another pane;
* **Left of / Right of / Below the Table**;
* **Float** in its own window.

**View ▸ Arrange Panes** lays out every docked pane at once: as tabs, side by side, stacked or
in a grid.

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
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QByteArray, Qt
from PySide6.QtWidgets import QDockWidget, QMenu

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


def is_open(dock: QDockWidget) -> bool:
    """Open, even when behind another tab (Qt reports a background tab as not visible)."""
    return dock.toggleViewAction().isChecked()


def docked(window: MainWindow) -> list[QDockWidget]:
    """Visible, docked panes in a stable order (fixed panes, then panels as opened)."""
    fixed = [window.filter_dock, window.search_dock, window.browser_dock]
    return [
        d
        for d in [*fixed, *window._panels]
        if is_open(d) and not d.isFloating() and d.widget() is not None
    ]


def arrangeable(window: MainWindow) -> list[QDockWidget]:
    """Docked panes the Arrange menus may move (all but the pinned Filters pane)."""
    return [d for d in docked(window) if d is not window.filter_dock]


def pin_filters(window: MainWindow) -> None:
    """Put the Filters pane (and anything tabbed with it) back at the window's right edge."""
    dock = window.filter_dock
    if not is_open(dock) or getattr(window, "_pinning", False):
        return
    window._pinning = True  # re-adding the pane fires visibility signals that call back here
    try:
        _pin(window, dock)
    finally:
        window._pinning = False


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


def arrange(window: MainWindow, how: str, docks: list[QDockWidget] | None = None) -> int:
    """Lay out every docked pane (or *docks*) right of the table; returns how many moved."""
    docks = list(docks if docks is not None else arrangeable(window))
    if not docks:
        return 0
    area = Qt.DockWidgetArea.RightDockWidgetArea
    for d in docks:
        _take(window, d)
    first = docks[0]
    window.addDockWidget(area, first, Qt.Orientation.Horizontal)
    _put(first)
    rest = docks[1:]
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
        _even(window, docks, orientation)
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
        _even(window, heads, Qt.Orientation.Horizontal)
    else:
        raise ValueError(f"unknown arrangement {how!r}")
    total = window.width() // 2 if how != TABS else max(360, window.width() // 2)
    window.resizeDocks([first], [max(360, total)], Qt.Orientation.Horizontal)
    if how == SIDE:
        _even(window, docks, Qt.Orientation.Horizontal)
    pin_filters(window)
    return len(docks)


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


def fill_arrange_all_menu(menu: QMenu, window: MainWindow) -> None:
    menu.clear()
    count = len(arrangeable(window))
    for how, title in ARRANGEMENTS:
        action = menu.addAction(title, lambda h=how: arrange(window, h))
        action.setEnabled(count > 1 or (count == 1 and how == TABS))


# -- saving ------------------------------------------------------------------------------------
def name_panels(window: MainWindow, restorable: list[QDockWidget]) -> None:
    """Name the panels that are saved ``Panel0``, ``Panel1``… in saving order, so a reopened
    document's panels match the saved arrangement."""
    for i, dock in enumerate(restorable):
        dock.setObjectName(f"Panel{i}")
    for dock in window._panels:
        if dock not in restorable:
            dock.setObjectName(f"Transient{id(dock)}")


def capture(window: MainWindow) -> dict[str, Any]:
    """The window's arrangement, column widths and search bar, as JSON-friendly data."""
    state = bytes(window.saveState(LAYOUT_VERSION).data())
    out: dict[str, Any] = {"docks": base64.b64encode(state).decode("ascii")}
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
