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

"""Arranging panes (split, stack, tabs, sides) and saving the arrangement with the document."""

from __future__ import annotations

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMenu

from .conftest import gui

pytestmark = gui

RIGHT = Qt.DockWidgetArea.RightDockWidgetArea
LEFT = Qt.DockWidgetArea.LeftDockWidgetArea
BOTTOM = Qt.DockWidgetArea.BottomDockWidgetArea


@pytest.fixture
def data(window, qtbot):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["A", "B", "C"], "MW": ["100", "200", "300"], "pIC50": ["5", "6", "7"]},
        undoable=False,
    )
    window._install(Workspace(table))
    window.show()
    qtbot.waitExposed(window)
    return window


def two_plots(window):
    window.open_tool("scatter").set_spec(x="MW", y="pIC50")
    window.open_tool("histogram").set_spec(x="MW")
    return window._panels


def geometry(dock):
    return dock.mapTo(dock.window(), dock.rect().topLeft()), dock.size()


def test_new_panels_tab_and_split_below_stacks_them(data, qtbot):
    from chemsuite.app import dock_layout

    first, second = two_plots(data)
    assert second in data.tabifiedDockWidgets(first)
    dock_layout.split(data, first, second, Qt.Orientation.Vertical)
    qtbot.wait(50)
    assert not data.tabifiedDockWidgets(first)
    assert first.isVisible() and second.isVisible()
    (top, _), (below, _) = geometry(first), geometry(second)
    assert below.y() > top.y() and abs(below.x() - top.x()) < 5  # stacked in one column


def test_split_beside_a_pane_that_has_tabs_keeps_its_tabs(data, qtbot):
    from chemsuite.app import dock_layout

    first, second = two_plots(data)
    data.open_tool("scatter")
    third = data._panels[-1]
    assert set(data.tabifiedDockWidgets(first)) == {second, third}
    dock_layout.split(data, first, third, Qt.Orientation.Horizontal)
    qtbot.wait(50)
    assert data.tabifiedDockWidgets(first) == [second]
    assert geometry(third)[0].x() > geometry(first)[0].x()


def test_arrange_all_and_move_to_a_side(data, qtbot):
    from chemsuite.app import dock_layout

    first, second = two_plots(data)
    assert dock_layout.arrange(data, dock_layout.SIDE) == 2
    qtbot.wait(50)
    assert geometry(second)[0].x() > geometry(first)[0].x()
    dock_layout.arrange(data, dock_layout.GRID)
    dock_layout.arrange(data, dock_layout.TABS)
    assert second in data.tabifiedDockWidgets(first)
    dock_layout.move_to_area(data, second, BOTTOM)
    assert data.dockWidgetArea(second) == BOTTOM
    dock_layout.move_to_area(data, first, LEFT)
    assert data.dockWidgetArea(first) == LEFT


def test_arrange_menu_lists_the_other_panes(data):
    from chemsuite.app import dock_layout

    first, second = two_plots(data)
    menu = QMenu()
    dock_layout.fill_arrange_menu(menu, data, second)
    titles = [a.text() for a in menu.actions()]
    assert titles[:3] == ["Beside", "Below", "Tab With"]
    below = menu.actions()[1].menu()
    assert [a.text() for a in below.actions()] == [first.windowTitle()]
    below.actions()[0].trigger()
    assert not data.tabifiedDockWidgets(first)
    second.title_bar.arrange_menu.aboutToShow.emit()
    assert second.title_bar.arrange_menu.actions()


def test_layout_widths_and_search_are_saved_with_the_document(data, qtbot, tmp_path):
    from chemsuite.app import dock_layout
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    first, second = two_plots(data)
    dock_layout.split(data, first, second, Qt.Orientation.Vertical)
    data.show_browser()
    dock_layout.move_to_area(data, data.browser_dock, LEFT)
    data.view.set_column_widths({"MW": 211})
    data.search.focus()
    data.search.edit.setText("B")
    data.search.apply_now()
    target = tmp_path / "layout.csd"
    assert data._save_to(target)
    assert [d.objectName() for d in data._panels] == ["Panel0", "Panel1"]

    data._install(Workspace(Table()))
    data.browser_dock.hide()
    data.open_path(target)
    qtbot.waitUntil(lambda: len(data._panels) == 2, timeout=10_000)
    qtbot.wait(50)
    first, second = data._panels
    assert not data.tabifiedDockWidgets(first)  # still stacked, not tabbed
    assert geometry(second)[0].y() > geometry(first)[0].y()
    assert data.browser_dock.isVisible() and data.dockWidgetArea(data.browser_dock) == LEFT
    assert data.view.column_widths()["MW"] == 211
    assert data.search.is_open and data.search.edit.text() == "B"
    assert len(data.workspace.table.visible_oids()) == 1


def test_damaged_layout_is_ignored(data):
    from chemsuite.app import dock_layout

    dock_layout.restore(data, {"docks": "%%%", "widths": {"MW": "x"}, "search": 3})
    dock_layout.restore(data, "junk")


def test_tabbed_panes_show_their_tabs_in_the_header(data, qtbot):
    from PySide6.QtCore import QPointF
    from PySide6.QtWidgets import QTabBar

    from chemsuite.app import dock_layout

    first, second = two_plots(data)
    qtbot.waitUntil(lambda: second.title_bar.tabs.isVisible())
    header = second.title_bar
    assert header.tabs.titles() == [first.windowTitle(), second.windowTitle()]
    assert not header.title.isVisible()
    # Qt's own tab bar under the panes takes no room
    bars = [b for b in data.findChildren(QTabBar) if b.parent() is data and b.isVisible()]
    assert bars and all(b.height() == 0 for b in bars)
    # a click on the other tab brings that pane forward
    rect = header.tabs._rects[0]
    qtbot.mouseClick(header.tabs, Qt.MouseButton.LeftButton, pos=rect.center().toPoint())
    qtbot.waitUntil(lambda: first.title_bar.tabs.isVisible())
    assert first.title_bar.tabs.titles() == header.tabs.titles()
    # double-clicking the current tab renames the pane
    current = first.title_bar.tabs._rects[0].center()
    qtbot.mouseDClick(
        first.title_bar.tabs, Qt.MouseButton.LeftButton, pos=QPointF(current).toPoint()
    )
    assert first.title_bar.title.editing and first.title_bar.title.isVisible()
    first.title_bar.title.cancel()
    qtbot.waitUntil(lambda: first.title_bar.tabs.isVisible())
    # split apart: plain titles again
    dock_layout.split(data, first, second, Qt.Orientation.Vertical)
    qtbot.waitUntil(lambda: not second.title_bar.tabs.isVisible())
    assert second.title_bar.title.isVisible() and first.title_bar.title.isVisible()


def test_header_tabs_share_the_width_when_crowded():
    from chemsuite.app.panel_dock import _fit

    assert _fit([50, 60], 200, 44) == [50, 60]
    fitted = _fit([200, 60, 100], 240, 44)
    assert abs(sum(fitted) - 240) < 1 and fitted[1] == 60  # the widest give way first
