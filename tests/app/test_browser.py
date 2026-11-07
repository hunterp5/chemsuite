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

import pytest
from PySide6.QtCore import Qt

from .conftest import gui

pytestmark = gui


@pytest.fixture
def browser(window, qtbot):
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": list("ABCD"), "SMILES": ["CCO", "c1ccccc1", "CC(=O)O", "not a smiles"]},
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    from chemsuite.app.main_window import Workspace

    window._install(Workspace(table))
    window.show()
    qtbot.waitExposed(window)
    return window.show_browser()


def ids(window, oids) -> list[str]:
    table = window.workspace.table
    return [table.cell("ID", int(o)) for o in oids]


def current_id(browser) -> str:
    return browser.table.cell("ID", browser.current_oid)


def test_browses_selected_rows_in_display_order(window, browser, qtbot) -> None:
    table = window.workspace.table
    window.workspace.selection.set(table.oids[[3, 0, 2]].tolist(), "test")
    table.set_sort([("ID", False)])
    browser.refresh()
    assert ids(window, browser.oids) == ["D", "C", "A"]
    browser.first()
    assert current_id(browser) == "D"
    browser.next()
    assert current_id(browser) == "C"
    assert browser.property_values() == {"ID": "C"}
    # The table cursor follows, without touching the selection.
    view_row = window.view.currentIndex().row()
    assert window.workspace.model.oid_at(view_row) == browser.current_oid
    assert len(window.workspace.selection) == 3


def test_scope_all_shown_and_keys(window, browser, qtbot) -> None:
    browser.set_scope("shown")
    assert len(browser.oids) == 4
    qtbot.keyClick(browser, Qt.Key.Key_End)
    assert current_id(browser) == "D"
    assert browser.preview.text() == "No structure"  # the invalid SMILES
    qtbot.keyClick(browser, Qt.Key.Key_Left)
    assert current_id(browser) == "C"
    assert not browser.preview.pixmap().isNull()
    qtbot.keyClick(browser, Qt.Key.Key_Space)
    assert browser.current_oid in window.workspace.selection
    qtbot.keyClick(browser, Qt.Key.Key_Home)
    assert current_id(browser) == "A" and not browser.prev_button.isEnabled()


def test_jump_box_and_follow_table_cursor(window, browser) -> None:
    browser.set_scope("shown")
    browser.position.setValue(3)
    assert current_id(browser) == "C"
    model = window.workspace.model
    window.view.setCurrentIndex(model.index(1, 0))
    assert current_id(browser) == "B"


def test_deselecting_while_browsing_selected_stays_in_place(window, browser, qtbot) -> None:
    table = window.workspace.table
    window.workspace.selection.set(table.oids[:3].tolist(), "test")
    browser.refresh()
    browser.go(1)
    browser.set_selected(False)  # B leaves the selected rows
    qtbot.waitUntil(lambda: len(browser.oids) == 2)
    assert current_id(browser) == "C"
    assert browser.count_label.text() == "of 2"


def test_empty_scope(window, browser) -> None:
    browser.set_scope("selected")
    assert browser.current_oid is None
    assert browser.preview.text() == "No rows selected"
    assert not browser.select_box.isEnabled()


def test_hidden_browser_catches_up_when_shown(window, browser, qtbot) -> None:
    window.browser_dock.hide()
    window.workspace.selection.set(window.workspace.table.oids[:2].tolist(), "test")
    assert len(browser.oids) == 0  # not refreshed while hidden
    window.browser_dock.show()
    qtbot.waitUntil(lambda: len(browser.oids) == 2)


def test_scope_is_a_check_box(browser) -> None:
    assert browser.only_selected.isChecked() and browser.scope_key() == "selected"
    browser.only_selected.setChecked(False)
    assert browser.scope_key() == "shown"
    browser.set_scope("selected")
    assert browser.only_selected.isChecked()
