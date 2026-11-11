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

from pathlib import Path

from PySide6.QtWidgets import QInputDialog, QMessageBox

from .conftest import gui
from .test_column_actions import answer_question, answer_text, menu_titles, trigger

pytestmark = gui


def answer_item(monkeypatch, text: str, ok: bool = True) -> None:
    monkeypatch.setattr(QInputDialog, "getItem", lambda *a, **k: (text, ok))


def select(window, ids: list[str]) -> None:
    table = window.workspace.table
    oids = [o for o in table.oids.tolist() if table.cell("ID", o) in ids]
    window.workspace.selection.set(oids, "test")


def selected_ids(window) -> list[str]:
    table = window.workspace.table
    return sorted(table.cell("ID", o) for o in window.workspace.selection.oids)


def lists_menu(window):
    menu = window._lists_menu
    menu.aboutToShow.emit()
    return menu


def test_add_selected_creates_list_and_marks_modified(window, monkeypatch) -> None:
    select(window, ["A", "C"])
    answer_item(monkeypatch, "hits")
    assert window.lists.add_selected() == "hits"
    assert window.workspace.lists.count("hits") == 2
    assert window.isWindowModified()
    select(window, ["B"])
    window._add_to_list.trigger()  # Ctrl+L; same name adds to it
    assert window.workspace.lists.count("hits") == 3
    assert "hits (3)" in menu_titles(lists_menu(window))


def test_nothing_selected_asks_nothing(window) -> None:
    assert window.lists.add_selected() is None  # no dialog (it would raise)
    assert window.lists.remove_selected() is None


def test_list_menu_selects_and_shows_only(window, monkeypatch) -> None:
    select(window, ["A", "C"])
    answer_item(monkeypatch, "hits")
    window.lists.add_selected()
    window.workspace.selection.clear("test")
    trigger(lists_menu(window), "Select Rows")
    assert selected_ids(window) == ["A", "C"]
    trigger(lists_menu(window), "Show Only These Rows")
    assert len(window.workspace.table.view()) == 2
    assert "Stop Showing Only a List" in menu_titles(lists_menu(window))
    window.show_all_rows()
    assert len(window.workspace.table.view()) == 3


def test_remove_selected_from_list(window, monkeypatch) -> None:
    select(window, ["A", "B", "C"])
    answer_item(monkeypatch, "hits")
    window.lists.add_selected()
    select(window, ["B"])
    window.lists.remove_selected()
    assert window.workspace.lists.count("hits") == 2


def test_manage_dialog_renames_and_deletes(window, monkeypatch) -> None:
    lists = window.workspace.lists
    lists.add("a", window.workspace.table.oids[:1])
    lists.add("b", window.workspace.table.oids[:2])
    dialog = window.lists.manage()
    assert dialog.view.count() == 2
    dialog.view.setCurrentRow(1)
    answer_text(monkeypatch, "a")
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    dialog.rename_current()  # taken: warned, unchanged
    assert lists.names == ["a", "b"]
    answer_text(monkeypatch, "best")
    dialog.rename_current()
    assert lists.names == ["a", "best"] and dialog.current_name() == "best"
    answer_question(monkeypatch, QMessageBox.StandardButton.Yes)
    dialog.delete_current()
    assert lists.names == ["a"] and dialog.view.count() == 1
    dialog.select_current()
    assert len(window.workspace.selection) == 1
    dialog.accept()


def test_lists_are_saved_with_the_document(window, monkeypatch, tmp_path: Path) -> None:
    from chemsuite.core.document import load_document

    select(window, ["B", "C"])
    answer_item(monkeypatch, "hits")
    window.lists.add_selected()
    target = tmp_path / "t.csd"
    assert window._save_to(target)
    assert not window.isWindowModified()
    window.workspace.lists.rename("hits", "best")
    assert window.isWindowModified()  # list changes alone need saving
    _table, extras = load_document(target)
    assert extras["lists"][0]["name"] == "hits"
    assert len(extras["lists"][0]["oids"]) == 2


def test_opening_a_document_restores_lists(qtbot, window, tmp_path: Path) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.document import save_document
    from chemsuite.core.table import Table

    table = window.workspace.table
    target = tmp_path / "t.csd"
    save_document(table, target, extras={"lists": [{"name": "x", "oids": table.oids[:2].tolist()}]})
    window._install(Workspace(Table()))
    window.open_path(target)
    qtbot.waitUntil(lambda: "x" in window.workspace.lists, timeout=10_000)
    assert window.workspace.lists.count("x") == 2
    assert not window.isWindowModified()
