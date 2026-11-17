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

from chemsuite.core.document import load_document
from chemsuite.core.table import Table
from tests.gui import gui

pytestmark = gui


def test_edit_marks_modified_and_undo_label(window) -> None:
    table = window.workspace.table
    assert not window.isWindowModified()
    table.set_cells("ID", [1], ["Z"], label="Edit ID")
    assert window.isWindowModified()
    assert window._undo.text() == "&Undo Edit ID" and window._undo.isEnabled()
    window.undo()
    assert not window.isWindowModified() and window._redo.isEnabled()


def test_status_bar_counts(window) -> None:
    window.workspace.selection.set([1, 2], "test")
    assert window._status.text() == "3 rows · 2 selected"


def test_save_marks_clean_and_writes_document(window, tmp_path: Path) -> None:
    window.workspace.table.set_cells("ID", [1], ["Z"])
    target = tmp_path / "doc.csd"
    assert window._save_to(target)  # waits for the background job
    assert not window.isWindowModified()
    assert window.workspace.path == target and "doc.csd" in window.windowTitle()
    loaded, _ = load_document(target)
    assert loaded.cell("ID", 1) == "Z"


def test_open_csv_parses_structures_in_background(qtbot, window, tmp_path: Path) -> None:
    csv = tmp_path / "mols.csv"
    csv.write_text("SMILES,ID\nCCO,a\nC1CC,b\nc1ccccc1,c\n")
    window.workspace.table.history.mark_clean()
    from chemsuite.app.document.main_window import Workspace

    window._install(Workspace(Table()))  # an empty window is reused by open_path
    window.open_path(csv)
    qtbot.waitUntil(lambda: len(window.workspace.table) == 3, timeout=10_000)
    table = window.workspace.table
    structure = next(n for n in table.column_names if table.kind(n) == "structure")
    column = table.column(structure)
    qtbot.waitUntil(lambda: not column.pending_mask().any(), timeout=20_000)
    assert column.invalid_mask().tolist() == [False, True, False]
    assert window.workspace.table.history.is_clean  # parsing is not an edit


def test_a_document_holding_only_a_protein_asks_to_save_and_is_not_reused(
    window, monkeypatch
) -> None:
    from PySide6.QtWidgets import QMessageBox

    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Asset, PutAsset

    table = Table()
    window._install(Workspace(table))
    assert window.workspace.is_empty
    table.edit("Add protein", [PutAsset("1ABC", None, Asset("1ABC", "protein", "1ABC", b""))])
    assert not window.workspace.is_empty  # File ▸ Open must not replace it
    asked: list[str] = []

    def question(_parent, _title, text, *_a, **_k):
        asked.append(text)
        return QMessageBox.StandardButton.Cancel

    monkeypatch.setattr(QMessageBox, "question", question)
    assert not window.close()
    assert asked == ["Save changes before closing?"]


def _drop(window, *paths: Path, enter: bool = False):
    """Drag *paths* onto *window*: just entering it, or dropping them."""
    from PySide6.QtCore import QMimeData, QPoint, QPointF, Qt, QUrl
    from PySide6.QtGui import QDragEnterEvent, QDropEvent

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(p)) for p in paths])
    args = (Qt.DropAction.CopyAction, mime, Qt.MouseButton.LeftButton)
    if enter:
        event = QDragEnterEvent(QPoint(10, 10), *args, Qt.KeyboardModifier.NoModifier)
        event.ignore()  # as Qt delivers it: a widget must accept it to take the drop
        window.dragEnterEvent(event)
    else:
        event = QDropEvent(QPointF(10, 10), *args, Qt.KeyboardModifier.NoModifier)
        event.ignore()
        window.dropEvent(event)
    return event


def test_files_dropped_on_an_empty_window_open_in_it(qtbot, window, tmp_path: Path) -> None:
    from chemsuite.app.document.main_window import Workspace

    smi = tmp_path / "mols.smi"
    smi.write_text("CCO a\nc1ccccc1 b\n")
    notes = tmp_path / "notes.docx"
    notes.write_text("not a compound file")
    assert window.acceptDrops()
    assert not _drop(window, smi, enter=True).isAccepted()  # the window holds a table already
    assert not _drop(window, smi).isAccepted()
    window.workspace.table.history.mark_clean()
    window._install(Workspace(Table()))
    assert not _drop(window, notes, enter=True).isAccepted()
    assert _drop(window, notes, smi, enter=True).isAccepted()
    assert _drop(window, notes, smi).isAccepted()
    qtbot.waitUntil(lambda: len(window.workspace.table) == 2, timeout=10_000)


def test_close_clears_the_last_window_and_keeps_it_open(window, monkeypatch) -> None:
    from PySide6.QtWidgets import QMessageBox

    from chemsuite.app.document import main_window

    assert main_window._windows == [window]
    window.workspace.table.set_cells("ID", [1], ["Z"])
    answers = [QMessageBox.StandardButton.Cancel, QMessageBox.StandardButton.Discard]
    monkeypatch.setattr(QMessageBox, "question", lambda *_a, **_k: answers.pop(0))
    before = window.workspace
    window._commands["file.close"].trigger()
    assert window.workspace is before  # Cancel: nothing changed
    window._commands["file.close"].trigger()
    assert window.workspace is not before and window.workspace.is_empty
    assert window.workspace.path is None and not window.isWindowModified()
    assert main_window._windows == [window]


def test_close_closes_a_window_when_another_is_open(window, qtbot) -> None:
    from chemsuite.app.document import main_window

    second = window.new_window()
    qtbot.addWidget(second)
    second._commands["file.close"].trigger()
    assert main_window._windows == [window]


def test_menus_have_no_quit_and_delete_selection_sits_with_the_selection(window) -> None:
    assert "file.quit" not in window._commands
    edit = next(a.menu() for a in window.menuBar().actions() if a.text() == "&Edit")
    names = [a.text() for a in edit.actions()]
    clear = names.index("Clear Selection")
    assert names[clear + 1] == "&Delete Selection"
    assert names[clear - 2 : clear] == ["Select &All", "Invert Selection"]
