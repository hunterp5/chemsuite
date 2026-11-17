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
"""Tools ▸ Notepad: notes kept with the document."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QTextCursor
from PySide6.QtWidgets import QMessageBox

from chemsuite.tools.notepad.notes import KEY, Notes
from chemsuite.tools.notepad.panel import NotepadPanel
from tests.gui import gui

pytestmark = gui


def notepads(window) -> list[NotepadPanel]:
    """The open note panes (a closed one's dock is hidden until Qt deletes it)."""
    docks = [d for d in window.tool_panels.docks if not d.isHidden()]
    return [d.widget() for d in docks if isinstance(d.widget(), NotepadPanel)]


def notes(window) -> Notes:
    return window.workspace.part(KEY, Notes)


def test_notes_part_round_trip_and_titles() -> None:
    part = Notes()
    first, second = part.new(), part.new()
    assert [part.title(first), part.title(second)] == ["Notepad", "Notepad 2"]
    part.set_title(second, "Assay plan")
    part.attach(first, lambda: "<p>live</p>")
    saved = part.capture()
    assert saved["notes"][first]["html"] == "<p>live</p>"
    again = Notes(saved)
    assert again.ids == [first, second] and again.title(second) == "Assay plan"
    assert Notes({"notes": {"x": "junk", 3: {}}}).ids == []  # damaged entries are skipped


def test_typing_marks_the_document_modified(window, qtbot) -> None:
    window.open_tool("notepad")
    panel = notepads(window)[0]
    workspace = window.workspace
    workspace.mark_saved(workspace.table.history.checkpoint(), workspace.lists.revision)
    assert not workspace.is_modified
    qtbot.keyClicks(panel.editor, "hello")
    assert workspace.is_modified
    html = notes(window).capture()["notes"][panel.note_id]["html"]
    assert "hello" in html
    body = html[html.index("<body") : html.index(">", html.index("<body"))]
    assert "font-family" not in body and "font-size" not in body  # each machine's own font


def test_formatting_survives_save_and_reopen(window, qtbot, tmp_path: Path) -> None:
    window.open_tool("notepad")
    panel = notepads(window)[0]
    editor = panel.editor
    qtbot.keyClicks(editor, "first")
    editor.selectAll()
    panel.actions["bold"].trigger()
    panel._set_size("18")
    cursor = editor.textCursor()
    cursor.movePosition(QTextCursor.MoveOperation.End)
    editor.setTextCursor(cursor)
    qtbot.keyClick(editor, Qt.Key.Key_Return)
    panel.actions["bullets"].trigger()
    qtbot.keyClicks(editor, "item")
    panel.set_custom_title("Plan")
    window.open_tool("notepad")
    second = notepads(window)[1]
    qtbot.keyClicks(second.editor, "other")
    assert second.title() == "Notepad"  # the first is now Plan

    path = tmp_path / "notes.csd"
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    assert window._save_to(path)
    window._install(Workspace(Table()))
    assert not notepads(window)
    window.open_path(path)
    qtbot.waitUntil(lambda: window.workspace.path == path, timeout=10_000)
    qtbot.waitUntil(lambda: len(notepads(window)) == 2, timeout=5_000)
    shown = {p.title(): p for p in notepads(window)}
    assert set(shown) == {"Plan", "Notepad"}
    document = shown["Plan"].editor.document()
    first_block = document.firstBlock()
    fmt = first_block.begin().fragment().charFormat()
    assert first_block.text() == "first" and fmt.fontWeight() >= 700 and fmt.fontPointSize() == 18
    item = first_block.next()
    assert item.text() == "item" and item.textList() is not None
    assert shown["Notepad"].editor.toPlainText() == "other"
    assert not window.workspace.is_modified


def test_closing_a_pane_deletes_its_note_asking_when_it_has_text(window, qtbot, monkeypatch):
    before = notes(window).revision
    window.open_tool("notepad")
    notepads(window)[0].parentWidget().title_bar.close_button.click()
    assert notes(window).ids == [] and notes(window).revision == before  # nothing to save
    window.open_tool("notepad")
    window.open_tool("notepad")
    empty, written = notepads(window)
    qtbot.keyClicks(written.editor, "keep me")
    asked = []

    def question(*args, **_kw):
        asked.append(args[2])
        return QMessageBox.StandardButton.Cancel

    monkeypatch.setattr(QMessageBox, "question", question)
    written.parentWidget().title_bar.close_button.click()
    assert asked == ["Delete note “Notepad 2”?"] and written.note_id in notes(window)
    empty.parentWidget().title_bar.close_button.click()
    assert len(asked) == 1 and empty.note_id not in notes(window)
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Discard)
    written.parentWidget().title_bar.close_button.click()
    assert notes(window).ids == []


def test_format_keys_and_undo_stay_in_the_note(window, qtbot) -> None:
    table = window.workspace.table
    table.set_cells("ID", [1], ["Z"], label="Edit ID")  # something the window could undo
    window.workspace.selection.set([1], "test")
    window.open_tool("notepad")
    panel = notepads(window)[0]
    window.show()
    window.activateWindow()
    qtbot.waitActive(window)
    panel.editor.setFocus()
    qtbot.waitUntil(panel.editor.hasFocus)
    qtbot.keyClicks(panel.editor, "ab")
    ctrl = Qt.KeyboardModifier.ControlModifier
    qtbot.keyClick(panel.editor, Qt.Key.Key_I, ctrl)
    assert panel.actions["italic"].isChecked()  # not the window's Invert Selection
    assert window.workspace.selection.as_array().tolist() == [1]
    qtbot.keyClick(panel.editor, Qt.Key.Key_Z, ctrl)
    assert panel.editor.toPlainText() == ""  # the note's typing is undone…
    assert table.cell("ID", 1) == "Z"  # …not the table's edit


def test_color_is_applied(window, qtbot, monkeypatch) -> None:
    from PySide6.QtWidgets import QColorDialog

    window.open_tool("notepad")
    panel = notepads(window)[0]
    monkeypatch.setattr(QColorDialog, "getColor", lambda *a, **k: QColor("#ff0000"))
    panel.actions["color"].trigger()
    qtbot.keyClicks(panel.editor, "red")
    fmt = panel.editor.document().firstBlock().begin().fragment().charFormat()
    assert fmt.foreground().color() == QColor("#ff0000")
