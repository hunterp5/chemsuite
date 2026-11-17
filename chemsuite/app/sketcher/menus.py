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
"""The standalone Sketcher's menu bar (Tools ▸ Sketcher…): File (Save…, then the window's
own action, Add to Table) and Edit, the main window's Edit items for the drawing.

The canvas already handles the Edit keys itself (Delete erases what the pointer is over), so
those shortcuts are shown in the menu but act only while the menu is open."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import QFileDialog, QMenu, QMenuBar, QMessageBox, QWidget

from chemsuite.app.sketcher.canvas import SketchCanvas
from chemsuite.sketch import convert, files, reaction


def invert_selection(canvas: SketchCanvas) -> None:
    """Select every atom, arrow and plus sign that was not selected, and nothing else."""
    sketch = canvas.sketch
    atoms = set(sketch.atoms) - canvas.selected_atoms
    objects = (set(sketch.arrows) | set(sketch.pluses)) - canvas.selected_objects
    canvas.select_atoms(atoms)
    canvas.select_objects(objects)


def undo(canvas: SketchCanvas) -> None:
    canvas.cancel_preview()
    canvas.history.undo()


def redo(canvas: SketchCanvas) -> None:
    canvas.cancel_preview()
    canvas.history.redo()


def ask_path(parent: QWidget, canvas: SketchCanvas) -> str:
    """Where to save the drawing ("" if cancelled), with the suffix of the type chosen."""
    kinds = files.file_types(canvas.sketch)
    filters = ";;".join(f"{title} (*{suffix})" for suffix, title in kinds.items())
    name = "reaction" if files.is_reaction(canvas.sketch) else "structure"
    path, chosen = QFileDialog.getSaveFileName(
        parent, "Save Drawing", name + next(iter(kinds)), filters
    )
    if path and Path(path).suffix.lower() not in kinds:  # no suffix typed: the filter's
        path += next((s for s, title in kinds.items() if chosen.startswith(title)), "")
    return path


def save_drawing(parent: QWidget, canvas: SketchCanvas, path: str = "") -> bool:
    """Write the drawing to *path* (asked for when not given); False if it was not saved."""
    path = path or ask_path(parent, canvas)
    if not path:
        return False
    try:
        text = files.file_text(canvas.sketch, Path(path).suffix)
        Path(path).write_text(text, encoding="utf-8", newline="\n")
    except (ValueError, OSError, convert.SketchError, reaction.ReactionError) as err:
        QMessageBox.warning(parent, "Save Drawing", f"The drawing was not saved: {err}")
        return False
    return True


def build(
    parent: QWidget, canvas: SketchCanvas, add_text: str, add: Callable[[], None]
) -> tuple[QMenuBar, QAction, QAction]:
    """The menu bar, with the File actions that depend on the drawing: Save and *add*."""
    bar = QMenuBar(parent)
    bar.setNativeMenuBar(False)  # in the window on macOS too: it belongs to this window
    key = QKeySequence.StandardKey

    def action(menu: QMenu, title: str, slot, keys=None, *, shown_only=False) -> QAction:
        act = QAction(title, parent)
        act.triggered.connect(lambda _checked=False: slot())
        if keys is not None:
            act.setShortcut(QKeySequence(keys))
        if shown_only:  # the canvas handles these keys; the menu only shows them
            act.setShortcutContext(Qt.ShortcutContext.WidgetShortcut)
        elif keys is not None:
            parent.addAction(act)  # its keys work anywhere in the window
        menu.addAction(act)
        return act

    file_menu = bar.addMenu("&File")
    save = action(file_menu, "&Save…", lambda: save_drawing(parent, canvas), key.Save)
    added = action(file_menu, add_text, add, "Ctrl+Return")

    edit = bar.addMenu("&Edit")
    action(edit, "&Undo", lambda: undo(canvas), key.Undo, shown_only=True)
    action(edit, "&Redo", lambda: redo(canvas), key.Redo, shown_only=True)
    edit.addSeparator()
    action(edit, "Cu&t", canvas.cut, key.Cut, shown_only=True)
    action(edit, "&Copy", canvas.copy, key.Copy, shown_only=True)
    action(edit, "&Paste", canvas.paste, key.Paste, shown_only=True)
    edit.addSeparator()
    select_all = lambda: canvas.select_atoms(canvas.sketch.atoms)  # noqa: E731
    action(edit, "Select &All", select_all, key.SelectAll, shown_only=True)
    action(edit, "Invert Selection", lambda: invert_selection(canvas), "Ctrl+I")
    action(edit, "Clear Selection", canvas.clear_selection, "Ctrl+Shift+A", shown_only=True)
    action(edit, "&Delete Selection", canvas.erase_selection, key.Delete, shown_only=True)
    return bar, save, added
