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

from .conftest import gui

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
    from chemsuite.app.main_window import Workspace

    window._install(Workspace(Table()))  # an empty window is reused by open_path
    window.open_path(csv)
    qtbot.waitUntil(lambda: len(window.workspace.table) == 3, timeout=10_000)
    table = window.workspace.table
    structure = next(n for n in table.column_names if table.kind(n) == "structure")
    column = table.column(structure)
    qtbot.waitUntil(lambda: not column.pending_mask().any(), timeout=20_000)
    assert column.invalid_mask().tolist() == [False, True, False]
    assert window.workspace.table.history.is_clean  # parsing is not an edit
