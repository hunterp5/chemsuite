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
"""File ▸ Append…, File ▸ Merge… and File ▸ Export ▸ Visible Rows…."""

from __future__ import annotations

from PySide6.QtWidgets import QFileDialog

from tests.gui import gui

pytestmark = gui


def cells(table, name):
    return [table.cell(name, int(o)) for o in table.oids]


def _pick(monkeypatch, name, path) -> None:
    monkeypatch.setattr(QFileDialog, name, lambda *a, **k: (str(path), ""))


def test_the_file_menu_has_append_merge_and_visible_rows(window) -> None:
    file_menu = next(a.menu() for a in window.menuBar().actions() if a.text() == "&File")
    titles = [a.text() for a in file_menu.actions()]
    assert titles.index("A&ppend…") == titles.index("Open &Recent") + 1
    assert titles.index("&Merge…") == titles.index("A&ppend…") + 1
    export = next(a.menu() for a in file_menu.actions() if a.text() == "&Export")
    assert [a.text() for a in export.actions()][:3] == [
        "&Table…",
        "&Visible Rows…",
        "&Selected Rows…",
    ]


def test_append_adds_the_file_rows_in_one_step(window, qtbot, tmp_path, monkeypatch) -> None:
    from chemsuite.app.document import file_actions

    path = tmp_path / "more.csv"
    path.write_text("ID,pIC50,Series\nD,8,X\nE,9,Y\n")
    _pick(monkeypatch, "getOpenFileName", path)
    table = window.workspace.table
    file_actions.choose(window, merge=False)
    qtbot.waitUntil(lambda: len(table) == 5, timeout=10_000)
    assert cells(table, "ID") == ["A", "B", "C", "D", "E"]
    assert cells(table, "Series") == ["", "", "", "X", "Y"]
    assert window.workspace.selection.oids == {int(o) for o in table.oids[3:]}
    assert table.history.undo_label == "Append more.csv"
    table.undo()
    assert len(table) == 3 and "Series" not in table.column_names


def test_merge_matches_rows_by_a_column(window, qtbot, tmp_path, monkeypatch) -> None:
    from chemsuite.app.document import file_actions

    path = tmp_path / "assay.csv"
    path.write_text("Compound,pIC50,LogD\nB,6.5,2.1\nZ,4,1\n")
    _pick(monkeypatch, "getOpenFileName", path)
    file_actions.choose(window, merge=True)
    qtbot.waitUntil(lambda: bool(window.findChildren(file_actions.MergeDialog)), timeout=10_000)
    dialog = window.findChildren(file_actions.MergeDialog)[-1]
    assert dialog.key.currentText() == "ID" and dialog.other_key.currentText() == "Compound"
    dialog.add_unmatched.setChecked(False)
    table = window.workspace.table
    dialog.accept()
    qtbot.waitUntil(lambda: "LogD" in table.column_names, timeout=10_000)
    assert cells(table, "pIC50") == ["5", "6.5", "7"] and cells(table, "LogD") == ["", "2.1", ""]
    assert len(table) == 3 and table.history.undo_label == "Merge assay.csv"


def test_export_visible_rows_writes_only_the_rows_shown(window, qtbot, tmp_path, monkeypatch):
    from chemsuite.app.document import file_actions
    from chemsuite.core.table.filters import CategoryFilter

    table = window.workspace.table
    table.set_filter("test", CategoryFilter("ID", frozenset({"B"})))
    path = tmp_path / "shown.csv"
    _pick(monkeypatch, "getSaveFileName", path)
    file_actions.export(window, "visible")
    qtbot.waitUntil(path.exists, timeout=10_000)
    qtbot.waitUntil(lambda: "B" in path.read_text(), timeout=10_000)
    assert path.read_text().splitlines()[1:] == ["B,6"]


def test_export_smiles_asks_for_the_delimiter_and_header(window, qtbot, tmp_path, monkeypatch):
    from chemsuite.app.document import file_actions
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": ["A", "B"], "SMILES": ["CCO", "CCN"]}, undoable=False)
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    path = tmp_path / "list.smi"
    _pick(monkeypatch, "getSaveFileName", path)
    file_actions.export(window)
    dialog = window.findChildren(file_actions.TextExportDialog)[-1]
    assert dialog.names.currentData() == "ID" and not dialog.header.isChecked()
    dialog.delimiter.setCurrentIndex(dialog.delimiter.findData(","))
    dialog.header.setChecked(True)
    dialog.accept()
    qtbot.waitUntil(path.exists, timeout=10_000)
    assert path.read_text().splitlines() == ["SMILES,ID", "CCO,A", "CCN,B"]


def test_export_text_can_leave_out_the_header(window, qtbot, tmp_path, monkeypatch):
    from chemsuite.app.document import file_actions

    path = tmp_path / "table.txt"
    _pick(monkeypatch, "getSaveFileName", path)
    file_actions.export(window)
    dialog = window.findChildren(file_actions.TextExportDialog)[-1]
    assert dialog.delimiter.findData(" ") == -1  # a space only for SMILES
    dialog.delimiter.setCurrentIndex(dialog.delimiter.findData(";"))
    dialog.header.setChecked(False)
    dialog.accept()
    qtbot.waitUntil(path.exists, timeout=10_000)
    assert path.read_text().splitlines() == ["A;5", "B;6", "C;7"]
