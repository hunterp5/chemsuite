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
"""File ▸ Append… and File ▸ Merge…, another file's rows brought into the open document as
one undoable step (:mod:`chemsuite.core.table.combine`), and File ▸ Export."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtWidgets import QCheckBox, QComboBox, QFileDialog, QMessageBox

from chemsuite.core.document import SUFFIX as DOCUMENT_SUFFIX
from chemsuite.core.document import load_document
from chemsuite.core.logbook import audit
from chemsuite.core.table import Table
from chemsuite.core.table.columns import STRUCTURE
from chemsuite.core.table.combine import append, keys_of, merge
from chemsuite.core.wording import plural
from chemsuite.io.load import load_table, save_table

from ..shell.jobs import Job
from ..widgets.dialogs import ToolDialog

if TYPE_CHECKING:
    from .main_window import MainWindow

#: words in a column name that make it the likely key to merge on
KEY_WORDS = ("id", "name", "compound", "cmpd", "number", "code")


def add_actions(window: MainWindow, menu, action: Callable) -> None:
    """Append… and Merge… in the File menu (*action* makes a registered command)."""
    menu.addAction(action("file.append", "A&ppend…", lambda: choose(window, merge=False)))
    menu.addAction(action("file.merge", "&Merge…", lambda: choose(window, merge=True)))


def export(window: MainWindow, rows: str = "all") -> None:
    """Export every row, the ``"visible"`` ones (in view order) or the ``"selected"``."""
    from .main_window import EXPORT_FILTER, log

    workspace = window.workspace
    picked = {"visible": workspace.table.visible_oids, "selected": workspace.selection.as_array}
    oids = picked[rows]() if rows in picked else None
    if oids is not None and not len(oids):
        QMessageBox.information(window, "Export", f"No rows are {rows}.")
        return
    path, _ = QFileDialog.getSaveFileName(window, "Export", window._start_folder(), EXPORT_FILTER)
    if not path:
        return
    snapshot = workspace.table.snapshot()
    pool = window._runner.process_pool()

    def work(job: Job) -> None:
        save_table(snapshot, path, oids=oids, executor=pool, progress=job.report,
                   cancel=job.cancel_token)  # fmt: skip

    def done(_result: None) -> None:
        count = len(oids) if oids is not None else len(snapshot)
        window.show_message(f"Exported {Path(path).name}", 4000)
        log.info("exported %s (%d rows)", path, count)
        audit("export", path=path, rows=count, which=rows)

    window.jobs.run(f"Exporting {Path(path).name}…", work, done)


def choose(window: MainWindow, *, merge: bool) -> None:
    """Ask for a file, read it in the background, then append it or open the Merge dialog."""
    from .main_window import OPEN_FILTER

    title = "Merge" if merge else "Append"
    path, _ = QFileDialog.getOpenFileName(window, title, window._start_folder(), OPEN_FILTER)
    if path:
        read(window, Path(path), lambda other: _merge_dialog(window, other, Path(path))
             if merge else append_file(window, other, Path(path)))  # fmt: skip


def read(window: MainWindow, path: Path, then: Callable[[Table], None]) -> None:
    def work(job: Job) -> Table:
        if path.suffix.lower() == DOCUMENT_SUFFIX:
            return load_document(path, progress=job.report, cancel=job.cancel_token)[0]
        pool = window._runner.process_pool()
        return load_table(path, executor=pool, progress=job.report, cancel=job.cancel_token)

    window.jobs.run(f"Reading {path.name}…", work, then)


def append_file(window: MainWindow, other: Table, path: Path) -> None:
    oids = append(window.workspace.table, other, label=f"Append {path.name}")
    if len(oids):
        window.workspace.selection.set([int(o) for o in oids], "append")
    window.show_message(f"Appended {plural(len(oids), 'row')} from {path.name}")
    audit("append", path=path, rows=len(oids))


def _merge_dialog(window: MainWindow, other: Table, path: Path) -> MergeDialog:
    dialog = MergeDialog(window, other, path)
    dialog.show()
    return dialog


def guess_key(names: list[str], kinds: dict[str, str], wanted: str = "") -> str:
    """*wanted* if it is there, else an ID-like column, else the first structure column."""
    if wanted in names:
        return wanted
    for name in names:
        if any(word in name.lower().split() or name.lower() == word for word in KEY_WORDS):
            return name
    return next((n for n in names if kinds[n] == STRUCTURE), names[0] if names else "")


class MergeDialog(ToolDialog):
    """Match the file's rows to the table's by a column in each, then merge."""

    def __init__(self, window: MainWindow, other: Table, path: Path) -> None:
        super().__init__(window, f"Merge {path.name}", guide="append-and-merge")
        self.window_, self.other, self.path = window, other, path
        table = window.workspace.table
        mine = table.column_names
        theirs = other.column_names

        inputs = self.section("Input")
        self.key = QComboBox()
        self.key.addItems(mine)
        self.key.setCurrentText(guess_key(mine, {n: table.kind(n) for n in mine}))
        inputs.row(
            "Match",
            self.key,
            help="The table's column whose values identify a row. Structures match whatever "
            "way they were written; other values match as text.",
        )
        self.other_key = QComboBox()
        self.other_key.addItems(theirs)
        self.other_key.setCurrentText(
            guess_key(theirs, {n: other.kind(n) for n in theirs}, self.key.currentText())
        )
        inputs.row("With", self.other_key, help=f"The column of {path.name} to match it to.")

        method = self.section("Method")
        self.replace = QComboBox()
        self.replace.addItem("Replace them", True)
        self.replace.addItem("Keep them, fill blanks", False)
        method.row(
            "Values in the table",
            self.replace,
            help="What a matched row does with values it has already. A blank cell in the file "
            "never overwrites one.",
        )
        self.add_unmatched = QCheckBox("Add them as new rows")
        self.add_unmatched.setChecked(True)
        method.row("Rows with no match", self.add_unmatched)

        for signal in (self.key.currentIndexChanged, self.other_key.currentIndexChanged):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.key.currentText():
            return "The table has no columns to match."
        if not self.other_key.currentText():
            return f"{self.path.name} has no columns."
        return ""

    def verb(self) -> str:
        return "Merge"

    def run(self):
        window = self.window_
        table, other = window.workspace.table, self.other
        key, other_key = self.key.currentText(), self.other_key.currentText()
        replace, add = self.replace.currentData(), self.add_unmatched.isChecked()
        mine, theirs = table.snapshot().column(key), other.column(other_key)

        def work(_job: Job) -> tuple[list[str], list[str]]:
            return keys_of(mine), keys_of(theirs)

        def done(found: tuple[list[str], list[str]]) -> None:
            if key not in table.column_names or len(found[0]) != len(table):
                window.show_message("The table changed while the file was read: merge again")
                return
            result = merge(table, other, key, other_key, replace=replace, add_unmatched=add,
                           keys=found[0], other_keys=found[1],
                           label=f"Merge {self.path.name}")  # fmt: skip
            parts = [f"updated {plural(result.updated, 'row')}"]
            if len(result.added):
                parts.append(f"added {plural(len(result.added), 'row')}")
            elif result.unmatched:
                parts.append(f"left out {plural(result.unmatched, 'unmatched row')}")
            if result.new_columns:
                parts.append(f"added {plural(len(result.new_columns), 'column')}")
            window.show_message(f"Merged {self.path.name}: " + ", ".join(parts))
            audit("merge", path=self.path, rows=result.updated + len(result.added))

        return window.jobs.run("Matching rows…", work, done)
