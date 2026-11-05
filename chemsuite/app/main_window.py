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

"""The main window: one window per open document.

The window wires a :class:`Workspace` (table, selection, model) to a :class:`TableView` and
menus. It holds no tool logic: tools will plug in through the registry (phase 4).
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from PySide6.QtCore import QEventLoop, QPoint, QTimer
from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QWidget,
)

from chemsuite.core.document import SUFFIX as DOCUMENT_SUFFIX
from chemsuite.core.document import load_document, save_document
from chemsuite.core.io.load import (
    apply_parsed,
    load_table,
    parse_chunks,
    pending_structure_chunks,
    save_table,
)
from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import STRUCTURE

from .column_actions import ColumnActions
from .jobs import Job, JobRunner
from .table_model import TableModel
from .table_view import StructureImages, TableView

log = logging.getLogger(__name__)

OPEN_FILTER = (
    "All supported (*.csd *.sdf *.sd *.mol *.csv *.tsv *.txt *.smi *.smiles *.xlsx *.gz);;"
    "chemsuite documents (*.csd);;SD files (*.sdf *.sd *.mol *.sdf.gz);;"
    "CSV/TSV (*.csv *.tsv *.txt *.csv.gz *.tsv.gz);;SMILES (*.smi *.smiles);;Excel (*.xlsx)"
)
EXPORT_FILTER = "SD file (*.sdf);;CSV (*.csv);;TSV (*.tsv);;SMILES (*.smi);;Excel (*.xlsx)"

_windows: list[MainWindow] = []


class Workspace:
    """Everything that belongs to one open table."""

    def __init__(self, table: Table, path: Path | None = None) -> None:
        self.table = table
        self.path = path
        self.selection = Selection(table)
        self.model = TableModel(table)

    @property
    def is_empty(self) -> bool:
        return len(self.table) == 0 and not self.table.column_names

    def close(self) -> None:
        self.model.close()
        self.selection.close()


class MainWindow(QMainWindow):
    def __init__(
        self,
        runner: JobRunner,
        images: StructureImages,
        workspace: Workspace | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._runner = runner
        self._images = images
        self._workspace: Workspace | None = None
        self._view: TableView | None = None
        self._parse_job: Job | None = None
        self._parse_again = False
        self._columns: ColumnActions | None = None
        self._unsubscribe: Callable[[], None] | None = None
        self._status = QLabel()
        self.statusBar().addPermanentWidget(self._status)
        self.resize(1400, 860)
        self._build_menus()
        self._install(workspace or Workspace(Table()))
        _windows.append(self)

    # -- workspace ---------------------------------------------------------------------------
    @property
    def workspace(self) -> Workspace:
        assert self._workspace is not None
        return self._workspace

    @property
    def view(self) -> TableView:
        assert self._view is not None
        return self._view

    def _install(self, workspace: Workspace) -> None:
        self._teardown()
        self._workspace = workspace
        self._view = TableView(workspace.model, workspace.selection, self._images, self)
        self.setCentralWidget(self._view)
        self._columns = ColumnActions(
            self, workspace.table, lambda message: self.statusBar().showMessage(message, 4000)
        )
        self._view.headerMenuRequested.connect(self._show_header_menu)
        self._unsubscribe = workspace.table.subscribe(self._on_table_event)
        workspace.selection.subscribe(lambda *_: self._update_status())
        self._update_title()
        self._update_actions()
        self._update_status()
        self._start_parsing()

    def _teardown(self) -> None:
        if self._parse_job is not None:
            self._parse_job.cancel()
            self._parse_job = None
        if self._unsubscribe is not None:
            self._unsubscribe()
            self._unsubscribe = None
        if self._view is not None:
            self._view.close_view()
        if self._workspace is not None:
            self._workspace.close()

    def _on_table_event(self, event: Any) -> None:
        if isinstance(event, Changed):
            self._update_title()
            self._update_actions()
            self._update_status()
            if self._may_leave_unparsed(event):
                QTimer.singleShot(0, self._start_parsing)

    def _may_leave_unparsed(self, event: Changed) -> bool:
        """Could this change have left structure cells waiting to be parsed?"""
        table = self.workspace.table
        for change in event.changes:
            if isinstance(change, ev.ColumnKindChanged) and change.kind == STRUCTURE:
                return True
            if isinstance(change, ev.RowsInserted | ev.ColumnInserted):
                return True
            if (
                isinstance(change, ev.CellsChanged)
                and change.column in table
                and table.kind(change.column) == STRUCTURE
            ):
                return True
        return False

    @property
    def columns(self) -> ColumnActions:
        assert self._columns is not None
        return self._columns

    def _show_header_menu(self, column: int, pos: QPoint) -> None:
        name = self.workspace.model.column_name(column)
        menu = self.columns.build_menu(name, self.view.selected_column_names(), column)
        menu.exec(pos)
        menu.deleteLater()

    # -- menus -------------------------------------------------------------------------------
    def _action(
        self, title: str, slot: Callable[[], Any], shortcut: QKeySequence | str | None = None
    ) -> QAction:
        action = QAction(title, self)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
        action.triggered.connect(lambda _checked=False: slot())
        return action

    def _build_menus(self) -> None:
        bar = self.menuBar()
        file_menu = bar.addMenu("&File")
        file_menu.addAction(self._action("&New", self.new_window, QKeySequence.StandardKey.New))
        file_menu.addAction(self._action("&Open…", self.open_dialog, QKeySequence.StandardKey.Open))
        file_menu.addSeparator()
        file_menu.addAction(self._action("&Save", self.save, QKeySequence.StandardKey.Save))
        file_menu.addAction(
            self._action("Save &As…", self.save_as, QKeySequence.StandardKey.SaveAs)
        )
        file_menu.addAction(self._action("&Export…", lambda: self.export_dialog(False)))
        file_menu.addAction(self._action("Export Selected…", lambda: self.export_dialog(True)))
        file_menu.addSeparator()
        file_menu.addAction(self._action("&Close", self.close, QKeySequence.StandardKey.Close))
        file_menu.addAction(self._action("&Quit", self._quit, QKeySequence.StandardKey.Quit))

        edit = bar.addMenu("&Edit")
        self._undo = self._action("&Undo", self.undo, QKeySequence.StandardKey.Undo)
        self._redo = self._action("&Redo", self.redo, QKeySequence.StandardKey.Redo)
        edit.addActions([self._undo, self._redo])
        edit.addSeparator()
        edit.addAction(
            self._action("&Copy", lambda: self.view.copy_selection(), QKeySequence.StandardKey.Copy)
        )
        edit.addAction(self._action("&Paste", self.paste, QKeySequence.StandardKey.Paste))
        edit.addAction(self._action("Clear Cells", lambda: self.view.clear_selected_cells()))
        edit.addSeparator()
        edit.addAction(
            self._action(
                "Select &All", lambda: self.view.selectAll(), QKeySequence.StandardKey.SelectAll
            )
        )
        edit.addAction(
            self._action(
                "Invert Selection", lambda: self.workspace.selection.invert("menu"), "Ctrl+I"
            )
        )
        edit.addAction(
            self._action(
                "Clear Selection", lambda: self.workspace.selection.clear("menu"), "Ctrl+Shift+A"
            )
        )
        edit.addSeparator()
        edit.addAction(self._action("&Delete Rows", self.delete_rows, "Ctrl+Backspace"))

        data = bar.addMenu("&Data")
        table_menu = data.addMenu("&Table")
        table_menu.addAction(self._action("Add &Rows…", self.add_rows))
        table_menu.addAction(self._action("Add &Column…", self.add_column))
        table_menu.addSeparator()
        table_menu.addAction(self._action("Clear Sort", lambda: self.workspace.table.set_sort([])))

    def _update_actions(self) -> None:
        if self._workspace is None:
            return
        history = self.workspace.table.history
        self._undo.setEnabled(history.can_undo)
        self._redo.setEnabled(history.can_redo)
        self._undo.setText(f"&Undo {history.undo_label}" if history.can_undo else "&Undo")
        self._redo.setText(f"&Redo {history.redo_label}" if history.can_redo else "&Redo")

    def _update_title(self) -> None:
        if self._workspace is None:
            return
        path = self.workspace.path
        self.setWindowTitle(f"{path.name if path else 'Untitled'}[*] — chemsuite")
        self.setWindowModified(not self.workspace.table.history.is_clean)

    def _update_status(self) -> None:
        if self._workspace is None:
            return
        table = self.workspace.table
        shown = len(table.view())
        parts = [f"{len(table):,} rows"]
        if shown != len(table):
            parts.append(f"{shown:,} shown")
        parts.append(f"{len(self.workspace.selection):,} selected")
        self._status.setText(" · ".join(parts))

    # -- editing -----------------------------------------------------------------------------
    def undo(self) -> None:
        label = self.workspace.table.undo()
        if label:
            self.statusBar().showMessage(f"Undid {label}", 3000)

    def redo(self) -> None:
        label = self.workspace.table.redo()
        if label:
            self.statusBar().showMessage(f"Redid {label}", 3000)

    def paste(self) -> None:
        self.statusBar().showMessage(self.view.paste(), 4000)

    def delete_rows(self) -> None:
        oids = self.view.selected_oids()
        if len(oids):
            self.workspace.table.remove_rows(oids, label=f"Delete {len(oids)} row(s)")

    def add_rows(self) -> None:
        count, ok = QInputDialog.getInt(self, "Add Rows", "Number of empty rows:", 1, 1, 100_000)
        if ok:
            self.workspace.table.add_rows(count=count, label=f"Add {count} row(s)")

    def add_column(self) -> None:
        table = self.workspace.table
        name, ok = QInputDialog.getText(
            self, "Add Column", "Column name:", text=table.unique_name("New column")
        )
        name = name.strip()
        if ok and name:
            if name in table:
                QMessageBox.warning(self, "Add Column", f"A column named “{name}” already exists.")
                return
            table.add_column(name)

    # -- files -------------------------------------------------------------------------------
    def new_window(self) -> MainWindow:
        window = MainWindow(self._runner, self._images)
        window.show()
        return window

    def open_dialog(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(self, "Open", "", OPEN_FILTER)
        for path in paths:
            self.open_path(Path(path))

    def open_path(self, path: Path) -> None:
        """Open *path* here if this window is empty and unmodified, else in a new window."""
        target = self if self.workspace.is_empty and self.workspace.table.history.is_clean else None
        is_document = path.suffix.lower() == DOCUMENT_SUFFIX

        def work(job: Job) -> tuple[Table, Path | None]:
            if is_document:
                table, _extras = load_document(path, progress=job.report, cancel=job.cancel_token)
                return table, path
            table = load_table(
                path,
                executor=self._runner.process_pool(),
                progress=job.report,
                cancel=job.cancel_token,
            )
            return table, None  # an imported file is saved as a new document

        def opened(result: tuple[Table, Path | None]) -> None:
            table, doc_path = result
            window = target if target is not None and target in _windows else self.new_window()
            window._install(Workspace(table, doc_path))
            window.statusBar().showMessage(f"Opened {path.name}", 4000)

        self._run_with_progress(f"Opening {path.name}…", work, opened)

    def save(self) -> bool:
        if self.workspace.path is None:
            return self.save_as()
        return self._save_to(self.workspace.path)

    def save_as(self) -> bool:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save As", "", f"chemsuite document (*{DOCUMENT_SUFFIX})"
        )
        if not path:
            return False
        target = Path(path)
        if target.suffix.lower() != DOCUMENT_SUFFIX:
            target = target.with_suffix(DOCUMENT_SUFFIX)
        return self._save_to(target)

    def _save_to(self, path: Path) -> bool:
        workspace = self.workspace
        snapshot = workspace.table.snapshot()
        checkpoint = workspace.table.history.checkpoint()
        saved: list[bool] = []

        def work(job: Job) -> None:
            save_document(snapshot, path, progress=job.report, cancel=job.cancel_token)

        def done(_result: None) -> None:
            # The saved state is the snapshot's; edits made while saving stay unsaved.
            workspace.table.history.mark_clean(checkpoint)
            workspace.path = path
            self._update_title()
            self.statusBar().showMessage(f"Saved {path.name}", 4000)
            saved.append(True)

        self._run_with_progress(f"Saving {path.name}…", work, done, wait=True)
        return bool(saved)

    def export_dialog(self, selected_only: bool) -> None:
        workspace = self.workspace
        oids = workspace.selection.as_array() if selected_only else None
        if selected_only and not len(oids):
            QMessageBox.information(self, "Export Selected", "No rows are selected.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export", "", EXPORT_FILTER)
        if not path:
            return
        snapshot = workspace.table.snapshot()

        def work(job: Job) -> None:
            save_table(snapshot, path, oids=oids, progress=job.report, cancel=job.cancel_token)

        self._run_with_progress(
            f"Exporting {Path(path).name}…",
            work,
            lambda _r: self.statusBar().showMessage(f"Exported {Path(path).name}", 4000),
        )

    # -- background work ---------------------------------------------------------------------
    def _run_with_progress(
        self,
        label: str,
        work: Callable[[Job], Any],
        on_done: Callable[[Any], None],
        *,
        wait: bool = False,
    ) -> Job:
        dialog = QProgressDialog(label, "Cancel", 0, 0, self)
        dialog.setWindowTitle("chemsuite")
        dialog.setMinimumDuration(400)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)
        job = self._runner.submit(label, work)
        timer = QTimer(dialog)

        def poll() -> None:
            done, total = job.progress
            if total:
                dialog.setMaximum(1000)
                dialog.setValue(int(1000 * min(done, total) / total))

        timer.timeout.connect(poll)
        timer.start(100)
        dialog.canceled.connect(job.cancel)
        job.finished.connect(on_done)
        job.failed.connect(lambda err: self._show_error(label, err))
        job.cancelled.connect(lambda: self.statusBar().showMessage("Cancelled", 3000))
        job.done.connect(timer.stop)
        job.done.connect(dialog.close)
        if wait:  # e.g. Save from the close prompt: the caller needs the outcome
            dialog.setModal(True)
            ended: list[bool] = []
            loop = QEventLoop()
            job.done.connect(lambda: (ended.append(True), loop.quit()))
            if not ended:  # `done` is delivered on this thread, so it cannot be missed
                loop.exec()
        return job

    def _show_error(self, label: str, error: BaseException) -> None:
        QMessageBox.critical(self, "chemsuite", f"{label.rstrip('…')} failed:\n\n{error}")

    def _start_parsing(self) -> None:
        """Parse text structures (CSV, SMILES, Excel) in the background after opening."""
        if self._workspace is None:
            return
        if self._parse_job is not None and self._parse_job.running:
            self._parse_again = True  # look again once the running job ends
            return
        self._parse_again = False
        table = self.workspace.table
        chunks = pending_structure_chunks(table)
        if not chunks:
            return
        total = sum(len(c.oids) for c in chunks)

        def work(job: Job) -> None:
            done = 0
            for parsed in parse_chunks(
                chunks, executor=self._runner.process_pool(), cancel=job.cancel_token
            ):
                job.emit_partial(parsed)
                done += len(parsed.oids)
                job.report(done, total)

        invalid: list[int] = []

        def apply(parsed) -> None:
            if self._workspace is not None and self.workspace.table is table:
                invalid.append(apply_parsed(table, parsed))

        def finished(_r: None) -> None:
            bad = sum(invalid)
            message = f"Read {total:,} structures"
            self.statusBar().showMessage(message + (f"; {bad:,} invalid" if bad else ""), 6000)

        job = self._runner.submit("Reading structures", work)
        job.partial.connect(apply)
        job.finished.connect(finished)
        job.failed.connect(lambda err: self._show_error("Reading structures", err))
        job.done.connect(lambda: self._parse_again and self._start_parsing())
        self._parse_job = job

    # -- closing -----------------------------------------------------------------------------
    def closeEvent(self, event: QCloseEvent) -> None:
        if not self.workspace.table.history.is_clean and len(self.workspace.table):
            answer = QMessageBox.question(
                self,
                "chemsuite",
                "Save changes before closing?",
                QMessageBox.StandardButton.Save
                | QMessageBox.StandardButton.Discard
                | QMessageBox.StandardButton.Cancel,
            )
            if answer == QMessageBox.StandardButton.Cancel or (
                answer == QMessageBox.StandardButton.Save and not self.save()
            ):
                event.ignore()
                return
        self._teardown()
        if self in _windows:
            _windows.remove(self)
        event.accept()

    def _quit(self) -> None:
        for window in list(_windows):
            if not window.close():
                return
        QApplication.quit()
