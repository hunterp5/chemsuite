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

from PySide6.QtCore import QEventLoop, QModelIndex, QPoint, Qt, QTimer
from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QFileDialog,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QVBoxLayout,
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
from chemsuite.core.lists import EXTRAS_KEY as LISTS_KEY
from chemsuite.core.lists import RowLists
from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import STRUCTURE

from .column_actions import ColumnActions
from .filters.cards import FilterCard
from .filters.panel import FilterPanel
from .find_replace import FindReplaceDialog
from .jobs import Job, JobRunner
from .list_actions import ListActions
from .search_bar import SearchBar
from .selection_browser import SelectionBrowser
from .split_join import JoinColumnsDialog, SplitColumnDialog
from .structure_editor import StructureEditor, copy_structure
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

    def __init__(
        self, table: Table, path: Path | None = None, extras: dict[str, Any] | None = None
    ) -> None:
        self.table = table
        self.path = path
        self.selection = Selection(table)
        self.model = TableModel(table)
        self.lists = RowLists(table)
        if extras and LISTS_KEY in extras:
            self.lists.load_json(extras[LISTS_KEY])
        self._saved_lists = self.lists.revision

    @property
    def is_empty(self) -> bool:
        return len(self.table) == 0 and not self.table.column_names

    @property
    def is_modified(self) -> bool:
        """Unsaved table edits or list changes."""
        return not self.table.history.is_clean or self.lists.revision != self._saved_lists

    def extras(self) -> dict[str, Any]:
        """Document state saved beside the table."""
        return {LISTS_KEY: self.lists.to_json()} if len(self.lists) else {}

    def mark_saved(self, checkpoint: Any, lists_revision: int) -> None:
        self.table.history.mark_clean(checkpoint)
        self._saved_lists = lists_revision

    def close(self) -> None:
        self.model.close()
        self.selection.close()
        self.lists.close()


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
        self._filters: FilterPanel | None = None
        self._search: SearchBar | None = None
        self._find: FindReplaceDialog | None = None
        self._lists: ListActions | None = None
        self._browser: SelectionBrowser | None = None
        self.filter_dock = QDockWidget("Filters", self)
        self.filter_dock.setObjectName("FilterDock")
        self.filter_dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.filter_dock)
        self.filter_dock.hide()
        self.browser_dock = QDockWidget("Selection Browser", self)
        self.browser_dock.setObjectName("BrowserDock")
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.browser_dock)
        self.browser_dock.hide()
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
        self._view = TableView(workspace.model, workspace.selection, self._images)
        self._search = SearchBar(workspace.table)
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._search)
        layout.addWidget(self._view, 1)
        self.setCentralWidget(central)
        self._filters = FilterPanel(workspace.table)
        self.filter_dock.setWidget(self._filters)
        self._columns = ColumnActions(
            self,
            workspace.table,
            lambda message: self.statusBar().showMessage(message, 4000),
            add_filter=self.filter_column,
            split=self.split_column,
            join=self.join_columns,
        )
        self._lists = ListActions(
            self,
            workspace.table,
            workspace.selection,
            workspace.lists,
            lambda message: self.statusBar().showMessage(message, 4000),
        )
        workspace.lists.subscribe(self._update_title)
        self._browser = SelectionBrowser(self._view, workspace.selection)
        self.browser_dock.setWidget(self._browser)
        self._view.headerMenuRequested.connect(self._show_header_menu)
        self._view.structureEditRequested.connect(self._on_structure_edit_requested)
        self._unsubscribe = workspace.table.subscribe(self._on_table_event)
        workspace.selection.subscribe(lambda *_: self._update_status())
        self._update_title()
        self._update_actions()
        self._update_status()
        self._start_parsing()

    def _teardown(self) -> None:
        if self._find is not None:
            self._find.close_dialog()
            self._find = None
        if self._browser is not None:
            self._browser.close_browser()
            self._browser.deleteLater()
            self._browser = None
        if self._filters is not None:
            self._filters.close()
            self._filters.deleteLater()
            self._filters = None
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
    def search(self) -> SearchBar:
        assert self._search is not None
        return self._search

    def show_find_replace(self) -> FindReplaceDialog:
        if self._find is None:
            self._find = FindReplaceDialog(self.view, self)
        self._find.show()
        self._find.raise_()
        self._find.activateWindow()
        self._find.find_edit.setFocus()
        self._find.find_edit.selectAll()
        return self._find

    @property
    def filters(self) -> FilterPanel:
        assert self._filters is not None
        return self._filters

    def filter_column(self, name: str, kind: str | None = None) -> FilterCard:
        """Show the filter panel with a new card for column *name*."""
        self.filter_dock.show()
        self.filter_dock.raise_()
        return self.filters.add_filter(name, kind)

    @property
    def lists(self) -> ListActions:
        assert self._lists is not None
        return self._lists

    @property
    def browser(self) -> SelectionBrowser:
        assert self._browser is not None
        return self._browser

    def show_browser(self) -> SelectionBrowser:
        self.browser_dock.show()
        self.browser_dock.raise_()
        self.browser.refresh(restart=True)
        self.browser.setFocus()
        return self.browser

    def show_all_rows(self) -> None:
        with self.workspace.table.batch():
            self.filters.remove_all()
            self.lists.show_all()

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
        edit.addAction(self._action("Copy Structure as Molfile", self.copy_structure))
        edit.addSeparator()
        edit.addAction(self._action("Edit &Structure…", self.edit_current_structure, "Ctrl+E"))
        edit.addSeparator()
        edit.addAction(
            self._action("&Search Rows", lambda: self.search.focus(), QKeySequence.StandardKey.Find)
        )
        find_replace = self._action("Find and &Replace…", self.show_find_replace)
        shortcuts = [*QKeySequence.keyBindings(QKeySequence.StandardKey.Replace)]
        shortcuts.append(QKeySequence("Ctrl+Alt+F"))
        find_replace.setShortcuts(shortcuts)
        edit.addAction(find_replace)
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
        table_menu.addAction(self._action("Add Co&mpound…", self.add_compound))
        table_menu.addAction(self._action("Add &Rows…", self.add_rows))
        table_menu.addAction(self._action("Add &Column…", self.add_column))
        table_menu.addSeparator()
        table_menu.addAction(self._action("&Split Column…", lambda: self.split_column()))
        table_menu.addAction(self._action("&Join Columns…", lambda: self.join_columns()))
        table_menu.addSeparator()
        table_menu.addAction(self._action("Clear Sort", lambda: self.workspace.table.set_sort([])))
        self._lists_menu = data.addMenu("&Lists")
        self._add_to_list = self._action(
            "&Add Selected to List…", lambda: self.lists.add_selected(), "Ctrl+L"
        )
        self._remove_from_list = self._action(
            "&Remove Selected from List…", lambda: self.lists.remove_selected()
        )
        self._lists_menu.aboutToShow.connect(
            lambda: self.lists.fill_menu(
                self._lists_menu, self._add_to_list, self._remove_from_list
            )
        )
        self._lists_menu.addAction(self._add_to_list)  # so the shortcut works before first show
        self.addAction(self._add_to_list)
        data.addSeparator()
        data.addAction(self._action("Show All Rows", self.show_all_rows))

        view = bar.addMenu("&View")
        toggle = self.filter_dock.toggleViewAction()
        toggle.setText("&Filters")
        toggle.setShortcut(QKeySequence("Ctrl+Shift+F"))
        view.addAction(toggle)
        browse = self._action("Selection &Browser", self.toggle_browser, "Ctrl+Shift+B")
        view.addAction(browse)

    def toggle_browser(self) -> None:
        if self.browser_dock.isVisible():
            self.browser_dock.hide()
        else:
            self.show_browser()

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
        self.setWindowModified(self.workspace.is_modified)

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

    # -- split / join ------------------------------------------------------------------------
    def _current_column(self) -> str | None:
        current = self.view.currentIndex()
        return self.workspace.model.column_name(current.column()) if current.isValid() else None

    def split_column(self, name: str | None = None) -> SplitColumnDialog | None:
        workspace = self.workspace
        if not workspace.table.column_names:
            self.statusBar().showMessage("The table has no columns", 4000)
            return None
        dialog = SplitColumnDialog(
            workspace.table, workspace.selection, name or self._current_column(), self
        )
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage(
                f"Added {len(dialog.result_names)} column(s)", 4000
            )
        )
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def join_columns(self, names: list[str] | None = None) -> JoinColumnsDialog | None:
        workspace = self.workspace
        if len(workspace.table.column_names) < 2:
            self.statusBar().showMessage("Joining needs at least two columns", 4000)
            return None
        if not names or len(names) < 2:
            names = self.view.selected_column_names() or names
        dialog = JoinColumnsDialog(workspace.table, workspace.selection, list(names or []), self)
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage(f"Added column {dialog.result_name}", 4000)
        )
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    # -- structures --------------------------------------------------------------------------
    def _structure_cell(self) -> tuple[int, str] | None:
        """(oid, column) of the current cell, or of the first structure cell in its row."""
        table = self.workspace.table
        model = self.workspace.model
        current = self.view.currentIndex()
        if not current.isValid():
            return None
        name = model.column_name(current.column())
        if table.kind(name) != STRUCTURE:
            name = next((n for n in table.column_names if table.kind(n) == STRUCTURE), "")
        return (model.oid_at(current.row()), name) if name else None

    def _on_structure_edit_requested(self, index: QModelIndex, typed: str) -> None:
        model = self.workspace.model
        self.edit_structure(model.oid_at(index.row()), model.column_name(index.column()), typed)

    def edit_current_structure(self) -> StructureEditor | None:
        cell = self._structure_cell()
        if cell is None:
            self.statusBar().showMessage("Pick a row in a table with a structure column", 4000)
            return None
        return self.edit_structure(*cell)

    def edit_structure(self, oid: int, name: str, typed: str = "") -> StructureEditor:
        """Open the structure editor on one cell; OK writes it back as one undoable edit."""
        table = self.workspace.table
        row = int(self.workspace.model.rows_of_oids([oid])[0])
        title = f"Edit Structure — {name}" + (f", row {row + 1}" if row >= 0 else "")
        dialog = StructureEditor(typed or table.cell(name, oid), title, self)
        if typed:  # typing into the cell starts a new structure; keep the cursor at the end
            cursor = dialog.edit.textCursor()
            cursor.movePosition(cursor.MoveOperation.End)
            dialog.edit.setTextCursor(cursor)

        def accepted() -> None:
            parsed = dialog.result_input
            if parsed is None or name not in table or oid not in table.oids:
                return
            if parsed.display != table.cell(name, oid) or parsed.kind == "molfile":
                table.set_cells(name, [oid], parsed.as_column(), label=f"Edit {name}")

        dialog.accepted.connect(accepted)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def add_compound(self) -> StructureEditor:
        """Enter a structure and append it as a new row."""
        table = self.workspace.table
        dialog = StructureEditor("", "Add Compound", self)

        def accepted() -> None:
            parsed = dialog.result_input
            if parsed is None or parsed.is_blank:
                return
            name = next((n for n in table.column_names if table.kind(n) == STRUCTURE), None)
            name = name or table.unique_name("Structure")
            oids = table.add_rows({name: parsed.as_column()}, label="Add compound")
            self.workspace.selection.set(oids.tolist(), "menu")
            row = int(self.workspace.model.rows_of_oids(oids)[0])
            if row >= 0:
                self.view.scrollTo(self.workspace.model.index(row, 0))

        dialog.accepted.connect(accepted)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def copy_structure(self) -> None:
        cell = self._structure_cell()
        table = self.workspace.table
        mol = None
        if cell is not None:
            oid, name = cell
            mol = table.column(name).mol(int(table.slots_of([oid])[0]))
        if copy_structure(mol):
            self.statusBar().showMessage("Copied structure as molfile", 3000)
        else:
            self.statusBar().showMessage("Pick a cell with a valid structure to copy", 4000)

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
        target = self if self.workspace.is_empty and not self.workspace.is_modified else None
        is_document = path.suffix.lower() == DOCUMENT_SUFFIX

        def work(job: Job) -> tuple[Table, Path | None, dict[str, Any]]:
            if is_document:
                table, extras = load_document(path, progress=job.report, cancel=job.cancel_token)
                return table, path, extras
            table = load_table(
                path,
                executor=self._runner.process_pool(),
                progress=job.report,
                cancel=job.cancel_token,
            )
            return table, None, {}  # an imported file is saved as a new document

        def opened(result: tuple[Table, Path | None, dict[str, Any]]) -> None:
            table, doc_path, extras = result
            window = target if target is not None and target in _windows else self.new_window()
            window._install(Workspace(table, doc_path, extras))
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
        extras, lists_revision = workspace.extras(), workspace.lists.revision
        saved: list[bool] = []

        def work(job: Job) -> None:
            save_document(
                snapshot, path, extras=extras, progress=job.report, cancel=job.cancel_token
            )

        def done(_result: None) -> None:
            # The saved state is the snapshot's; edits made while saving stay unsaved.
            workspace.mark_saved(checkpoint, lists_revision)
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
        if self.workspace.is_modified and len(self.workspace.table):
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
