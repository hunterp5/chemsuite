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

import contextlib
import functools
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import shiboken6
from PySide6.QtCore import QEvent, QModelIndex, QPoint, Qt, QTimer, QUrl
from PySide6.QtGui import QAction, QCloseEvent, QDesktopServices, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QDockWidget,
    QFileDialog,
    QInputDialog,
    QMainWindow,
    QMenu,
    QMessageBox,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.structure_search import StructureQuery
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
from chemsuite.core.logbook import audit
from chemsuite.core.mols import smiles_of
from chemsuite.core.selection import Selection
from chemsuite.core.settings import Settings, add_recent
from chemsuite.core.structure_input import molblock_of, read_structure
from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import REACTION, STRUCTURE

from . import dock_layout, session
from .advanced_search import AdvancedSearchPanel
from .column_actions import ColumnActions
from .commands import apply_shortcuts, register
from .dock_layout import LAYOUT_KEY
from .document_windows import DocumentPanels
from .error_report import show_log_folder
from .filters.cards import FilterCard
from .filters.panel import FilterPanel
from .find_replace import FindReplaceDialog
from .job_host import JobHost
from .jobs import Job, JobRunner
from .jobs_ui import JobsIndicator, JobsWindow
from .list_actions import ListActions
from .log_ui import LogWindow
from .panel_dock import HIDE_DOCK_TABS, PanelDock, PanelHeader
from .search_bar import SearchBar
from .settings_hub import settings_hub
from .split_join import JoinColumnsDialog, SplitColumnDialog
from .status_bar import StatusBar
from .structure_clipboard import copy_structure
from .structure_copy import FORMATS as COPY_FORMATS
from .structure_copy import REACTION_FORMATS, copy_reactions, copy_structures
from .structure_search import StructureSearchDialog, StructureSearcher
from .table_model import TableModel
from .table_view import StructureImages, TableView
from .tools import BROWSER_TOOL, Tool, ToolContext, load_tools
from .widgets.panel import ToolPanel

log = logging.getLogger(__name__)

OPEN_FILTER = (
    "All supported (*.csd *.sdf *.sd *.mol *.csv *.tsv *.txt *.smi *.smiles *.xlsx *.rdf *.rd"
    " *.rxn *.gz);;"
    "chemsuite documents (*.csd);;SD files (*.sdf *.sd *.mol *.sdf.gz);;"
    "CSV/TSV (*.csv *.tsv *.txt *.csv.gz *.tsv.gz);;SMILES (*.smi *.smiles);;Excel (*.xlsx);;"
    "Reactions: RD and RXN files (*.rdf *.rd *.rxn *.rdf.gz)"
)
EXPORT_FILTER = (
    "SD file (*.sdf);;CSV (*.csv);;TSV (*.tsv);;SMILES (*.smi);;Excel (*.xlsx);;"
    "RD file, reactions (*.rdf);;RXN file, one reaction (*.rxn)"
)

_windows: list[MainWindow] = []
_jobs_windows: dict[int, JobsWindow] = {}


_log_windows: dict[int, LogWindow] = {}


def _log_window(runner: JobRunner) -> LogWindow:
    """One Log window per application; a job's entry opens it in the Jobs window."""
    window = _log_windows.get(id(runner))
    if window is None or getattr(window, "runner", None) is not runner:  # ids get reused

        def show_job(job_id: int) -> None:
            jobs = _jobs_window(runner)
            jobs.show()
            jobs.raise_()
            job = next((j for j in runner.jobs if j.id == job_id), None)
            if job is not None:
                jobs.select(job)

        window = _log_windows[id(runner)] = LogWindow(show_job=show_job)
        window.runner = runner
        runner.destroyed.connect(lambda *_: _log_windows.pop(id(runner), None))
    return window


def _jobs_window(runner: JobRunner) -> JobsWindow:
    """One Jobs window per runner (the application's), created on first use."""
    window = _jobs_windows.get(id(runner))
    if window is None or window.runner is not runner:
        window = _jobs_windows[id(runner)] = JobsWindow(runner)
        runner.destroyed.connect(lambda *_: _jobs_windows.pop(id(runner), None))
    return window


@functools.cache
def registered_tools() -> tuple[Tool, ...]:
    from chemsuite.tools import TOOL_MODULES

    return tuple(load_tools(TOOL_MODULES))


class Workspace:
    """Everything that belongs to one open table."""

    def __init__(
        self,
        table: Table,
        path: Path | None = None,
        extras: dict[str, Any] | None = None,
        *,
        unsaved: bool = False,
    ) -> None:
        self.table = table
        #: new content that was never saved (e.g. rows copied from another document)
        self._unsaved = unsaved
        self.path = path
        self.selection = Selection(table)
        self.model = TableModel(table)
        self.lists = RowLists(table)
        if extras and LISTS_KEY in extras:
            self.lists.load_json(extras[LISTS_KEY])
        #: the window session saved with the document (see ``session.PARTS``), put back when
        #: the window installs this workspace
        self.saved_session: dict[str, Any] = {
            key: value for key, value in (extras or {}).items() if key in session.KEYS
        }
        self._saved_lists = self.lists.revision

    @property
    def is_empty(self) -> bool:
        return not (len(self.table) or self.table.column_names or self.table.assets)

    @property
    def is_modified(self) -> bool:
        """Unsaved table edits or list changes, or content that was never saved."""
        return (
            self._unsaved
            or not self.table.history.is_clean
            or self.lists.revision != self._saved_lists
        )

    def extras(self) -> dict[str, Any]:
        """Document state saved beside the table."""
        return {LISTS_KEY: self.lists.to_json()} if len(self.lists) else {}

    def mark_saved(self, checkpoint: Any, lists_revision: int) -> None:
        self.table.history.mark_clean(checkpoint)
        self._saved_lists = lists_revision
        self._unsaved = False

    def close(self) -> None:
        self.model.close()
        self.selection.close()
        self.lists.close()


#: Panes meet at a 1 px line in the middle of a 5 px separator that is easy to grab.
_LINE = "stop:0 transparent, stop:0.4 transparent, stop:0.41 palette(mid), stop:0.59 palette(mid), stop:0.6 transparent, stop:1 transparent"  # noqa: E501
SEPARATOR_STYLE = (
    "QMainWindow::separator { width: 5px; height: 5px; }"
    "QMainWindow::separator:vertical { background: qlineargradient(x1:0, y1:0, x2:1, y2:0, "
    f"{_LINE}); }}"
    "QMainWindow::separator:horizontal { background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
    f"{_LINE}); }}"
)

#: a new panel takes half the window, but never less than this
PANEL_MIN_WIDTH = 360


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
        self.jobs = JobHost(self, runner, current=lambda: self._workspace,
                            name=lambda: self.document_name, message=self.show_message)  # fmt: skip
        self._images = images
        self._workspace: Workspace | None = None
        self._view: TableView | None = None
        self._parse_job: Job | None = None
        self._parse_again = False
        self._columns: ColumnActions | None = None
        self._filters: FilterPanel | None = None
        self._advanced: AdvancedSearchPanel | None = None
        self._search: SearchBar | None = None
        self._find: FindReplaceDialog | None = None
        self._searcher: StructureSearcher | None = None
        self._last_structure_query: StructureQuery | None = None
        self._lists: ListActions | None = None
        self.tool_panels = DocumentPanels(self)  # docked panels and tool windows
        self._dialogs: list[QDialog] = []  # open tool dialogs (they do not block the table)
        self.filter_dock = QDockWidget("Filters", self)
        self.filter_dock.setObjectName("FilterDock")
        # Filters is pinned to the window's right edge at a fixed width: it cannot be dragged,
        # floated or arranged elsewhere (dock_layout.pin_filters keeps it there).
        self.filter_dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea)
        self.filter_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable)
        self.filter_dock.setFixedWidth(dock_layout.FILTER_WIDTH)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.filter_dock)
        self.filter_dock.hide()
        # panes added while Filters was closed went right of it; opening it moves it back
        self.filter_dock.toggleViewAction().toggled.connect(
            lambda on: on and dock_layout.pin_filters(self)
        )
        self.search_dock = QDockWidget("Advanced Search", self)
        self.search_dock.setObjectName("AdvancedSearchDock")
        # Advanced Search is always its own floating window, so it never crowds the Filters
        # pane: it cannot be docked (dock_layout.float_search keeps it floating).
        self.search_dock.setAllowedAreas(Qt.DockWidgetArea.NoDockWidgetArea)
        self.search_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.search_dock)
        self.search_dock.hide()
        self.search_dock.setFloating(True)
        self.search_dock.toggleViewAction().toggled.connect(
            lambda on: on and dock_layout.float_search(self)
        )
        for dock in (self.filter_dock, self.search_dock):
            dock.setTitleBarWidget(PanelHeader(dock, self.tool_panels))  # one style for all
        self._unsubscribe: Callable[[], None] | None = None
        self.jobs_indicator = JobsIndicator(runner)
        self.jobs_indicator.clicked.connect(self.show_jobs)
        self.setStatusBar(StatusBar(self.jobs_indicator))
        self._status = self.statusBar().document
        self.resize(1400, 860)
        # Panes can sit side by side, stacked or as tabs in any dock area (dock_layout.py).
        self.setDockOptions(
            QMainWindow.DockOption.AnimatedDocks
            | QMainWindow.DockOption.AllowNestedDocks
            | QMainWindow.DockOption.AllowTabbedDocks
        )
        # Panes meet at a hairline, as in native macOS split views; the separator is wider than
        # the line so it is easy to grab. Tabbed panes show their tabs in their headers, so
        # Qt's tab bar below them is hidden.
        self.setStyleSheet(SEPARATOR_STYLE + HIDE_DOCK_TABS)
        self.tabifiedDockWidgetActivated.connect(
            lambda *_: self.tool_panels.schedule_header_refresh()
        )
        self._build_menus()
        hub = settings_hub()
        hub.changed.connect(self._apply_settings)
        self._install(workspace or Workspace(Table()))
        self._apply_settings(hub.settings)
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
        self._searcher = StructureSearcher(
            workspace.table,
            self.jobs.run,
            self._runner.process_pool,
            lambda column, highlight: self.view.set_highlight(column, highlight),
        )
        # hidden until Search Rows opens it
        self._search = SearchBar(workspace.table, searcher=self._searcher)
        self._search.openChanged.connect(self._sync_search_action)
        self._sync_search_action(False)
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._search)
        layout.addWidget(self._view, 1)
        self.setCentralWidget(central)
        self._filters = FilterPanel(workspace.table, searcher=self._searcher)
        self.filter_dock.setWidget(self._filters)
        self.filter_dock.titleBarWidget().add_tools(self._filters.toolbar)
        self._advanced = AdvancedSearchPanel(
            workspace.table,
            workspace.selection,
            workspace.lists,
            self._searcher,
            lambda message: self.statusBar().showMessage(message, 6000),
        )
        self.search_dock.setWidget(self._advanced)
        self._columns = ColumnActions(
            self,
            workspace.table,
            lambda message: self.statusBar().showMessage(message, 4000),
            add_filter=self.filter_column,
            search=lambda name: self.search.focus(name),
            selection=workspace.selection,
            freeze=lambda count: self.view.set_frozen(count),
            frozen=lambda: self.view.frozen_count,
            align=self.align_drawings,
            unalign=lambda name: self.view.set_alignment(name, None),
            aligned=lambda name: name in self.view.alignments(),
        )
        self._lists = ListActions(
            self,
            workspace.table,
            workspace.selection,
            workspace.lists,
            lambda message: self.statusBar().showMessage(message, 4000),
        )
        workspace.lists.subscribe(self._update_title)
        self._view.apply_appearance(settings_hub().settings.appearance)
        self._view.headerMenuRequested.connect(self._show_header_menu)
        self._view.structureEditRequested.connect(self._on_structure_edit_requested)
        self._view.structureDrawRequested.connect(self._on_structure_draw_requested)
        self._view.cellMenuRequested.connect(self._show_cell_menu)
        self._view.rowHeaderMenuRequested.connect(self._show_row_menu)
        self._view.deleteColumnsRequested.connect(lambda names: self.columns.delete(names))
        self._unsubscribe = workspace.table.subscribe(self._on_table_event)
        workspace.selection.subscribe(lambda *_: self._update_status())
        self._update_title()
        self._update_actions()
        self._update_status()
        self._start_parsing()
        saved = workspace.saved_session
        session.restore(self, saved)
        # The new table and panels are laid out on the next pass, which shares the space out
        # afresh: apply the saved sizes and scroll position once more after that pass.
        QTimer.singleShot(0, lambda ws=workspace: self._settle_session(ws))
        if LAYOUT_KEY not in saved:
            self._show_hiding_panes(saved)

    def _settle_session(self, workspace: Workspace) -> None:
        if self._workspace is workspace:
            session.settle(self, workspace.saved_session)

    def _show_hiding_panes(self, saved: dict[str, Any]) -> None:
        """A document saved without a layout: open the panes that hide rows, to say why."""
        if self._filters.cards:
            self.filter_dock.show()
        search = saved.get(session.ADVANCED_KEY)
        if isinstance(search, dict) and search.get("showing"):
            dock_layout.float_search(self)

    def show_dialog(self, dialog: QDialog) -> QDialog:
        """Show a tool dialog without blocking the window (see ``ToolContext.show_dialog``)."""
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self._dialogs.append(dialog)

        def forget(_result: int = 0) -> None:
            if dialog in self._dialogs:
                self._dialogs.remove(dialog)

        dialog.finished.connect(forget)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
        return dialog

    def _teardown(self) -> None:
        for dialog in list(self._dialogs):
            if shiboken6.isValid(dialog):
                dialog.reject()
        self._dialogs.clear()
        self.tool_panels.close_all()
        if self._find is not None:
            self._find.close_dialog()
            self._find = None
        if self._search is not None:
            self._search.close_bar()
        if self._advanced is not None:
            self._advanced.close_panel()
            self._advanced.deleteLater()
            self._advanced = None
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
            if isinstance(change, ev.ColumnKindChanged) and change.kind in (STRUCTURE, REACTION):
                return True
            if isinstance(change, ev.RowsInserted | ev.ColumnInserted):
                return True
            if (
                isinstance(change, ev.CellsChanged)
                and change.column in table
                and table.kind(change.column) in (STRUCTURE, REACTION)
            ):
                return True
        return False

    def _sync_search_action(self, is_open: bool) -> None:
        self._search_action.blockSignals(True)
        self._search_action.setChecked(is_open)
        self._search_action.blockSignals(False)
        if not is_open and self._view is not None:
            self.view.setFocus()  # back to the table after Esc or ✕

    @property
    def search(self) -> SearchBar:
        assert self._search is not None
        return self._search

    def show_find_replace(self) -> FindReplaceDialog:
        if self._find is None:
            self._find = FindReplaceDialog(self.view, self)
        self._find.rescope()  # the selection may have changed since it was last used
        self._find.show()
        self._find.raise_()
        self._find.activateWindow()
        self._find.find_edit.setFocus()
        self._find.find_edit.selectAll()
        return self._find

    @property
    def structure_searcher(self) -> StructureSearcher:
        assert self._searcher is not None
        return self._searcher

    def structure_search(self) -> StructureSearchDialog | None:
        """Edit ▸ Structure Search…"""
        table = self.workspace.table
        if not any(table.kind(n) == STRUCTURE for n in table.column_names):
            self.show_message("The table has no structure column")
            return None
        dialog = StructureSearchDialog(self, self._last_structure_query)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def remember_structure_query(self, query: StructureQuery) -> None:
        """The query the next Structure Search dialog starts with."""
        self._last_structure_query = query

    def clear_structure_highlights(self) -> None:
        for column in self.view.highlights():
            self.view.set_highlight(column, None)

    def current_column_name(self) -> str | None:
        return self._current_column()

    def current_structure_text(self, column: str | None = None) -> str:
        """SMILES of the current row's structure in *column* (default: the current cell's
        column if it holds structures, else the first structure column); "" if none."""
        table = self.workspace.table
        current = self.view.currentIndex()
        if not current.isValid():
            return ""
        if column is None or column not in table or table.kind(column) != STRUCTURE:
            cell = self._structure_cell()
            if cell is None:
                return ""
            column = cell[1]
        oid = self.workspace.model.oid_at(current.row())
        mol = table.column(column).mol(int(table.slots_of([oid])[0]))
        return smiles_of(mol) if mol is not None else ""

    def statusBar(self) -> StatusBar:
        return super().statusBar()  # type: ignore[return-value]

    def event(self, event: QEvent) -> bool:
        # Menu items' status tips: Qt would call the C++ QStatusBar.showMessage, which hides
        # the status bar's widgets; show them in our message label instead.
        if event.type() == QEvent.Type.StatusTip:
            self.statusBar().showMessage(event.tip())
            return True
        return super().event(event)

    def show_message(self, text: str, ms: int = 5000) -> None:
        self.statusBar().showMessage(text, ms)

    @property
    def advanced_search(self) -> AdvancedSearchPanel:
        assert self._advanced is not None
        return self._advanced

    def show_advanced_search(self) -> AdvancedSearchPanel:
        """Edit ▸ Advanced Search…: its own window, over the table."""
        dock_layout.float_search(self)
        return self.advanced_search

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

    def show_browser(self, oid: int | None = None) -> QWidget:
        """Data ▸ Browser: the Structure Viewer in 2D; a cell's Browse: at *oid*, as it is."""
        tool = next(t for t in registered_tools() if t.id == BROWSER_TOOL)
        assert tool.open_at is not None
        return tool.open_at(ToolContext(self), oid)

    def show_all_rows(self) -> None:
        with self.workspace.table.batch():
            self.filters.remove_all()
            self.lists.show_all()
            self.advanced_search.show_all()

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
        self,
        command_id: str,
        title: str,
        slot: Callable[[], Any],
        shortcut: QKeySequence | QKeySequence.StandardKey | str | list | None = None,
    ) -> QAction:
        """A menu action registered as a command, so Settings can change its shortcut."""
        defaults = (
            shortcut if isinstance(shortcut, list) else ([] if shortcut is None else [shortcut])
        )
        register(command_id, title, self._building, defaults)
        action = QAction(title, self)
        action.triggered.connect(lambda _checked=False: slot())
        self._commands[command_id] = action
        return action

    def _build_menus(self) -> None:
        bar = self.menuBar()
        self._commands: dict[str, QAction] = {}
        A = self._action

        self._building = "File"
        file_menu = bar.addMenu("&File")
        file_menu.addAction(A("file.new", "&New", self.new_window, QKeySequence.StandardKey.New))
        file_menu.addAction(
            A(
                "data.copy_shown_rows",
                "New from Shown Rows",
                self.copy_shown_rows_to_new_document,
            )
        )
        file_menu.addAction(
            A("file.open", "&Open…", self.open_dialog, QKeySequence.StandardKey.Open)
        )
        self._recent_menu = file_menu.addMenu("Open &Recent")
        self._recent_menu.aboutToShow.connect(self._fill_recent_menu)
        file_menu.addSeparator()
        file_menu.addAction(A("file.save", "&Save", self.save, QKeySequence.StandardKey.Save))
        file_menu.addAction(
            A("file.save_as", "Save &As…", self.save_as, QKeySequence.StandardKey.SaveAs)
        )
        export_menu = file_menu.addMenu("&Export")
        export_menu.setToolTipsVisible(True)
        export_menu.addAction(A("file.export", "&Table…", lambda: self.export_dialog(False)))
        export_menu.addAction(
            A("file.export_selected", "&Selected Rows…", lambda: self.export_dialog(True))
        )
        file_menu.addSeparator()
        file_menu.addAction(A("file.close", "&Close", self.close, QKeySequence.StandardKey.Close))
        quit_action = A("file.quit", "&Quit", self._quit, QKeySequence.StandardKey.Quit)
        quit_action.setMenuRole(QAction.MenuRole.QuitRole)
        file_menu.addAction(quit_action)

        self._building = "Edit"
        edit = bar.addMenu("&Edit")
        self._undo = A("edit.undo", "&Undo", self.undo, QKeySequence.StandardKey.Undo)
        self._redo = A("edit.redo", "&Redo", self.redo, QKeySequence.StandardKey.Redo)
        edit.addActions([self._undo, self._redo])
        edit.addSeparator()
        edit.addAction(
            A(
                "edit.copy",
                "&Copy",
                lambda: self.view.copy_selection(),
                QKeySequence.StandardKey.Copy,
            )
        )
        edit.addAction(A("edit.paste", "&Paste", self.paste, QKeySequence.StandardKey.Paste))
        copy_as = edit.addMenu("Copy Structure As")
        for fmt, title in COPY_FORMATS:
            command = "edit.copy_molfile" if fmt == "molfile" else f"edit.copy_{fmt}"
            copy_as.addAction(A(command, title, lambda f=fmt: self.copy_current_structure_as(f)))
        edit.addAction(
            A("edit.clear_cells", "Clear Cells", lambda: self.view.clear_selected_cells())
        )
        edit.addAction(
            A("edit.structure", "Edit &Structure…", self.edit_current_structure, "Ctrl+E")
        )
        edit.addSeparator()
        edit.addAction(
            A(
                "edit.select_all",
                "Select &All",
                lambda: self.view.selectAll(),
                QKeySequence.StandardKey.SelectAll,
            )
        )
        edit.addAction(
            A(
                "edit.invert_selection",
                "Invert Selection",
                lambda: self.workspace.selection.invert("menu"),
                "Ctrl+I",
            )
        )
        edit.addAction(
            A(
                "edit.clear_selection",
                "Clear Selection",
                lambda: self.workspace.selection.clear("menu"),
                "Ctrl+Shift+A",
            )
        )
        self._building = "Edit ▸ Lists"
        self._lists_menu = edit.addMenu("&Lists")
        self._add_to_list = A(
            "lists.add_selected",
            "&Add Selected to List…",
            lambda: self.lists.add_selected(),
            "Ctrl+L",
        )
        self._remove_from_list = A(
            "lists.remove_selected",
            "&Remove Selected from List…",
            lambda: self.lists.remove_selected(),
        )
        self._lists_menu.aboutToShow.connect(
            lambda: self.lists.fill_menu(
                self._lists_menu, self._add_to_list, self._remove_from_list
            )
        )
        self._lists_menu.addAction(self._add_to_list)  # so the shortcut works before first show
        self.addAction(self._add_to_list)
        self._building = "Edit"
        edit.addSeparator()
        edit.addAction(A("edit.delete_rows", "&Delete Rows", self.delete_rows, "Ctrl+Backspace"))
        edit.addSeparator()
        settings_action = A(
            "app.settings",
            "Settings…",
            self.open_settings,
            [QKeySequence.StandardKey.Preferences, "Ctrl+,"],
        )
        settings_action.setMenuRole(QAction.MenuRole.PreferencesRole)  # the app menu on macOS
        edit.addAction(settings_action)

        self._building = "Data ▸ Search"
        search = QMenu("&Search", self)  # placed in the Data menu below
        search.setToolTipsVisible(True)
        search.addAction(
            A(
                "edit.search",
                "&Search Rows",
                lambda: self.search.focus(),
                QKeySequence.StandardKey.Find,
            )
        )
        search.addAction(
            A(
                "edit.find_replace",
                "Find and &Replace…",
                self.show_find_replace,
                [QKeySequence.StandardKey.Replace, "Ctrl+Alt+F"],
            )
        )
        search.addAction(A("edit.advanced_search", "Ad&vanced Search…", self.show_advanced_search))
        search.addSeparator()
        search.addAction(
            A("edit.structure_search", "Structure Searc&h…", self.structure_search, "Ctrl+Shift+M")
        )
        search_tools = search.addSeparator()  # Reaction Search goes here
        search_tools.setVisible(False)
        search.addAction(
            A("edit.clear_highlight", "Clear Highlighting", self.clear_structure_highlights)
        )
        search.addSeparator()
        search.addAction(A("data.show_all_rows", "Show &All Rows", self.show_all_rows))

        self._building = "View"
        view = bar.addMenu("&View")
        self._search_action = QAction("&Search Bar", self, checkable=True)
        self._search_action.toggled.connect(
            lambda on: self.search.focus() if on else self.search.dismiss()
        )
        register("view.search_bar", "Search Bar", "View")
        self._commands["view.search_bar"] = self._search_action
        view.addAction(self._search_action)
        view.addSeparator()
        arrange = view.addMenu("&Arrange Panes")
        arrange.aboutToShow.connect(lambda: dock_layout.fill_arrange_all_menu(arrange, self))
        dock_layout.fill_arrange_all_menu(arrange, self)
        view.addSeparator()
        view.addAction(A("view.jobs", "&Jobs", self.show_jobs))
        view.addAction(A("view.log", "&Log", self.show_log))
        status = QAction("Status &Bar", self, checkable=True)
        status.toggled.connect(self._set_status_bar)
        register("view.status_bar", "Status Bar", "View")
        self._commands["view.status_bar"] = status
        self._status_bar_action = status
        view.addAction(status)

        self._building = "Tools ▸ Table"
        tools_menu = bar.addMenu("&Tools")
        tools_menu.setToolTipsVisible(True)
        table_menu = tools_menu.addMenu("&Table")
        table_menu.setToolTipsVisible(True)
        table_menu.addAction(A("data.add_compound", "Add Co&mpound…", self.add_compound))
        table_menu.addAction(A("data.add_rows", "Add &Rows…", self.add_rows))
        table_menu.addAction(A("data.add_column", "Add &Column…", self.add_column))
        table_menu.addSeparator()
        table_menu.addAction(A("data.split_column", "&Split Column…", lambda: self.split_column()))
        table_menu.addAction(A("data.join_columns", "&Join Columns…", lambda: self.join_columns()))
        table_menu.addSeparator()  # Calculator, Random… follow

        # the tool menus, in order (tools fill them; see chemsuite.tools.TOOL_MODULES)
        self._tool_actions: dict[str, QAction] = {}
        chemistry = bar.addMenu("&Chemistry")
        chemistry.setToolTipsVisible(True)
        self._building = "Data"
        data = bar.addMenu("&Data")  # Search, Plot and the panes first, then the analysis tools
        data.setToolTipsVisible(True)
        data.addMenu(search)
        plot = data.addMenu("&Plot")
        plot.setToolTipsVisible(True)
        toggle = self.filter_dock.toggleViewAction()
        toggle.setText("&Filters")
        register("view.filters", "Filters", "Data", ["Ctrl+Shift+F"])
        self._commands["view.filters"] = toggle
        data.addAction(toggle)
        data.addAction(A("view.browser", "&Browser", self.show_browser, "Ctrl+Shift+B"))
        data.addSeparator()
        menus: dict[tuple[str, ...], Any] = {
            ("File",): file_menu,
            ("File", "Export"): export_menu,
            ("Edit",): edit,
            ("View",): view,
            ("Tools",): tools_menu,
            ("Tools", "Table"): table_menu,
            ("Chemistry",): chemistry,
            ("Data",): data,
            ("Data", "Search"): search,
            ("Data", "Plot"): plot,
        }
        #: tools in these menus go before this action instead of at the end
        self._tool_slots = {("Data", "Search"): search_tools}
        for tool in registered_tools():
            self._add_tool(tool, menus)

        self._building = "Help"
        help_menu = bar.addMenu("&Help")
        help_menu.addAction(A("help.guide", "chemsuite &User Guide", self.open_user_guide))
        help_menu.addSeparator()
        help_menu.addAction(A("help.log_folder", "Show &Log Folder", show_log_folder))
        about = A("help.about", "&About chemsuite", self.show_about)
        about.setMenuRole(QAction.MenuRole.AboutRole)  # the app menu on macOS
        help_menu.addAction(about)

    def open_user_guide(self) -> None:
        from .dialogs.kit import GUIDE

        if GUIDE.exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(GUIDE)))
        else:
            self.show_message("The user guide is not installed with this copy of chemsuite")

    def show_about(self) -> None:
        import chemsuite

        QMessageBox.about(
            self,
            "About chemsuite",
            f"<b>chemsuite {chemsuite.__version__}</b><br>"
            "Compound tables, chemistry and SAR tools.<br><br>"
            "Free software under the GNU General Public License, version 3 or later.<br>"
            "Built on RDKit, Qt (PySide6), NumPy and pyqtgraph.",
        )

    def _add_tool(self, tool: Tool, menus: dict[tuple[str, ...], Any]) -> None:
        """Put *tool* in its menu, creating the menu path as needed."""
        for depth in range(1, len(tool.menu) + 1):
            path = tool.menu[:depth]
            if path not in menus:
                parent = self.menuBar() if depth == 1 else menus[path[:-1]]
                menus[path] = parent.addMenu(f"&{path[-1]}" if depth == 1 else path[-1])
                menus[path].setToolTipsVisible(True)
        action = QAction(tool.title, self)
        register(
            f"tool.{tool.id}",
            tool.title,
            " ▸ ".join(tool.menu),
            [tool.shortcut] if tool.shortcut else [],
        )
        missing = tool.missing()
        action.setToolTip(
            f"Needs {', '.join(missing)}, which is not installed" if missing else tool.tip
        )
        action.setEnabled(not missing)
        action.triggered.connect(lambda _checked=False, t=tool: self.open_tool(t))
        menu = menus[tool.menu]
        slot = self._tool_slots.get(tool.menu)
        if slot is not None:  # among the menu's own commands
            menu.insertAction(slot, action)
            if tool.separator_after:
                menu.insertSeparator(slot)
        else:
            menu.addAction(action)
            if tool.separator_after:
                menu.addSeparator()
        self._tool_actions[tool.id] = action
        self._commands[f"tool.{tool.id}"] = action

    # -- settings ----------------------------------------------------------------------------
    def _apply_settings(self, settings: Settings) -> None:
        apply_shortcuts(self._commands, settings)
        visible = settings.appearance.status_bar
        self.statusBar().setVisible(visible)
        self._status_bar_action.blockSignals(True)
        self._status_bar_action.setChecked(visible)
        self._status_bar_action.setEnabled(not settings_hub().locked("appearance.status_bar"))
        self._status_bar_action.blockSignals(False)
        if self._view is not None:
            self._view.apply_appearance(settings.appearance)

    def _set_status_bar(self, visible: bool) -> None:
        if settings_hub().user.appearance.status_bar != visible:
            settings_hub().update(lambda s: setattr(s.appearance, "status_bar", visible))

    def open_settings(self):
        from .settings_dialog import SettingsDialog

        dialog = SettingsDialog(self, self.jobs)
        dialog.open()
        return dialog

    def command_action(self, command_id: str) -> QAction:
        return self._commands[command_id]

    def _remember_path(self, path: Path) -> None:
        def change(settings: Settings) -> None:
            add_recent(settings, str(path))
            settings.files.last_directory = str(path.parent)

        settings_hub().update(change)

    def _start_folder(self) -> str:
        folder = settings_hub().settings.files.last_directory
        return folder if folder and Path(folder).is_dir() else ""

    def _fill_recent_menu(self) -> None:
        menu = self._recent_menu
        menu.clear()
        recent = settings_hub().settings.files.recent
        for text in recent:
            path = Path(text)
            action = menu.addAction(path.name)
            action.setToolTip(text)
            action.setStatusTip(text)
            if path.exists():
                action.triggered.connect(lambda _c=False, p=path: self.open_path(p))
            else:
                action.setEnabled(False)
                action.setText(f"{path.name} (missing)")
        if not recent:
            menu.addAction("No Recent Files").setEnabled(False)
        menu.addSeparator()
        clear = menu.addAction("Clear Recent")
        clear.setEnabled(bool(recent))
        clear.triggered.connect(
            lambda: settings_hub().update(lambda s: setattr(s.files, "recent", []))
        )

    def tool_action(self, tool_id: str) -> QAction:
        return self._tool_actions[tool_id]

    def open_tool(self, tool: Tool | str) -> Any:
        """Open a tool (by declaration or id) on this window's document."""
        if isinstance(tool, str):
            tool = next(t for t in registered_tools() if t.id == tool)
        return tool.open(ToolContext(self))

    def add_panel(self, widget: QWidget, title: str) -> PanelDock:
        # Beside the table (left | right), not under it: tables and plots both want height.
        # Later panels join the last one as tabs.
        previous = [d for d in self.tool_panels.docks if d.isVisible() and not d.isFloating()]
        dock = self.tool_panels.new_dock(title)
        dock.setObjectName(f"Panel{id(dock)}")
        dock.setWidget(widget)
        if isinstance(widget, ToolPanel):
            widget.titleChanged.connect(dock.setWindowTitle)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock, Qt.Orientation.Horizontal)
        if previous:
            self.tabifyDockWidget(previous[-1], dock)
        dock.show()
        dock_layout.pin_filters(self)  # the new panel goes left of Filters
        if not previous:  # half the window, sized once Filters is back in place
            dock_layout.widen(self, dock, max(PANEL_MIN_WIDTH, self.width() // 2))
        dock.raise_()
        return dock

    def panel_states(self) -> list[dict]:
        return self.tool_panels.states()

    def restore_panels(self, saved: list[dict]) -> None:
        self.tool_panels.restore(saved, registered_tools())

    @property
    def panels(self) -> list[QWidget]:
        return self.tool_panels.widgets()

    def _update_actions(self) -> None:
        if self._workspace is None:
            return
        table = self.workspace.table
        has_structures = any(table.kind(n) == STRUCTURE for n in table.column_names)
        has_reactions = any(table.kind(n) == REACTION for n in table.column_names)
        for tool in registered_tools():
            ready = (
                not tool.missing()
                and (has_structures or not tool.needs_structures)
                and (has_reactions or not tool.needs_reactions)
            )
            self._tool_actions[tool.id].setEnabled(ready)
        history = table.history
        self._undo.setEnabled(history.can_undo)
        self._redo.setEnabled(history.can_redo)
        self._undo.setText(f"&Undo {history.undo_label}" if history.can_undo else "&Undo")
        self._redo.setText(f"&Redo {history.redo_label}" if history.can_redo else "&Redo")

    @property
    def document_name(self) -> str:
        path = self._workspace.path if self._workspace is not None else None
        return path.name if path else "Untitled"

    def _update_title(self) -> None:
        if self._workspace is None:
            return
        self.setWindowTitle(f"{self.document_name}[*] — chemsuite")
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

    def _on_structure_edit_requested(self, index: QModelIndex, _typed: str) -> None:
        """F2, or typing in a structure cell, opens the sketcher (structures are only edited by
        drawing)."""
        model = self.workspace.model
        oid, name = model.oid_at(index.row()), model.column_name(index.column())
        if self.workspace.table.kind(name) == REACTION:
            self.draw_reaction(oid, name)
        else:
            self.draw_structure(oid, name)

    def _on_structure_draw_requested(self, index: QModelIndex) -> None:
        model = self.workspace.model
        oid, name = model.oid_at(index.row()), model.column_name(index.column())
        if self.workspace.table.kind(name) == REACTION:
            self.draw_reaction(oid, name)
        else:
            self.draw_structure(oid, name)

    def align_drawings(self, column: str):
        """Align Drawings… for *column*: ask for the core, then draw the column aligned."""
        from .align_drawings import AlignDrawingsDialog

        def current_structure() -> str:
            index = self.view.currentIndex()
            if not index.isValid():
                return ""
            oid = self.workspace.model.oid_at(index.row())
            return self.workspace.table.cell(column, oid)

        dialog = AlignDrawingsDialog(
            self,
            column,
            self.view.set_alignment,
            current_structure,
            self.view.alignments().get(column, ""),
        )
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def _show_cell_menu(self, index: QModelIndex, pos: QPoint) -> None:
        menu = self.cell_menu(index)
        menu.exec(pos)
        menu.deleteLater()

    def _show_row_menu(self, row: int, pos: QPoint) -> None:
        menu = self.row_menu(row)
        menu.exec(pos)
        menu.deleteLater()

    def row_menu(self, row: int) -> QMenu:
        """The right-click menu of the row numbers: for the selected rows, and freezing."""
        view = self.view
        count = len(view.selected_oids())
        rows = "Row" if count == 1 else f"{count:,} Rows"
        menu = QMenu(self)
        menu.addAction(f"Delete {rows}", self.delete_rows).setEnabled(count > 0)
        menu.addAction(f"Duplicate {rows}", self.duplicate_rows).setEnabled(count > 0)
        menu.addSeparator()
        freeze = menu.addAction("Freeze Rows up to Here", lambda: view.set_frozen_rows(row + 1))
        freeze.setToolTip("Keep this row and the ones above it in place while the rest scroll")
        freeze.setEnabled(view.frozen_row_count != row + 1)
        if view.frozen_row_count:
            menu.addAction("Unfreeze Rows", lambda: view.set_frozen_rows(0))
        menu.setToolTipsVisible(True)
        return menu

    def duplicate_rows(self) -> None:
        """Add a copy of each selected row (new rows go at the end, or by the sort) and select
        the copies."""
        oids = self.view.selected_oids()
        if not len(oids):
            return
        copies = self.workspace.table.duplicate_rows(oids, label=f"Duplicate {len(oids)} row(s)")
        self.workspace.selection.set(copies.tolist(), "menu")

    def cell_menu(self, index: QModelIndex) -> QMenu:
        """The right-click menu of a cell. Structure and reaction cells: the sketcher and Copy
        As (no plain Copy: Copy As says what you get); other cells: Copy."""
        model = self.workspace.model
        table = self.workspace.table
        name = model.column_name(index.column())
        oid = model.oid_at(index.row())
        menu = QMenu(self)
        menu.addAction("Browse", lambda: self.show_browser(oid))
        menu.addSeparator()
        kind = table.kind(name)
        if kind == STRUCTURE:
            menu.addAction("Edit in Sketcher…", lambda: self.draw_structure(oid, name))
            menu.addSeparator()
        elif kind == REACTION:
            menu.addAction("Edit in Sketcher…", lambda: self.draw_reaction(oid, name))
            menu.addSeparator()
        else:
            menu.addAction("Copy", self.view.copy_selection)
        if kind == STRUCTURE:
            copy_as = menu.addMenu("Copy As")
            for fmt, title in COPY_FORMATS:
                copy_as.addAction(title, lambda f=fmt: self.copy_structures_as(index, f))
        elif kind == REACTION:
            copy_as = menu.addMenu("Copy As")
            for fmt, title in REACTION_FORMATS:
                copy_as.addAction(title, lambda f=fmt: self.copy_reactions_as(index, f))
        menu.addAction("Paste", self.paste)
        menu.addSeparator()
        menu.addAction("Clear Cells", self.view.clear_selected_cells)
        return menu

    def _cell_oids(self, index: QModelIndex) -> np.ndarray:
        """The rows selected in *index*'s column, in view order, or just *index*'s row."""
        model = self.workspace.model
        name = model.column_name(index.column())
        oids = self.view.selected_cells_by_column().get(name)
        if oids is None or model.oid_at(index.row()) not in oids:
            return np.array([model.oid_at(index.row())])
        return oids[np.argsort(model.rows_of_oids(oids))]

    def copy_reactions_as(self, index: QModelIndex, fmt: str) -> str:
        """Copy the reactions selected in *index*'s column (or just that cell) as *fmt*."""
        table = self.workspace.table
        column = table.column(self.workspace.model.column_name(index.column()))
        slots = table.slots_of(self._cell_oids(index))
        message = copy_reactions([column.reaction(int(s)) for s in slots], fmt)
        self.statusBar().showMessage(message, 4000)
        return message

    def copy_structures_as(self, index: QModelIndex, fmt: str) -> str:
        """Copy the structures selected in *index*'s column (or just that cell) as *fmt*."""
        table = self.workspace.table
        column = table.column(self.workspace.model.column_name(index.column()))
        mols = [column.mol(int(slot)) for slot in table.slots_of(self._cell_oids(index))]
        message = copy_structures(mols, fmt)
        self.statusBar().showMessage(message, 4000)
        return message

    def draw_structure(self, oid: int, name: str):
        """Open the sketcher on one structure cell; OK writes it back (one undo step).

        A cell whose text is not a valid structure opens an empty sketcher; its text stays in the
        cell until a drawing replaces it.
        """
        from chemsuite.app.sketcher import SketcherDialog

        table = self.workspace.table
        column = table.column(name)
        mol = column.mol(int(table.slots_of([oid])[0]))
        row = int(self.workspace.model.rows_of_oids([oid])[0])
        title = f"Draw Structure — {name}" + (f", row {row + 1}" if row >= 0 else "")
        # open from the stored molecule: it keeps drawn coordinates and query atoms the cell
        # text (SMILES) cannot hold
        start = molblock_of(mol) if mol is not None else ""
        dialog = SketcherDialog(start, title, self)

        def accepted() -> None:
            if name not in table or oid not in table.oids:
                return
            if not dialog.canvas.history.can_undo:
                return  # nothing was drawn or changed
            parsed = read_structure(dialog.molblock)
            if parsed.is_blank and not table.cell(name, oid):
                return
            table.set_cells(name, [oid], parsed.as_column(), label=f"Draw {name}")

        dialog.accepted.connect(accepted)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def draw_reaction(self, oid: int, name: str):
        """Open the sketcher on one reaction cell; OK writes it back (one undo step). A cell
        whose text is not a valid reaction opens with that text in an inline editor instead."""
        from chemsuite.app.sketcher import SketcherDialog
        from chemsuite.core.io.rdf import rxn_blocks
        from chemsuite.core.reactions import parse_reaction
        from chemsuite.core.table.columns import ReactionColumn

        table = self.workspace.table
        column = table.column(name)
        slot = int(table.slots_of([oid])[0])
        if column.reaction(slot) is None and table.cell(name, oid).strip():
            row = int(self.workspace.model.rows_of_oids([oid])[0])
            if row >= 0:
                self.view.edit(
                    self.workspace.model.index(row, self.workspace.model.column_index(name))
                )
            return None
        row = int(self.workspace.model.rows_of_oids([oid])[0])
        title = f"Draw Reaction — {name}" + (f", row {row + 1}" if row >= 0 else "")
        cells = (column.blobs(np.array([slot])), column.displays(np.array([slot])))
        start = rxn_blocks(cells)[0] if table.cell(name, oid) else ""
        dialog = SketcherDialog(start, title, self, reaction_only=True)

        def accepted() -> None:
            if name not in table or oid not in table.oids:
                return
            if not dialog.canvas.history.can_undo:
                return  # nothing was drawn or changed
            if not dialog.rxn and not table.cell(name, oid):
                return
            # parsed now: the cell shows reaction SMILES and keeps the drawing's coordinates
            blob, shown = parse_reaction(dialog.rxn) if dialog.rxn else (b"", "")
            cell = ReactionColumn.from_parsed([blob], [shown])
            table.set_cells(name, [oid], cell, label=f"Draw {name}")

        dialog.accepted.connect(accepted)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def edit_current_structure(self):
        """Edit ▸ Edit Structure (Ctrl+E): the sketcher on the current structure cell."""
        cell = self._structure_cell()
        if cell is None:
            self.statusBar().showMessage("Pick a row in a table with a structure column", 4000)
            return None
        return self.draw_structure(*cell)

    def add_compound(self):
        """Draw a structure in the sketcher and append it as a new row."""
        from chemsuite.app.sketcher import SketcherDialog

        dialog = SketcherDialog("", "Add Compound", self, ok_text="Add")

        def accepted() -> None:
            if dialog.molblock:
                self.append_compound(dialog.molblock)

        dialog.accepted.connect(accepted)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def append_compound(self, text: str, parsed=None) -> int | None:
        """Add a structure as a new row (selected and scrolled to); its oid, or None if blank.

        It goes in the first structure column, or a new "Structure" column if there is none.
        """
        table = self.workspace.table
        parsed = parsed if parsed is not None else read_structure(text)
        if parsed.is_blank:
            return None
        name = next((n for n in table.column_names if table.kind(n) == STRUCTURE), None)
        name = name or table.unique_name("Structure")
        oids = table.add_rows({name: parsed.as_column()}, label="Add compound")
        self.workspace.selection.set(oids.tolist(), "menu")
        row = int(self.workspace.model.rows_of_oids(oids)[0])
        if row >= 0:
            self.view.scrollTo(self.workspace.model.index(row, 0))
        return int(oids[0])

    def copy_current_structure_as(self, fmt: str) -> None:
        """Copy As for the current structure cell (or the first structure in its row)."""
        cell = self._structure_cell()
        if cell is None:
            self.statusBar().showMessage("Pick a cell with a valid structure to copy", 4000)
            return
        oid, name = cell
        model = self.workspace.model
        row = int(model.rows_of_oids([oid])[0])
        if row < 0:
            return
        self.copy_structures_as(model.index(row, model.column_index(name)), fmt)

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
    def copy_shown_rows_to_new_document(self) -> MainWindow | None:
        """A new window holding the rows the filters and searches leave shown, in view order."""
        table = self.workspace.table
        shown = table.view()
        if not len(shown):
            self.show_message("No rows are shown")
            return None
        copy = table.rows_copy(shown)
        window = MainWindow(self._runner, self._images, Workspace(copy, unsaved=True))
        window.show()
        window.show_message(f"{len(shown):,} of {len(table):,} rows from {self.document_name}")
        return window

    def open_table(self, table: Table, message: str = "") -> MainWindow:
        """A new, unsaved document window holding *table* (a tool's results)."""
        window = MainWindow(self._runner, self._images, Workspace(table, unsaved=True))
        window.show()
        if message:
            window.show_message(message)
        return window

    def new_window(self) -> MainWindow:
        window = MainWindow(self._runner, self._images)
        window.show()
        return window

    def open_dialog(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(self, "Open", self._start_folder(), OPEN_FILTER)
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
            usable = (
                target is not None
                and target in _windows
                and target.workspace.is_empty
                and not target.workspace.is_modified
            )  # it may have been edited while the file was read
            window = target if usable else self.new_window()
            window._install(Workspace(table, doc_path, extras))
            window.statusBar().showMessage(f"Opened {path.name}", 4000)
            log.info("opened %s (%d rows, %d columns)", path, len(table), len(table.column_names))
            audit("open", path=path, rows=len(table))
            self._remember_path(path)

        self.jobs.run(f"Opening {path.name}…", work, opened)

    def save(self) -> bool:
        if self.workspace.path is None:
            return self.save_as()
        return self._save_to(self.workspace.path)

    def save_as(self) -> bool:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save As", self._start_folder(), f"chemsuite document (*{DOCUMENT_SUFFIX})"
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
        extras.update(session.capture(self))
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
            log.info("saved %s (%d rows)", path, len(snapshot))
            audit("save", path=path, rows=len(snapshot))
            self._remember_path(path)
            saved.append(True)

        self.jobs.run(f"Saving {path.name}…", work, done, wait=True)
        return bool(saved)

    def export_dialog(self, selected_only: bool) -> None:
        workspace = self.workspace
        oids = workspace.selection.as_array() if selected_only else None
        if selected_only and not len(oids):
            QMessageBox.information(self, "Export Selected", "No rows are selected.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export", self._start_folder(), EXPORT_FILTER)
        if not path:
            return
        snapshot = workspace.table.snapshot()

        def work(job: Job) -> None:
            save_table(
                snapshot,
                path,
                oids=oids,
                executor=self._runner.process_pool(),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(_result: None) -> None:
            rows = len(oids) if oids is not None else len(snapshot)
            self.statusBar().showMessage(f"Exported {Path(path).name}", 4000)
            log.info("exported %s (%d rows)", path, rows)
            audit("export", path=path, rows=rows, selected_only=selected_only)

        self.jobs.run(f"Exporting {Path(path).name}…", work, done)

    # -- background work ---------------------------------------------------------------------
    def show_jobs(self) -> JobsWindow:
        window = _jobs_window(self._runner)
        window.show()
        window.raise_()
        window.activateWindow()
        for w in _windows:
            w.jobs_indicator.acknowledge()
        return window

    def show_log(self) -> LogWindow:
        window = _log_window(self._runner)
        window.show()
        window.raise_()
        window.activateWindow()
        return window

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

        # The window may be gone by the time a chunk or the end arrives (closed mid-parse).
        def apply(parsed) -> None:
            if not shiboken6.isValid(self):
                return
            if self._workspace is not None and self.workspace.table is table:
                invalid.append(apply_parsed(table, parsed))

        def finished(_r: None) -> None:
            if not shiboken6.isValid(self):
                return
            bad = sum(invalid)
            message = f"Read {total:,} structures"
            self.statusBar().showMessage(message + (f"; {bad:,} invalid" if bad else ""), 6000)

        job = self._runner.submit("Reading structures", work)
        job.partial.connect(apply)
        job.finished.connect(finished)
        job.failed.connect(
            lambda err: shiboken6.isValid(self) and self.jobs.show_error("Reading structures", err)
        )
        job.done.connect(
            lambda: shiboken6.isValid(self) and self._parse_again and self._start_parsing()
        )
        self._parse_job = job

    # -- closing -----------------------------------------------------------------------------
    def closeEvent(self, event: QCloseEvent) -> None:
        jobs = [job for job in self._runner.jobs_of(self) if not job.background]
        if jobs:
            answer = QMessageBox.question(
                self,
                "chemsuite",
                f"{len(jobs)} job(s) for this document are still running or queued.\n\n"
                "Cancel them and close?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            )
            if answer != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
        if self.workspace.is_modified and not self.workspace.is_empty:
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
        self._runner.cancel_all(owner=self)
        if self._parse_job is not None:  # not owned (it never blocks closing), so stop it here
            self._parse_job.cancel()
        self._teardown()
        if self in _windows:  # closing twice must not disconnect twice
            _windows.remove(self)
            with contextlib.suppress(RuntimeError, TypeError):
                settings_hub().changed.disconnect(self._apply_settings)
        event.accept()

    def _quit(self) -> None:
        for window in list(_windows):
            if not window.close():
                return
        QApplication.quit()
