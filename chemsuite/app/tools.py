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
"""Tools: self-contained features that plug into the main window's menus.

Each tool declares a :class:`Tool` in its own module; the window builds its menus from the
declarations and hands each tool a :class:`ToolContext` when it opens. The context is the
whole API a tool gets: read the table, see the selection, run background work, write result
columns. Tools never reach into the window, so adding one never means editing it.
"""

from __future__ import annotations

import contextlib
import importlib
import importlib.util
from collections.abc import Callable, Collection, Mapping, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from PySide6.QtCore import QItemSelectionModel, QSize
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QAbstractItemView, QDialog, QDockWidget, QWidget

from chemsuite.core.lists import RowLists
from chemsuite.core.selection import Selection
from chemsuite.core.table import Table
from chemsuite.core.table.columns import REACTION, STRUCTURE, Column

from .jobs import QUEUE, Job

if TYPE_CHECKING:
    from .main_window import MainWindow


@dataclass(frozen=True)
class Tool:
    """A menu entry and what it opens.

    ``menu`` is the menu path, e.g. ``("Chemistry",)`` or ``("Chemistry", "Calculate")``.
    ``requires`` names importable modules; a tool whose module is missing is shown disabled
    with the reason. ``needs_structures`` disables the tool while the table has no structure
    column, ``needs_reactions`` while it has no reaction column. ``restore`` reopens a panel
    saved with the document: a panel that should be saved has ``state() -> dict`` and a
    ``tool_id`` naming the tool whose ``restore`` takes that dict. ``aliases`` are former ids
    whose saved panels this tool restores (a tool that was renamed or merged into another).
    """

    id: str
    title: str
    menu: tuple[str, ...]
    open: Callable[[ToolContext], Any]
    requires: tuple[str, ...] = ()
    needs_structures: bool = False
    needs_reactions: bool = False
    shortcut: str | None = None
    tip: str = ""
    restore: Callable[[ToolContext, dict], Any] | None = None
    aliases: tuple[str, ...] = ()
    #: a separator follows this tool in its menu (ends a group)
    separator_after: bool = False

    def missing(self) -> list[str]:
        """Required modules that cannot be imported."""
        return [name for name in self.requires if importlib.util.find_spec(name) is None]


class ToolContext:
    """What a tool may do with the open document (and nothing else)."""

    def __init__(self, window: MainWindow) -> None:
        self._window = window
        self._finishing: Job | None = None  # the job whose results are being applied

    # -- reading ------------------------------------------------------------------------------
    @property
    def table(self) -> Table:
        return self._window.workspace.table

    @property
    def selection(self) -> Selection:
        return self._window.workspace.selection

    @property
    def lists(self) -> RowLists:
        """The document's named row lists."""
        return self._window.workspace.lists

    @property
    def parent(self) -> QWidget:
        """Parent for the tool's dialogs."""
        return self._window

    def structure_columns(self) -> list[str]:
        table = self.table
        return [name for name in table.column_names if table.kind(name) == STRUCTURE]

    def reaction_columns(self) -> list[str]:
        table = self.table
        return [name for name in table.column_names if table.kind(name) == REACTION]

    def default_reaction_column(self) -> str | None:
        """The current cell's column when it holds reactions, else the first reaction column."""
        current = self.current_column()
        names = self.reaction_columns()
        return current if current in names else (names[0] if names else None)

    def current_column(self) -> str | None:
        """The column of the table's current cell, if any."""
        index = self._window.view.currentIndex()
        return self._window.workspace.model.column_name(index.column()) if index.isValid() else None

    def current_oid(self) -> int | None:
        """The row of the table's current cell, if any."""
        index = self._window.view.currentIndex()
        return self._window.workspace.model.oid_at(index.row()) if index.isValid() else None

    def set_current_oid(self, oid: int) -> bool:
        """Make *oid* the table's current row (the selection is left alone), scrolled into
        view. False when the row is hidden by a filter or gone."""
        view, model = self._window.view, self._window.workspace.model
        row = int(model.rows_of_oids([oid])[0])
        if row < 0:
            return False
        current = view.currentIndex()
        target = model.index(row, current.column() if current.isValid() else 0)
        if target != current:
            view.selectionModel().setCurrentIndex(
                target, QItemSelectionModel.SelectionFlag.NoUpdate
            )
        view.scrollTo(target, QAbstractItemView.ScrollHint.EnsureVisible)
        return True

    def watch_current(self, callback: Callable[[int | None], None]) -> Callable[[], None]:
        """Call *callback* with the new current row whenever the table's current row changes;
        returns the function that stops it."""
        selection_model = self._window.view.selectionModel()

        def changed(current, _previous) -> None:
            model = self._window.workspace.model
            callback(model.oid_at(current.row()) if current.isValid() else None)

        selection_model.currentChanged.connect(changed)

        def stop() -> None:
            # RuntimeError: the document (and its model) closed first
            with contextlib.suppress(RuntimeError, TypeError):
                selection_model.currentChanged.disconnect(changed)

        return stop

    def default_structure_column(self) -> str | None:
        """The current cell's column when it holds structures, else the first structure column."""
        current = self.current_column()
        names = self.structure_columns()
        return current if current in names else (names[0] if names else None)

    def target_oids(self, selected_only: bool) -> np.ndarray:
        """The rows a tool works on: the selected rows, or every row (storage order)."""
        return self.selection.as_array() if selected_only else self.table.oids.copy()

    # -- work ---------------------------------------------------------------------------------
    def run(
        self,
        label: str,
        work: Callable[[Job], Any],
        on_done: Callable[[Any], None],
        *,
        lane: str = QUEUE,
        keep_partial: bool = False,
        background: bool = False,
    ) -> Job:
        """Run *work* in the background; *on_done* runs on the GUI thread when it finishes.

        Heavy work goes in the job queue (one job at a time, in order); pass ``lane=QUICK``
        for short interactive work that should never wait behind it. The window stays usable
        throughout. A message shown from *on_done* becomes the job's result line.

        With *keep_partial*, a cancelled job whose work raised ``Cancelled(partial=…)`` still
        calls *on_done* with that partial result, so finished rows are written, not thrown away.
        *on_done* must then accept a result covering only the first rows.

        Every tool that writes per-row results (or picks rows one at a time) keeps partial
        results this way: its core function raises ``Cancelled(partial)`` with what it
        finished (see ``compute_descriptors``, ``generate_conformers``, ``max_similarity``,
        ``decompose``, ``_maxmin``), and its dialog passes ``keep_partial=True``. Only work
        whose result means nothing until complete (a model fit, an embedding, clustering, a
        pair search) leaves it off.

        *background* marks work the user did not ask for (a panel making things ahead):
        closing the document cancels it without asking.
        """
        holder: list[Job] = []

        def finish(result: Any) -> None:
            self._finishing = holder[0]
            try:
                on_done(result)
            finally:
                self._finishing = None

        job = self._window._run_with_progress(
            label, work, finish, lane=lane, guard=True, keep_partial=keep_partial
        )
        job.background = background
        holder.append(job)
        return job

    def process_pool(self) -> Executor:
        """The shared worker processes, for CPU-heavy work inside a job."""
        return self._window._runner.process_pool()

    # -- writing ------------------------------------------------------------------------------
    def apply_results(
        self,
        results: Mapping[str, tuple[np.ndarray, Sequence[object] | Column]],
        *,
        label: str,
        overwrite: bool | Collection[str] = False,
    ) -> dict[str, str]:
        """Write result columns as one undoable step, skipping rows deleted meanwhile.

        Returns ``{requested name: actual name}``.
        """
        table = self.table
        kept: dict[str, tuple[np.ndarray, Sequence[object] | Column]] = {}
        for name, (oids, values) in results.items():
            oids = np.asarray(oids, dtype=np.int64)
            present = np.isin(oids, table.oids)
            if not present.all():  # the user deleted rows while the job ran
                index = np.flatnonzero(present)
                values = (
                    values.take(index)
                    if isinstance(values, Column)
                    else [values[i] for i in index.tolist()]
                )
                oids = oids[present]
            kept[name] = (oids, values)
        return table.apply_results(kept, label=label, overwrite=overwrite)

    # -- showing ------------------------------------------------------------------------------
    def add_panel(self, widget: QWidget, title: str) -> QDockWidget:
        """Show *widget* in a dock beside the table; it closes with the document.

        If the widget has ``close_panel()``, it is called when the dock is closed. A widget
        with a ``titleChanged(str)`` signal keeps the dock's title up to date.
        """
        return self._window.add_panel(widget, title)

    def find_panel(self, tool_id: str) -> QWidget | None:
        """The open panel of *tool_id*, if there is one (raised to the front).

        A tool that produces more results of the same kind (further retrosynthetic routes,
        say) fills the panel it already has instead of opening a second one.
        """
        for dock in self._window._panels:
            widget = dock.widget()
            if getattr(widget, "tool_id", None) == tool_id:
                dock.show()
                dock.raise_()
                return widget
        return None

    def show_dialog(self, dialog: QDialog) -> QDialog:
        """Show a tool's dialog without blocking the window.

        The table stays usable (select rows, change filters) while the dialog is open. The
        dialog closes with the document.
        """
        return self._window.show_dialog(dialog)

    def structure_image(
        self, column: str, oid: int, size: QSize, *, fit: bool = False
    ) -> QImage | None:
        """A drawing of one structure cell (cached; None for a blank or invalid cell). With
        *fit*, *size* is the largest box: a smaller molecule gets a smaller picture, drawn at
        the normal bond length so it is always legible."""
        table = self.table
        if column not in table or table.kind(column) != STRUCTURE or oid not in table.oids:
            return None
        window = self._window
        dark = window.palette().base().color().lightness() < 128
        slot = int(table.slots_of([oid])[0])
        structures = table.column(column)
        if fit:
            from chemsuite.core.depict import fitted_size

            mol = structures.obj(slot)
            if mol is None:
                return None
            size = QSize(*fitted_size(mol, size.width(), size.height()))
        return window._images.image(structures, slot, size, window.devicePixelRatioF(), dark)

    def add_compound(self, text: str) -> int | None:
        """Add a structure (SMILES, molfile, InChI) as a new, selected row, one undo step.

        It goes in the first structure column (a new "Structure" column if there is none).
        Returns the new row's oid, or None when *text* is blank.
        """
        return self._window.append_compound(text)

    def open_table(self, table: Table, message: str = "") -> Any:
        """Show *table* (e.g. an enumerated library) in a new, unsaved document window."""
        return self._window.open_table(table, message)

    def open_tool(self, tool_id: str) -> Any:
        """Open another tool (e.g. a scatter plot from a correlation) and return what it opened."""
        return self._window.open_tool(tool_id)

    def show_message(self, text: str, ms: int = 5000) -> None:
        self._window.statusBar().showMessage(text, ms)
        if self._finishing is not None:
            self._finishing.set_summary(text)


def load_tools(modules: Sequence[str]) -> list[Tool]:
    """Import each module and collect its ``TOOL`` (a module may also define ``TOOLS``)."""
    tools: list[Tool] = []
    for name in modules:
        module = importlib.import_module(name)
        tools.extend(getattr(module, "TOOLS", ()))
        if hasattr(module, "TOOL"):
            tools.append(module.TOOL)
    return tools
