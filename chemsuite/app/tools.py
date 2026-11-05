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

import importlib
import importlib.util
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from PySide6.QtWidgets import QWidget

from chemsuite.core.selection import Selection
from chemsuite.core.table import Table
from chemsuite.core.table.columns import STRUCTURE, Column

from .jobs import Job

if TYPE_CHECKING:
    from .main_window import MainWindow


@dataclass(frozen=True)
class Tool:
    """A menu entry and what it opens.

    ``menu`` is the menu path, e.g. ``("Chemistry",)`` or ``("Chemistry", "Calculate")``.
    ``requires`` names importable modules; a tool whose module is missing is shown disabled
    with the reason. ``needs_structures`` disables the tool while the table has no structure
    column.
    """

    id: str
    title: str
    menu: tuple[str, ...]
    open: Callable[[ToolContext], Any]
    requires: tuple[str, ...] = ()
    needs_structures: bool = False
    shortcut: str | None = None
    tip: str = ""

    def missing(self) -> list[str]:
        """Required modules that cannot be imported."""
        return [name for name in self.requires if importlib.util.find_spec(name) is None]


class ToolContext:
    """What a tool may do with the open document (and nothing else)."""

    def __init__(self, window: MainWindow) -> None:
        self._window = window

    # -- reading ------------------------------------------------------------------------------
    @property
    def table(self) -> Table:
        return self._window.workspace.table

    @property
    def selection(self) -> Selection:
        return self._window.workspace.selection

    @property
    def parent(self) -> QWidget:
        """Parent for the tool's dialogs."""
        return self._window

    def structure_columns(self) -> list[str]:
        table = self.table
        return [name for name in table.column_names if table.kind(name) == STRUCTURE]

    def current_column(self) -> str | None:
        """The column of the table's current cell, if any."""
        index = self._window.view.currentIndex()
        return self._window.workspace.model.column_name(index.column()) if index.isValid() else None

    def default_structure_column(self) -> str | None:
        """The current cell's column when it holds structures, else the first structure column."""
        current = self.current_column()
        names = self.structure_columns()
        return current if current in names else (names[0] if names else None)

    def target_oids(self, selected_only: bool) -> np.ndarray:
        """The rows a tool works on: the selected rows, or every row (storage order)."""
        return self.selection.as_array() if selected_only else self.table.oids.copy()

    # -- work ---------------------------------------------------------------------------------
    def run(self, label: str, work: Callable[[Job], Any], on_done: Callable[[Any], None]) -> Job:
        """Run *work* in the background with a progress dialog; *on_done* runs on the GUI thread."""
        return self._window._run_with_progress(label, work, on_done)

    def process_pool(self) -> Executor:
        """The shared worker processes, for CPU-heavy work inside a job."""
        return self._window._runner.process_pool()

    # -- writing ------------------------------------------------------------------------------
    def apply_results(
        self,
        results: Mapping[str, tuple[np.ndarray, Sequence[object] | Column]],
        *,
        label: str,
        overwrite: bool = False,
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

    def show_message(self, text: str, ms: int = 5000) -> None:
        self._window.statusBar().showMessage(text, ms)


def load_tools(modules: Sequence[str]) -> list[Tool]:
    """Import each module and collect its ``TOOL`` (a module may also define ``TOOLS``)."""
    tools: list[Tool] = []
    for name in modules:
        module = importlib.import_module(name)
        tools.extend(getattr(module, "TOOLS", ()))
        if hasattr(module, "TOOL"):
            tools.append(module.TOOL)
    return tools
