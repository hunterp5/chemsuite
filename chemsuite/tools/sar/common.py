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
"""Widgets and helpers shared by the SAR tools."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.selectivity import looks_like_assay
from chemsuite.app.tools import ToolContext
from chemsuite.core.mols import mol_from_blob, mol_from_text
from chemsuite.plotting.data import number_columns
from chemsuite.sar.rgroups import Decomposition, decompose, parse_core, suggest_core

MCS_SAMPLE = 50  # structures used to suggest a core


def activity_combo(ctx: ToolContext) -> QComboBox:
    combo = QComboBox()
    names = number_columns(ctx.table)
    combo.addItems(names)
    first = next((n for n in names if looks_like_assay(n)), None)
    if first:
        combo.setCurrentText(first)
    return combo


def structure_cells(ctx: ToolContext, column: str, oids: np.ndarray) -> list:
    return ctx.structure_inputs(column, oids)


def saved_rows(ctx: ToolContext, state: dict, *columns: str) -> list[int]:
    """The rows a saved analysis covered that are still in the table. Raises KeyError when
    a column it used has gone (the panel is then not reopened)."""
    for name in columns:
        if name not in ctx.table:
            raise KeyError(name)
    oids = np.asarray(state.get("oids", []), dtype=np.int64)
    return oids[np.isin(oids, ctx.table.oids)].tolist()


class RerunPanel(QWidget):
    """Stands in for an analysis panel reopened with its document while the analysis reruns,
    then holds the result. The dock is in place at once, so the saved arrangement applies."""

    titleChanged = Signal(str)

    def __init__(self, tool_id: str, title: str, state: dict) -> None:
        super().__init__()
        self.tool_id = tool_id
        self._title = title
        self._state = state
        self.panel: QWidget | None = None
        self.job = None
        self.waiting = QLabel("Recalculating…")
        self.waiting.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.waiting.setEnabled(False)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.addWidget(self.waiting)

    def title(self) -> str:
        return self.panel.title() if self.panel is not None else self._title

    def fill(self, panel: QWidget) -> None:
        self.panel = panel
        self.job = None
        self.waiting.hide()
        self._layout.addWidget(panel)
        if hasattr(panel, "titleChanged"):
            panel.titleChanged.connect(self.titleChanged)
        self.titleChanged.emit(panel.title())

    def state(self) -> dict:
        return self.panel.state() if self.panel is not None else self._state

    def close_panel(self) -> None:
        if self.job is not None and self.job.running:
            self.job.cancel()
        if self.panel is not None and hasattr(self.panel, "close_panel"):
            self.panel.close_panel()


def cell_mol(cell):
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


class CoreChooser(QWidget):
    """A structure column and a core (typed, or suggested from the selected compounds).

    Its rows go into a dialog *section* (see ``chemsuite.app.dialogs``).
    """

    def __init__(self, ctx: ToolContext, section) -> None:
        from chemsuite.app.dialogs import inline

        super().__init__()
        self.ctx = ctx
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        self.core = QLineEdit()
        self.core.setPlaceholderText("SMILES or SMARTS")
        suggest = QPushButton("Suggest")
        suggest.setAutoDefault(False)
        suggest.setToolTip(
            "The largest substructure the selected compounds share (MCS); with fewer than two "
            "selected, the first rows shown."
        )
        suggest.clicked.connect(self.suggest)
        self.exhaustive = QCheckBox("Exhaustive")
        section.row("Structures", self.structure)
        core_help = (
            "The common core, for example [*:1]c1ccc(C(N)=O)cc1[*:2]. Number the attachment "
            "points [*:1], [*:2]… to fix which position is R1, R2…; without numbers, every "
            "open position on the core becomes an R-group."
        )
        section.row("Core", inline(self.core, suggest, stretch=0), help=core_help)
        self.core.setToolTip(core_help)
        section.row(
            "Matching",
            self.exhaustive,
            help="Exhaustive matching tries every way the core fits: slower, for symmetric cores.",
        )
        self.message = ""

    def suggest(self) -> str:
        ctx = self.ctx
        oids = ctx.selection.as_array()
        if len(oids) < 2:
            oids = ctx.table.oids[ctx.table.view()][:MCS_SAMPLE]
        oids = oids[:MCS_SAMPLE]
        column = self.structure.currentText()
        mols = [cell_mol(c) for c in structure_cells(ctx, column, oids)] if column else []
        try:
            smarts = suggest_core(mols)
        except ValueError as err:
            self.message = str(err)
            ctx.show_message(self.message)
            return ""
        self.core.setText(smarts)
        self.message = ""
        return smarts

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column."
        try:
            parse_core(self.core.text())
        except ValueError as err:
            return str(err)
        return ""

    def recipe(self) -> dict:
        """The column, core and matching, as saved with an analysis panel."""
        return {
            "structure": self.structure.currentText(),
            "core": self.core.text(),
            "exhaustive": self.exhaustive.isChecked(),
        }

    def task(self, oids: np.ndarray):
        """A function for a job: ``work(job) -> Decomposition`` for *oids*."""
        return decompose_task(self.ctx, self.recipe(), oids)


def decompose_task(ctx: ToolContext, recipe: dict, oids: np.ndarray):
    """``work(job) -> Decomposition`` of *oids* onto the core *recipe* names
    (see :meth:`CoreChooser.recipe`)."""
    core = parse_core(recipe["core"])
    cells = structure_cells(ctx, recipe["structure"], oids)
    exhaustive = bool(recipe["exhaustive"])

    def work(job) -> Decomposition:
        mols = [cell_mol(c) for c in cells]
        return decompose(
            mols, core, exhaustive=exhaustive, progress=job.report, cancel=job.cancel_token
        )

    return work


def number_text(value: float, digits: int = 3) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.{digits}f}"


class ResultTable(QTableWidget):
    """A read-only, sortable table of results."""

    def __init__(self, headers: Sequence[str]) -> None:
        super().__init__(0, len(headers))
        self.setHorizontalHeaderLabels(list(headers))
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.verticalHeader().hide()
        self.verticalHeader().setDefaultSectionSize(22)
        self.horizontalHeader().setStretchLastSection(False)

    def fill(self, rows: Sequence[Sequence[object]], data: Sequence[object] | None = None) -> None:
        """One list of cells per row; numbers align right and sort as numbers.

        *data* attaches a value per row (read back with :meth:`row_data`), so sorting by a
        column does not lose track of which result a row shows.
        """
        self.setSortingEnabled(False)
        self.setRowCount(len(rows))
        for r, cells in enumerate(rows):
            for c, value in enumerate(cells):
                item = QTableWidgetItem()
                if isinstance(value, int | float | np.integer | np.floating):
                    number = float(value)
                    item.setData(
                        Qt.ItemDataRole.DisplayRole, number if math.isfinite(number) else None
                    )
                    if isinstance(value, float | np.floating) and math.isfinite(number):
                        item.setData(Qt.ItemDataRole.DisplayRole, round(number, 3))
                    item.setTextAlignment(
                        Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                    )
                else:
                    item.setText(str(value))
                if c == 0 and data is not None:
                    item.setData(Qt.ItemDataRole.UserRole, r)
                self.setItem(r, c, item)
        self._data = list(data) if data is not None else None
        self.resizeColumnsToContents()
        for c in range(self.columnCount()):  # room for the sort arrow; long text stays readable
            self.setColumnWidth(c, min(self.columnWidth(c) + 18, 360))
        # Keep the given order (e.g. best first) until the user clicks a header.
        self.horizontalHeader().setSortIndicator(-1, Qt.SortOrder.AscendingOrder)
        self.setSortingEnabled(True)

    def row_data(self, row: int):
        if self._data is None or row < 0:
            return None
        return self._data[self.item(row, 0).data(Qt.ItemDataRole.UserRole)]

    def selected_data(self) -> list:
        rows = sorted({i.row() for i in self.selectedIndexes()})
        return [self.row_data(r) for r in rows]

    def copy_text(self) -> str:
        headers = [self.horizontalHeaderItem(c).text() for c in range(self.columnCount())]
        lines = ["\t".join(headers)]
        for r in range(self.rowCount()):
            cells = [
                self.item(r, c).text() if self.item(r, c) else "" for c in range(self.columnCount())
            ]
            lines.append("\t".join(cells))
        text = "\n".join(lines)
        QGuiApplication.clipboard().setText(text)
        return text
