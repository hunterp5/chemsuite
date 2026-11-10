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
"""Self-Organizing Map…: train a map, write each compound's cell, and show the grid."""

from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QHeaderView,
    QLabel,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.som import MAX_SIDE, Som, default_size, som
from chemsuite.app.dialogs import RowsChoice, ToolDialog, inline
from chemsuite.app.flow_layout import FlowRow, compact, field
from chemsuite.app.tools import ToolContext
from chemsuite.core.table import Changed
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import number_columns, numbers, shown_oids
from chemsuite.sar.matrix import AGGREGATES
from chemsuite.tools.sar.common import activity_combo
from chemsuite.tools.sar.matrix_panel import heat
from chemsuite.tools.similarity.common import FeatureChooser

from .common import features_task

ROW_COLUMN = "SOM Row"
COLUMN_COLUMN = "SOM Column"
COUNT = "Count"
DELAY_MS = 300
OIDS_ROLE = Qt.ItemDataRole.UserRole


class SomDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Self-Organizing Map", guide="self-organizing-map")
        self.ctx = ctx
        self.panel: SomPanel | None = None
        inputs = self.section("Input")
        self.features = FeatureChooser(ctx, metric=False, section=inputs)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.grid_rows = QSpinBox()
        self.grid_columns = QSpinBox()
        for spin in (self.grid_rows, self.grid_columns):
            spin.setRange(2, MAX_SIDE)
        self.auto = QCheckBox("Automatic")
        self.auto.setChecked(True)
        method.row(
            "Grid",
            inline(self.grid_rows, QLabel("×"), self.grid_columns, self.auto),
            help="Rows × columns of cells. Automatic: a size that suits the number of compounds.",
        )
        self.epochs = QSpinBox()
        self.epochs.setRange(2, 200)
        self.epochs.setValue(20)
        method.row("Training passes", self.epochs, help="Passes over the data while training.")

        self.auto.toggled.connect(lambda *_: self._update())
        self.rows.changed.connect(lambda *_: self._update())
        self.features.changed.connect(self._update)
        self._update()

    def size(self) -> tuple[int, int] | None:
        if self.auto.isChecked():
            return None
        return self.grid_rows.value(), self.grid_columns.value()

    def _update(self) -> None:
        if self.auto.isChecked():
            rows, columns = default_size(len(self.rows.oids()))
            self.grid_rows.setValue(rows)
            self.grid_columns.setValue(columns)
        for spin in (self.grid_rows, self.grid_columns):
            spin.setEnabled(not self.auto.isChecked())
        self.refresh()

    def problem(self) -> str:
        return self.features.problem()

    def verb(self) -> str:
        return f"Train on {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.train()

    def train(self):
        ctx = self.ctx
        oids = self.rows.oids()
        build, remember = features_task(ctx, self.features, oids)
        size = self.size() or (None, None)
        epochs = self.epochs.value()

        def work(job) -> Som:
            features = build(job)
            return som(features, *size, epochs=epochs, progress=job.report, cancel=job.cancel_token)

        def done(result: Som) -> None:
            remember()
            self.finish(oids, result)

        return ctx.run("Training a self-organising map…", work, done)

    def finish(self, oids: np.ndarray, result: Som):
        ctx = self.ctx
        r, c = result.cell()
        ctx.apply_results(
            {
                ROW_COLUMN: (oids, NumberColumn.from_floats(r, 0)),
                COLUMN_COLUMN: (oids, NumberColumn.from_floats(c, 0)),
            },
            label="Self-organising map",
        )
        counts = result.counts()
        from .tool import _som_panel

        self.panel = _som_panel(ctx)
        ctx.show_message(
            f"{int(counts.sum()):,} compound(s) on a {result.rows} × {result.columns} map; "
            f"{int((counts == 0).sum()):,} empty cell(s)",
            10_000,
        )
        return self.panel


class SomPanel(QWidget):
    """The map's grid: each cell's count, or the mean (etc.) of a number column."""

    titleChanged = Signal(str)
    tool_id = "som"

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__()
        self.ctx = ctx
        self.table = ctx.table
        self.value = QComboBox()
        self.aggregate = QComboBox()
        for key, title in AGGREGATES.items():
            self.aggregate.addItem(title, key)
        self.higher_is_better = QCheckBox("Higher is better")
        self.higher_is_better.setChecked(True)
        self.grid = QTableWidget()
        self.grid.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        for header in (self.grid.horizontalHeader(), self.grid.verticalHeader()):
            header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            header.setMinimumSectionSize(18)
        self.status = QLabel()
        self.status.setWordWrap(True)
        compact(self.value)
        controls = FlowRow()
        controls.addWidget(field("Show:", self.value))
        controls.addWidget(self.aggregate)
        controls.addWidget(self.higher_is_better)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(controls)
        layout.addWidget(self.grid, 1)
        layout.addWidget(self.status)
        self._fill_values(initial=True)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(DELAY_MS)
        self._timer.timeout.connect(self.rebuild)
        self.value.currentIndexChanged.connect(lambda *_: self.rebuild())
        self.aggregate.currentIndexChanged.connect(lambda *_: self.rebuild())
        self.higher_is_better.toggled.connect(lambda *_: self.rebuild())
        self.grid.cellClicked.connect(self.select_cell)
        self._unsubscribe = self.table.subscribe(self._on_table_event)
        self.rebuild()

    def title(self) -> str:
        return "Self-Organizing Map"

    def _fill_values(self, *, initial: bool = False) -> None:
        current = self.value.currentText()
        names = [n for n in number_columns(self.table) if n not in (ROW_COLUMN, COLUMN_COLUMN)]
        self.value.blockSignals(True)
        self.value.clear()
        self.value.addItems([COUNT, *names])
        if current:
            self.value.setCurrentText(current)
        elif initial:
            guess = activity_combo(self.ctx).currentText()
            if guess in names:
                self.value.setCurrentText(guess)
        self.value.blockSignals(False)

    def _on_table_event(self, event: object) -> None:
        if isinstance(event, Changed):
            if event.structural:
                self._fill_values()
            self._timer.start()

    def rebuild(self) -> None:
        self._timer.stop()
        table = self.table
        grid = self.grid
        if ROW_COLUMN not in table or COLUMN_COLUMN not in table:
            grid.setRowCount(0)
            grid.setColumnCount(0)
            self.status.setText("Train a map first (Data ▸ Chemical Space ▸ Self-Organizing Map).")
            return
        oids = shown_oids(table)
        r = numbers(table, ROW_COLUMN, oids)
        c = numbers(table, COLUMN_COLUMN, oids)
        ok = np.isfinite(r) & np.isfinite(c)
        oids, r, c = oids[ok], r[ok].astype(int) - 1, c[ok].astype(int) - 1
        rows = int(r.max()) + 1 if len(r) else 0
        columns = int(c.max()) + 1 if len(c) else 0
        counting = self.value.currentText() == COUNT
        self.aggregate.setEnabled(not counting)
        values = None if counting else numbers(table, self.value.currentText(), oids)
        aggregate = self.aggregate.currentData()
        cell = r * max(columns, 1) + c
        order = np.argsort(cell, kind="stable")
        bounds = np.searchsorted(cell[order], np.arange(rows * columns + 1))
        shown = np.full(rows * columns, np.nan)
        for k in range(rows * columns):
            members = order[bounds[k] : bounds[k + 1]]
            if not len(members):
                continue
            if counting:
                shown[k] = len(members)
            else:
                v = values[members]
                v = v[np.isfinite(v)]
                if len(v):
                    shown[k] = float(getattr(np, aggregate)(v))
        finite = shown[np.isfinite(shown)]
        low, high = (float(finite.min()), float(finite.max())) if len(finite) else (0.0, 0.0)
        better = self.higher_is_better.isChecked() or counting
        grid.clear()
        grid.setRowCount(rows)
        grid.setColumnCount(columns)
        grid.setHorizontalHeaderLabels([str(i + 1) for i in range(columns)])
        grid.setVerticalHeaderLabels([str(i + 1) for i in range(rows)])
        for k in range(rows * columns):
            members = order[bounds[k] : bounds[k + 1]]
            item = QTableWidgetItem()
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            if len(members):
                value = shown[k]
                text = (
                    str(len(members))
                    if counting
                    else (f"{value:.2f}" if math.isfinite(value) else "–")
                )
                item.setText(text)
                item.setBackground(QBrush(heat(value, low, high, better)))
                item.setData(OIDS_ROLE, oids[members].tolist())
                item.setToolTip(
                    f"Cell {k // columns + 1}, {k % columns + 1}: {len(members)} compound(s)"
                )
            else:
                item.setBackground(QBrush(QColor(235, 235, 235)))
            grid.setItem(k // columns, k % columns, item)
        filled = int((bounds[1:] > bounds[:-1]).sum())
        self.status.setText(
            f"{len(oids):,} compound(s) in {filled:,} of {rows * columns:,} cell(s). "
            "Click a cell to select its compounds."
        )

    def select_cell(self, row: int, column: int) -> list[int]:
        item = self.grid.item(row, column)
        oids = item.data(OIDS_ROLE) if item is not None else None
        if oids:
            self.ctx.selection.set(oids, "som")
        return oids or []

    def state(self) -> dict:
        return {
            "value": self.value.currentText(),
            "aggregate": self.aggregate.currentData(),
            "higher_is_better": self.higher_is_better.isChecked(),
        }

    def set_state(self, state: dict) -> None:
        index = self.value.findText(str(state.get("value", "")))
        if index >= 0:
            self.value.setCurrentIndex(index)
        self.aggregate.setCurrentIndex(max(0, self.aggregate.findData(state.get("aggregate"))))
        self.higher_is_better.setChecked(bool(state.get("higher_is_better", True)))
        self.rebuild()

    def close_panel(self) -> None:
        self._timer.stop()
        self._unsubscribe()
