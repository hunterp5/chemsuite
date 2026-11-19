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
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QLabel,
    QSpinBox,
    QVBoxLayout,
)

from chemsuite.analysis.som import (
    END_RADIUS,
    INITS,
    MAX_SIDE,
    NEIGHBORHOODS,
    Som,
    default_size,
    som,
)
from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.widgets.dialogs import FitChoice, RowsChoice, Segmented, ToolDialog, inline
from chemsuite.app.widgets.flow_layout import FlowRow, compact, field
from chemsuite.core.table import Changed
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural
from chemsuite.plotting.data import number_columns, numbers, shown_oids
from chemsuite.sar.matrix import AGGREGATES
from chemsuite.tools.sar.common import activity_combo
from chemsuite.tools.sar.matrix_panel import heat
from chemsuite.tools.similarity.common import FeatureChooser

from .common import features_task
from .som_grid import Cell, SomGrid

ROW_COLUMN = "SOM Row"
COLUMN_COLUMN = "SOM Column"
COUNT = "Count"
DELAY_MS = 300


class SomDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Self-Organizing Map", guide="self-organizing-map")
        self.ctx = ctx
        self.panel: SomPanel | None = None
        inputs = self.section("Input")
        method = self.section("Method")
        more = self.section("More Options", collapsed=False)
        self._inputs, self._more_section = inputs, more
        self.features = FeatureChooser(ctx, metric=False, section=inputs, more=more)
        self.rows = inputs.row("Rows", RowsChoice(ctx))
        self.fit = inputs.row("Fit on", FitChoice(ctx),
                              help="The rows the map is trained on; the other rows are then "
                              "placed in their closest cell.")  # fmt: skip

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
        self.cells = Segmented((("square", "Square"), ("hexagonal", "Hexagonal")), "square")
        method.row("Cells", self.cells,
                   help="Hexagonal cells have six neighbors at the same distance, so groups "
                   "spread more evenly than on squares.")  # fmt: skip
        self.epochs = QSpinBox()
        self.epochs.setRange(2, 200)
        self.epochs.setValue(20)
        method.row("Training passes", self.epochs, help="Passes over the data while training.")

        self.neighborhood = QComboBox()
        for key, title in NEIGHBORHOODS.items():
            self.neighborhood.addItem(title, key)
        more.row("Neighborhood", self.neighborhood,
                 help="How a compound pulls on the cells around its best cell. Gaussian: less "
                 "with distance. Bubble: every cell within the radius alike.")  # fmt: skip
        self.start_radius = QDoubleSpinBox()
        self.start_radius.setRange(0.0, float(MAX_SIDE))
        self.start_radius.setSingleStep(0.5)
        self.start_radius.setSpecialValueText("Automatic")
        more.row("Start radius", self.start_radius,
                 help="How far the pull reaches at the start, in cells. Automatic: half the "
                 "grid.")  # fmt: skip
        self.end_radius = QDoubleSpinBox()
        self.end_radius.setRange(0.1, float(MAX_SIDE))
        self.end_radius.setSingleStep(0.1)
        self.end_radius.setValue(END_RADIUS)
        more.row("End radius", self.end_radius,
                 help="How far it reaches at the end: small lets each cell settle on its own "
                 "compounds.")  # fmt: skip
        self.init = QComboBox()
        for key, title in INITS.items():
            self.init.addItem(title, key)
        more.row("Start from", self.init,
                 help="The cells' first values: spread over the plane of the first two "
                 "principal components, or copied from random compounds.")  # fmt: skip
        self.wrap = QCheckBox("Join opposite edges")
        more.row("Wrap edges", self.wrap,
                 help="The grid wraps around like a torus, so edge cells have neighbors on "
                 "the far side and nothing crowds into the corners.")  # fmt: skip
        self.seed = QSpinBox()
        self.seed.setRange(0, 2_147_483_647)
        more.row("Seed", self.seed, help="The same seed gives the same map every run.")

        for signal in (self.auto.toggled, self.rows.changed, self.fit.currentIndexChanged,
                       self.start_radius.valueChanged, self.end_radius.valueChanged):  # fmt: skip
            signal.connect(lambda *_: self._update())
        self.features.changed.connect(self._update)
        self._update()
        self.remember(ctx, "som", {
            "auto": self.auto, "grid_rows": self.grid_rows, "grid_columns": self.grid_columns,
            "cells": self.cells, "epochs": self.epochs, "fit": self.fit,
            "neighborhood": self.neighborhood, "start_radius": self.start_radius,
            "end_radius": self.end_radius, "init": self.init, "wrap": self.wrap,
            "seed": self.seed, **self.features.remembered(),
        })  # fmt: skip

    def size(self) -> tuple[int, int] | None:
        if self.auto.isChecked():
            return None
        return self.grid_rows.value(), self.grid_columns.value()

    def _update(self) -> None:
        if self.auto.isChecked():
            oids = self.rows.oids()
            fit = self.fit.rows(oids)
            rows, columns = default_size(len(oids) if fit is None else len(fit))
            self.grid_rows.setValue(rows)
            self.grid_columns.setValue(columns)
        for spin in (self.grid_rows, self.grid_columns):
            spin.setEnabled(not self.auto.isChecked())
        self.refresh()

    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        fit = self.fit.rows(self.rows.oids())
        if fit is not None and len(fit) < 2:
            return "Train on at least 2 of the rows mapped."
        start = self.start_radius.value()
        if start and start < self.end_radius.value():
            return "The start radius must not be smaller than the end radius."
        return ""

    def verb(self) -> str:
        return f"Train on {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.train()

    def options(self) -> dict:
        """What ``som`` takes besides the features and the grid size."""
        oids = self.rows.oids()
        return {
            "epochs": self.epochs.value(),
            "fit": self.fit.rows(oids),
            "hexagonal": self.cells.value() == "hexagonal",
            "wrap": self.wrap.isChecked(),
            "neighborhood": self.neighborhood.currentData(),
            "start_radius": self.start_radius.value() or None,
            "end_radius": self.end_radius.value(),
            "init": self.init.currentData(),
            "seed": self.seed.value(),
        }

    def train(self):
        ctx = self.ctx
        oids = self.rows.oids()
        build, remember = features_task(ctx, self.features, oids)
        size = self.size() or (None, None)
        options = self.options()

        def work(job) -> Som:
            features = build(job)
            return som(features, *size, progress=job.report, cancel=job.cancel_token, **options)

        def done(result: Som) -> None:
            remember()
            self.finish(oids, result)

        return ctx.run("Calculating Self-Organizing Map…", work, done)

    def finish(self, oids: np.ndarray, result: Som):
        ctx = self.ctx
        r, c = result.cell()
        ctx.apply_results(
            {
                ROW_COLUMN: (oids, NumberColumn.from_floats(r, 0)),
                COLUMN_COLUMN: (oids, NumberColumn.from_floats(c, 0)),
            },
            label="Self-organizing map",
        )
        counts = result.counts()
        from .tool import _som_panel

        layout = {"hexagonal": result.hexagonal, "wrap": result.wrap}
        layout |= {"rows": result.rows, "columns": result.columns}
        self.panel = _som_panel(ctx, layout)
        ctx.show_message(
            f"Placed {plural(int(counts.sum()), 'compound')} on a {result.rows} × "
            f"{result.columns} map; {plural(int((counts == 0).sum()), 'empty cell')}"
        )
        return self.panel


class SomPanel(ToolPanel):
    """The map's grid: each cell's count, or the mean (etc.) of a number column."""

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
        self.hexagonal = False  # how the map was trained: drawn the same way
        self.wrap = False
        self.trained = (0, 0)  # the trained grid, which may have empty edge rows
        self.grid = SomGrid()
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
            grid.clear()
            self.status.setText(
                "Train a map first (Data ▸ Dimensionality Reduction ▸ Self-Organizing Map)."
            )
            return
        oids = shown_oids(table)
        r = numbers(table, ROW_COLUMN, oids)
        c = numbers(table, COLUMN_COLUMN, oids)
        ok = np.isfinite(r) & np.isfinite(c)
        oids, r, c = oids[ok], r[ok].astype(int) - 1, c[ok].astype(int) - 1
        rows = max(self.trained[0], int(r.max()) + 1 if len(r) else 0)
        columns = max(self.trained[1], int(c.max()) + 1 if len(c) else 0)
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
        cells = []
        for k in range(rows * columns):
            members = order[bounds[k] : bounds[k + 1]]
            if not len(members):
                cells.append(Cell())
                continue
            value = shown[k]
            text = (
                str(len(members)) if counting else (f"{value:.2f}" if math.isfinite(value) else "–")
            )
            tip = f"Cell {k // columns + 1}, {k % columns + 1}: {plural(len(members), 'compound')}"
            if not counting:
                tip += f"\n{self.value.currentText()}: {text}"
            cells.append(Cell(text, heat(value, low, high, better), oids[members].tolist(), tip))
        grid.set_cells(rows, columns, cells, self.hexagonal)
        filled = int((bounds[1:] > bounds[:-1]).sum())
        self.status.setText(
            f"{len(oids):,} compound(s) in {filled:,} of {rows * columns:,} cell(s). "
            "Click a cell to select its compounds."
        )

    def select_cell(self, row: int, column: int) -> list[int]:
        cell = self.grid.cell(row, column)
        oids = cell.oids if cell is not None else []
        if oids:
            self.ctx.selection.set(oids, "som")
        return list(oids)

    def state(self) -> dict:
        return {
            "value": self.value.currentText(),
            "aggregate": self.aggregate.currentData(),
            "higher_is_better": self.higher_is_better.isChecked(),
            "hexagonal": self.hexagonal,
            "wrap": self.wrap,
            "rows": self.trained[0],
            "columns": self.trained[1],
        }

    def set_state(self, state: dict) -> None:
        index = self.value.findText(str(state.get("value", "")))
        if index >= 0:
            self.value.setCurrentIndex(index)
        self.aggregate.setCurrentIndex(max(0, self.aggregate.findData(state.get("aggregate"))))
        self.higher_is_better.setChecked(bool(state.get("higher_is_better", True)))
        self.hexagonal = state.get("hexagonal") is True  # older documents: square cells
        self.wrap = state.get("wrap") is True
        size = state.get("rows"), state.get("columns")
        ok = all(isinstance(n, int) and 0 < n <= MAX_SIDE for n in size)
        self.trained = (size[0], size[1]) if ok else (0, 0)
        self.rebuild()

    def close_panel(self) -> None:
        self._timer.stop()
        self._unsubscribe()
