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
"""SAR Matrix: activity by substituent pairs, with Free-Wilson estimates for the gaps."""

from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QBrush, QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.flow_layout import FlowRow, compact, field
from chemsuite.app.tools import ToolContext
from chemsuite.core.table import Changed
from chemsuite.plotting.data import number_columns, numbers, shown_oids
from chemsuite.sar.matrix import AGGREGATES, SarMatrix, build_matrix

from .common import activity_combo, number_text

DELAY_MS = 300
NONE = "(none)"
LOW = QColor(214, 96, 77)
MID = QColor(247, 247, 247)
HIGH = QColor(77, 172, 120)
OIDS_ROLE = Qt.ItemDataRole.UserRole


def label_columns(ctx: ToolContext) -> list[str]:
    """Text columns: R-group labels and other categories."""
    table = ctx.table
    numeric = set(number_columns(table))
    return [n for n in table.column_names if n not in numeric and table.kind(n) != "structure"]


def heat(value: float, low: float, high: float, higher_is_better: bool = True) -> QColor:
    """Red (worst) through white to green (best)."""
    if not math.isfinite(value) or high <= low:
        return MID
    t = (value - low) / (high - low) * 2 - 1
    if not higher_is_better:
        t = -t
    end = HIGH if t > 0 else LOW
    t = min(1.0, abs(t))
    mix = [round(a + (b - a) * t) for a, b in zip(MID.getRgb()[:3], end.getRgb()[:3], strict=True)]
    return QColor(*mix)


class SarMatrixPanel(QWidget):
    titleChanged = Signal(str)
    tool_id = "sar_matrix"

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__()
        self.ctx = ctx
        self.table = ctx.table
        self.matrix: SarMatrix | None = None
        self.oids = np.empty(0, dtype=np.int64)
        self.rows_axis = QComboBox()
        self.columns_axis = QComboBox()
        self.pages_axis = QComboBox()
        self.activity = activity_combo(ctx)
        self.aggregate = QComboBox()
        for key, title in AGGREGATES.items():
            self.aggregate.addItem(title, key)
        self.page = QComboBox()
        self.predict = QCheckBox("Estimate empty cells (Free-Wilson)")
        self.predict.setChecked(True)
        self.higher_is_better = QCheckBox("Higher is better")
        self.higher_is_better.setChecked(True)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.grid = QTableWidget()
        self.grid.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.grid.horizontalHeader().setDefaultSectionSize(76)
        self.grid.verticalHeader().setDefaultSectionSize(26)
        copy = QPushButton("Copy")
        copy.setToolTip("Copy the matrix as tab-separated text, for Excel.")
        copy.clicked.connect(self.copy_text)

        compact(self.rows_axis, self.columns_axis, self.pages_axis, self.activity, self.page)
        controls = FlowRow()
        for title, widget in (
            ("Rows:", self.rows_axis),
            ("Columns:", self.columns_axis),
            ("Pages:", self.pages_axis),
            ("Activity:", self.activity),
            ("Cell value:", self.aggregate),
            ("Page:", self.page),
        ):
            controls.addWidget(field(title, widget))
        options = FlowRow()
        options.addWidget(self.predict)
        options.addWidget(self.higher_is_better)
        options.addWidget(copy)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(controls)
        layout.addWidget(options)
        layout.addWidget(self.grid, 1)
        layout.addWidget(self.status)

        self._fill_axes()
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(DELAY_MS)
        self._timer.timeout.connect(self.rebuild)
        for combo in (self.rows_axis, self.columns_axis, self.pages_axis, self.activity):
            combo.currentIndexChanged.connect(lambda *_: self.invalidate())
        self.predict.toggled.connect(lambda *_: self.invalidate())
        self.aggregate.currentIndexChanged.connect(lambda *_: self.show_page())
        self.higher_is_better.toggled.connect(lambda *_: self.show_page())
        self.page.currentIndexChanged.connect(lambda *_: self.show_page())
        self.grid.cellClicked.connect(self.select_cell)
        self._unsubscribe = self.table.subscribe(self._on_table_event)
        self.rebuild()

    def title(self) -> str:
        return "SAR Matrix"

    def _fill_axes(self) -> None:
        names = label_columns(self.ctx)
        for combo, default, blank in (
            (self.rows_axis, "R1", False),
            (self.columns_axis, "R2", False),
            (self.pages_axis, NONE, True),
        ):
            current = combo.currentText() or default
            combo.blockSignals(True)
            combo.clear()
            if blank:
                combo.addItem(NONE)
            combo.addItems(names)
            index = combo.findText(current)
            if index < 0 and not blank and names:
                index = min(len(names) - 1, 0 if combo is self.rows_axis else 1)
            combo.setCurrentIndex(max(0, index))
            combo.blockSignals(False)
        current = self.activity.currentText()
        self.activity.blockSignals(True)
        self.activity.clear()
        self.activity.addItems(number_columns(self.table))
        if current:
            self.activity.setCurrentText(current)
        self.activity.blockSignals(False)

    def _on_table_event(self, event: object) -> None:
        if isinstance(event, Changed):
            if event.structural:
                self._fill_axes()
            self.invalidate()

    def invalidate(self) -> None:
        self._timer.start()

    def problem(self) -> str:
        names = (self.rows_axis.currentText(), self.columns_axis.currentText())
        if not all(names):
            return "Choose two text columns, such as R1 and R2 from R-Group Decomposition."
        if names[0] == names[1]:
            return "Choose different columns for rows and columns."
        if not self.activity.currentText():
            return "There is no number column for activity."
        return ""

    def rebuild(self) -> None:
        self._timer.stop()
        self.matrix = None
        problem = self.problem()
        if problem:
            self.status.setText(problem)
            self.grid.clear()
            self.grid.setRowCount(0)
            self.grid.setColumnCount(0)
            return
        table = self.table
        oids = shown_oids(table)
        slots = table.slots_of(oids)

        def labels(name: str) -> list[str]:
            return table.column(name).displays(slots)

        pages_name = self.pages_axis.currentText()
        pages = labels(pages_name) if pages_name != NONE else None
        try:
            self.matrix = build_matrix(
                labels(self.rows_axis.currentText()),
                labels(self.columns_axis.currentText()),
                numbers(table, self.activity.currentText(), oids),
                pages,
                predict=self.predict.isChecked(),
            )
        except ValueError as err:
            self.status.setText(str(err))
            self.grid.setRowCount(0)
            self.grid.setColumnCount(0)
            return
        self.oids = oids
        page = self.page.currentText()
        self.page.blockSignals(True)
        self.page.clear()
        self.page.addItems(self.matrix.page_labels if pages is not None else [])
        self.page.setEnabled(pages is not None)
        if page:
            self.page.setCurrentText(page)
        self.page.blockSignals(False)
        self.show_page()

    def show_page(self) -> None:
        m = self.matrix
        if m is None:
            return
        page = self.page.currentText() if self.page.isEnabled() else ""
        aggregate = self.aggregate.currentData()
        values = m.values(aggregate)
        finite = values[np.isfinite(values)]
        low, high = (float(finite.min()), float(finite.max())) if len(finite) else (0.0, 0.0)
        better = self.higher_is_better.isChecked()
        grid = self.grid
        grid.clear()
        grid.setRowCount(len(m.row_labels))
        grid.setColumnCount(len(m.column_labels))
        grid.setVerticalHeaderLabels(m.row_labels)
        grid.setHorizontalHeaderLabels(m.column_labels)
        italic = QFont()
        italic.setItalic(True)
        made = estimated = 0
        for r, row in enumerate(m.row_labels):
            for c, column in enumerate(m.column_labels):
                item = QTableWidgetItem()
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                cell = m.cell(page, row, column)
                if cell is not None:
                    value = cell.value(aggregate)
                    count = len(cell.rows)
                    text = number_text(value, 2) if math.isfinite(value) else "–"
                    item.setText(text + (f" ({count})" if count > 1 else ""))
                    item.setBackground(QBrush(heat(value, low, high, better)))
                    item.setData(OIDS_ROLE, self.oids[cell.rows].tolist())
                    shown = ", ".join(f"{v:.2f}" for v in cell.values[:12])
                    item.setToolTip(f"{row} × {column}: {count} compound(s)\n{shown}")
                    made += 1
                else:
                    guess = m.predicted(page, row, column)
                    if math.isfinite(guess):
                        item.setText(f"{guess:.2f}")
                        item.setFont(italic)
                        item.setForeground(QBrush(QColor(130, 130, 130)))
                        item.setToolTip(f"{row} × {column}: not made; Free-Wilson estimate")
                        estimated += 1
                grid.setItem(r, c, item)
        cells = len(m.row_labels) * len(m.column_labels)
        text = f"{made:,} of {cells:,} cell(s) made"
        if estimated:
            text += f"; {estimated:,} grey estimate(s)"
        if m.model is not None:
            text += f"; additive model R² {m.model.fit.r2:.2f}"
        if m.dropped:
            text += f"; {m.dropped:,} row(s) without labels left out"
        self.status.setText(text + ". Click a cell to select its compounds.")

    def select_cell(self, row: int, column: int) -> list[int]:
        item = self.grid.item(row, column)
        oids = item.data(OIDS_ROLE) if item is not None else None
        if oids:
            self.ctx.selection.set(oids, "sar_matrix")
        return oids or []

    def copy_text(self) -> str:
        from PySide6.QtGui import QGuiApplication

        grid = self.grid
        heads = [grid.horizontalHeaderItem(c).text() for c in range(grid.columnCount())]
        lines = ["\t".join([""] + heads)]
        for r in range(grid.rowCount()):
            cells = [grid.item(r, c).text() if grid.item(r, c) else "" for c in range(len(heads))]
            lines.append("\t".join([grid.verticalHeaderItem(r).text(), *cells]))
        text = "\n".join(lines)
        QGuiApplication.clipboard().setText(text)
        return text

    def state(self) -> dict:
        return {
            "rows": self.rows_axis.currentText(),
            "columns": self.columns_axis.currentText(),
            "pages": self.pages_axis.currentText(),
            "activity": self.activity.currentText(),
            "aggregate": self.aggregate.currentData(),
            "predict": self.predict.isChecked(),
            "higher_is_better": self.higher_is_better.isChecked(),
        }

    def set_state(self, state: dict) -> None:
        for combo, key in (
            (self.rows_axis, "rows"),
            (self.columns_axis, "columns"),
            (self.pages_axis, "pages"),
            (self.activity, "activity"),
        ):
            index = combo.findText(str(state.get(key, "")))
            if index >= 0:
                combo.setCurrentIndex(index)
        self.aggregate.setCurrentIndex(max(0, self.aggregate.findData(state.get("aggregate"))))
        self.predict.setChecked(bool(state.get("predict", True)))
        self.higher_is_better.setChecked(bool(state.get("higher_is_better", True)))
        self.rebuild()

    def close_panel(self) -> None:
        self._timer.stop()
        self._unsubscribe()
