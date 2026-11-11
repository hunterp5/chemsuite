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
"""The Statistics panel: one tab per analysis, all over the same rows (shown, all or selected).

Tabs recompute when the table, filters or selection change, a moment after the last change,
and only while visible. Correlations on big tables run as a background job on request.
"""

from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QBrush, QColor, QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.statistics import (
    CORRELATION_METHODS,
    CURVE_MODELS,
    OUTLIER_METHODS,
    SUMMARY_FIELDS,
    TESTS,
    Correlation,
    correlation_matrix,
    find_outliers,
    fit_curve,
    run_test,
    summarize,
)
from chemsuite.app.flow_layout import FlowRow, compact, field
from chemsuite.app.jobs import QUICK
from chemsuite.app.tools import ToolContext, ToolPanel
from chemsuite.core.table import Changed
from chemsuite.core.table.columns import TEXT
from chemsuite.plotting.data import number_columns, numbers, shown_oids

SCOPES = {"shown": "Shown rows", "all": "All rows", "selected": "Selected rows"}
DELAY_MS = 300
#: Correlations above this many cells (rows × columns) wait for the Calculate button.
AUTO_CORRELATION_CELLS = 2_000_000
MAX_OUTLIERS_LISTED = 1000
ROW_LABEL = "(row number)"


def _number(value: float, digits: int = 4) -> str:
    if isinstance(value, int | np.integer):
        return f"{value:,}"
    return "" if not math.isfinite(value) else f"{value:.{digits}g}"


def _read_only(table: QTableWidget) -> None:
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
    table.verticalHeader().setDefaultSectionSize(22)


def copy_table(widget: QTableWidget) -> str:
    """The table as tab-separated text with headers (also put on the clipboard)."""
    columns = range(widget.columnCount())
    head = [""] + [widget.horizontalHeaderItem(c).text() for c in columns]
    lines = ["\t".join(head)]
    for r in range(widget.rowCount()):
        label = widget.verticalHeaderItem(r)
        cells = [(widget.item(r, c).text() if widget.item(r, c) else "") for c in columns]
        lines.append("\t".join([label.text() if label else "", *cells]))
    text = "\n".join(lines)
    QGuiApplication.clipboard().setText(text)
    return text


def correlation_color(r: float) -> QColor:
    """Blue for negative, white for none, red for positive."""
    if not math.isfinite(r):
        return QColor(235, 235, 235)
    t = min(1.0, abs(r))
    if r >= 0:
        return QColor(255, int(255 - 155 * t), int(255 - 175 * t))
    return QColor(int(255 - 175 * t), int(255 - 125 * t), 255)


class StatisticsPanel(ToolPanel):
    tool_id = "statistics"

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__()
        self.ctx = ctx
        self.table = ctx.table
        self.scope = QComboBox()
        for key, title in SCOPES.items():
            self.scope.addItem(title, key)
        self.rows_label = QLabel()
        self.tabs = QTabWidget()
        self.summary = SummaryTab(self)
        self.correlation = CorrelationTab(self)
        self.outliers = OutlierTab(self)
        self.fit = FitTab(self)
        self.tests = TestTab(self)
        for tab, title in (
            (self.summary, "Summary"),
            (self.correlation, "Correlation"),
            (self.outliers, "Outliers"),
            (self.fit, "Curve Fit"),
            (self.tests, "Tests"),
        ):
            self.tabs.addTab(tab, title)
        top = QHBoxLayout()
        top.addWidget(QLabel("Rows:"))
        top.addWidget(self.scope)
        top.addWidget(self.rows_label, 1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(top)
        layout.addWidget(self.tabs, 1)

        self._stale: set[QWidget] = set()
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(DELAY_MS)
        self._timer.timeout.connect(self.refresh)
        self.scope.currentIndexChanged.connect(lambda _i: self.invalidate())
        self.tabs.currentChanged.connect(lambda _i: self.refresh())
        self._unsubscribe_table = self.table.subscribe(self._on_table_event)
        self._unsubscribe_selection = ctx.selection.subscribe(self._on_selection)
        self.columns_changed()
        self.invalidate(now=True)

    # -- rows and values ----------------------------------------------------------------------
    def oids(self) -> np.ndarray:
        scope = self.scope.currentData()
        if scope == "selected":
            return self.ctx.selection.as_array()
        if scope == "all":
            return self.table.oids.copy()
        return shown_oids(self.table)

    def values(self, name: str, oids: np.ndarray | None = None) -> np.ndarray:
        return numbers(self.table, name, self.oids() if oids is None else oids)

    def number_columns(self) -> list[str]:
        return number_columns(self.table)

    def title(self) -> str:
        return "Statistics"

    # -- refreshing ---------------------------------------------------------------------------
    def _on_table_event(self, event: object) -> None:
        if isinstance(event, Changed):
            if event.structural:
                self.columns_changed()
            self.invalidate()

    def _on_selection(self, _oids: frozenset[int], _origin: str) -> None:
        if self.scope.currentData() == "selected":
            self.invalidate()

    def columns_changed(self) -> None:
        for tab in (self.outliers, self.fit, self.tests):
            tab.fill_columns()

    def invalidate(self, *, now: bool = False) -> None:
        """Mark every tab out of date; the visible one recomputes shortly."""
        self._stale = {self.summary, self.correlation, self.outliers, self.fit, self.tests}
        if now:
            self.refresh()
        else:
            self._timer.start()

    def refresh(self) -> None:
        self._timer.stop()
        self.rows_label.setText(f"{len(self.oids()):,} row(s)")
        tab = self.tabs.currentWidget()
        # A panel in a hidden dock tab catches up when shown (window off screen: scripts, tests).
        on_screen = self.isVisible() or not self.window().isVisible()
        if tab in self._stale and on_screen:
            self._stale.discard(tab)
            tab.refresh()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self.refresh()

    # -- persistence --------------------------------------------------------------------------
    def state(self) -> dict:
        return {
            "scope": self.scope.currentData(),
            "tab": self.tabs.currentIndex(),
            "correlation": self.correlation.method.currentData(),
        }

    def set_state(self, state: dict) -> None:
        self.scope.setCurrentIndex(max(0, self.scope.findData(state.get("scope"))))
        method = self.correlation.method
        method.setCurrentIndex(max(0, method.findData(state.get("correlation"))))
        tab = state.get("tab")
        if isinstance(tab, int) and 0 <= tab < self.tabs.count():
            self.tabs.setCurrentIndex(tab)
        self.invalidate(now=True)

    def close_panel(self) -> None:
        self._timer.stop()
        self._unsubscribe_table()
        self._unsubscribe_selection()
        self.correlation.cancel()


class SummaryTab(QWidget):
    """One row per number column: N, missing, mean, SD, percentiles, skew, kurtosis."""

    def __init__(self, panel: StatisticsPanel) -> None:
        super().__init__()
        self.panel = panel
        self.grid = QTableWidget(0, len(SUMMARY_FIELDS))
        self.grid.setHorizontalHeaderLabels([f[0] for f in SUMMARY_FIELDS])
        for i, (_h, _a, help_text) in enumerate(SUMMARY_FIELDS):
            if help_text:
                self.grid.horizontalHeaderItem(i).setToolTip(help_text)
        _read_only(self.grid)
        self.grid.setToolTip("Double-click a column to see its histogram.")
        self.grid.cellDoubleClicked.connect(self._histogram)
        self.note = QLabel()
        copy = QPushButton("Copy Table")
        copy.clicked.connect(lambda: copy_table(self.grid))
        bottom = QHBoxLayout()
        bottom.addWidget(self.note, 1)
        bottom.addWidget(copy)
        layout = QVBoxLayout(self)
        layout.addWidget(self.grid, 1)
        layout.addLayout(bottom)

    def refresh(self) -> None:
        names = self.panel.number_columns()
        oids = self.panel.oids()
        self.grid.setRowCount(len(names))
        self.grid.setVerticalHeaderLabels(names)
        for row, name in enumerate(names):
            summary = summarize(name, self.panel.values(name, oids))
            for col, (_h, attr, _help) in enumerate(SUMMARY_FIELDS):
                item = QTableWidgetItem(_number(getattr(summary, attr)))
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.grid.setItem(row, col, item)
        self.grid.resizeColumnsToContents()
        self.note.setText(
            ""
            if names
            else "No number columns. Calculate descriptors or convert a column to numbers."
        )

    def _histogram(self, row: int, _column: int) -> None:
        name = self.grid.verticalHeaderItem(row).text()
        plot = self.panel.ctx.open_tool("histogram")
        plot.set_spec(x=name)


class CorrelationTab(QWidget):
    """A coloured matrix of pairwise correlations; double-click a cell for its scatter plot."""

    def __init__(self, panel: StatisticsPanel) -> None:
        super().__init__()
        self.panel = panel
        self.result: Correlation | None = None
        self.job = None
        self.method = QComboBox()
        for key, title in CORRELATION_METHODS.items():
            self.method.addItem(title, key)
        self.calculate_button = QPushButton("Calculate")
        self.note = QLabel()
        self.grid = QTableWidget(0, 0)
        _read_only(self.grid)
        self.grid.horizontalHeader().setDefaultSectionSize(64)
        self.grid.cellDoubleClicked.connect(self._scatter)
        copy = QPushButton("Copy Table")
        copy.clicked.connect(lambda: copy_table(self.grid))
        top = FlowRow()
        top.addWidget(field("Method:", self.method))
        top.addWidget(self.calculate_button)
        top.addWidget(copy)
        self.note.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.addWidget(top)
        layout.addWidget(self.note)
        layout.addWidget(self.grid, 1)
        self.method.currentIndexChanged.connect(lambda _i: self.refresh())
        self.calculate_button.clicked.connect(lambda: self.calculate())

    def refresh(self) -> None:
        names = self.panel.number_columns()
        cells = len(self.panel.oids()) * len(names)
        if cells <= AUTO_CORRELATION_CELLS:
            self.calculate(background=False)
        else:
            self.note.setText("The data changed: press Calculate.")

    def calculate(self, *, background: bool = True):
        self.cancel()
        panel = self.panel
        names = panel.number_columns()
        oids = panel.oids()
        columns = [panel.values(n, oids) for n in names]
        method = self.method.currentData()
        if not background:
            self.display(correlation_matrix(names, columns, method))
            return None

        def work(job):
            return correlation_matrix(
                names, columns, method, progress=job.report, cancel=job.cancel_token
            )

        self.job = panel.ctx.run("Calculating correlations…", work, self.display, lane=QUICK)
        return self.job

    def cancel(self) -> None:
        if self.job is not None and self.job.running:
            self.job.cancel()
        self.job = None

    def display(self, result: Correlation) -> None:
        self.result = result
        names = list(result.names)
        grid = self.grid
        grid.setRowCount(len(names))
        grid.setColumnCount(len(names))
        grid.setHorizontalHeaderLabels(names)
        grid.setVerticalHeaderLabels(names)
        for i in range(len(names)):
            for j in range(len(names)):
                r = float(result.r[i, j])
                item = QTableWidgetItem("" if not math.isfinite(r) else f"{r:.2f}")
                item.setBackground(QBrush(correlation_color(r)))
                item.setForeground(QBrush(QColor(20, 20, 20)))
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                item.setToolTip(f"{names[i]} vs {names[j]}\nr = {r:.4f}\nn = {result.n[i, j]:,}")
                grid.setItem(i, j, item)
        grid.resizeColumnsToContents()
        self.note.setText(
            "Double-click a cell to plot the pair." if names else "No number columns."
        )

    def _scatter(self, row: int, column: int) -> None:
        if self.result is None or row == column:
            return
        plot = self.panel.ctx.open_tool("scatter")
        plot.set_spec(x=self.result.names[column], y=self.result.names[row])


def _form() -> QFormLayout:
    """A form whose fields drop below their labels when the panel is docked narrow (Windows'
    wider controls would otherwise set its minimum width)."""
    form = QFormLayout()
    form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
    return form


def _column_combo() -> QComboBox:
    combo = QComboBox()
    compact(combo, chars=10)
    return combo


def _refill(combo: QComboBox, names: list[str], *, blank: str | None = None) -> None:
    """Replace the items, keeping the current choice when it still exists."""
    current = combo.currentText()
    combo.blockSignals(True)
    combo.clear()
    if blank is not None:
        combo.addItem(blank)
    combo.addItems(names)
    if current and combo.findText(current) >= 0:
        combo.setCurrentText(current)
    combo.blockSignals(False)


class OutlierTab(QWidget):
    def __init__(self, panel: StatisticsPanel) -> None:
        super().__init__()
        self.panel = panel
        self.flagged = np.empty(0, dtype=np.int64)
        self.column = _column_combo()
        self.method = QComboBox()
        for key, (title, _label, _default) in OUTLIER_METHODS.items():
            self.method.addItem(title, key)
        compact(self.method, chars=12)
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.1, 100.0)
        self.threshold.setSingleStep(0.5)
        self.label_column = _column_combo()
        self.description = QLabel()
        self.description.setWordWrap(True)
        self.grid = QTableWidget(0, 2)
        self.grid.setHorizontalHeaderLabels(["Row", "Value"])
        _read_only(self.grid)
        self.grid.horizontalHeader().setStretchLastSection(True)
        self.select_button = QPushButton("Select These Rows")
        self.select_button.clicked.connect(self.select_rows)
        form = _form()
        form.addRow("Column:", self.column)
        form.addRow("Method:", self.method)
        self.threshold_label = QLabel()
        form.addRow(self.threshold_label, self.threshold)
        form.addRow("Label rows by:", self.label_column)
        bottom = QHBoxLayout()
        bottom.addStretch(1)
        bottom.addWidget(self.select_button)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.description)
        layout.addWidget(self.grid, 1)
        layout.addLayout(bottom)
        self.method.currentIndexChanged.connect(self._method_changed)
        for signal in (
            self.column.currentIndexChanged,
            self.threshold.valueChanged,
            self.label_column.currentIndexChanged,
        ):
            signal.connect(lambda *_: self.refresh())
        self._method_changed()

    def fill_columns(self) -> None:
        _refill(self.column, self.panel.number_columns())
        table = self.panel.table
        texts = [n for n in table.column_names if table.kind(n) == TEXT]
        _refill(self.label_column, texts, blank=ROW_LABEL)
        if self.label_column.currentText() == ROW_LABEL:
            preferred = next((n for n in ("ID", "Name", "Title") if n in texts), None)
            if preferred:
                self.label_column.setCurrentText(preferred)

    def _method_changed(self, *_args) -> None:
        _title, label, default = OUTLIER_METHODS[self.method.currentData()]
        self.threshold_label.setText(f"{label[0].upper()}{label[1:]}:")
        self.threshold.blockSignals(True)
        self.threshold.setValue(default)
        self.threshold.blockSignals(False)
        self.refresh()

    def refresh(self) -> None:
        name = self.column.currentText()
        self.grid.setRowCount(0)
        self.flagged = np.empty(0, dtype=np.int64)
        if not name:
            self.description.setText("No number columns.")
            self.select_button.setEnabled(False)
            return
        panel = self.panel
        oids = panel.oids()
        values = panel.values(name, oids)
        result = find_outliers(values, self.method.currentData(), self.threshold.value())
        hits = np.flatnonzero(result.mask)
        hits = hits[np.argsort(-np.abs(values[hits] - np.nanmedian(values)))]  # most extreme first
        self.flagged = oids[hits]
        shown = hits[:MAX_OUTLIERS_LISTED]
        label = self.label_column.currentText()
        table = panel.table
        if label and label != ROW_LABEL and label in table:
            labels = table.column(label).displays(table.slots_of(oids[shown]))
        else:  # the row number the table shows (hidden rows have none)
            position = np.full(len(table), -1, dtype=np.int64)
            view = table.view()
            position[view] = np.arange(len(view))
            rows = position[table.slots_of(oids[shown])]
            labels = [f"Row {r + 1:,}" if r >= 0 else "(hidden row)" for r in rows.tolist()]
        self.grid.setRowCount(len(shown))
        for i, (text, value) in enumerate(zip(labels, values[shown], strict=True)):
            self.grid.setItem(i, 0, QTableWidgetItem(text))
            self.grid.setItem(i, 1, QTableWidgetItem(_number(float(value), 6)))
        more = f" (listing {MAX_OUTLIERS_LISTED:,})" if len(hits) > MAX_OUTLIERS_LISTED else ""
        self.description.setText(f"{result.description}\n{len(hits):,} outlier(s){more}")
        self.select_button.setEnabled(len(hits) > 0)

    def select_rows(self) -> None:
        self.panel.ctx.selection.set(self.flagged.tolist(), "statistics")


class FitTab(QWidget):
    def __init__(self, panel: StatisticsPanel) -> None:
        super().__init__()
        self.panel = panel
        self.x = _column_combo()
        self.y = _column_combo()
        self.model = QComboBox()
        for key, title in CURVE_MODELS.items():
            self.model.addItem(title, key)
        compact(self.model, chars=12)
        self.degree = QSpinBox()
        self.degree.setRange(1, 6)
        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        form = _form()
        form.addRow("X:", self.x)
        form.addRow("Y:", self.y)
        form.addRow("Model:", self.model)
        form.addRow("Degree:", self.degree)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.report, 1)
        for signal in (
            self.x.currentIndexChanged,
            self.y.currentIndexChanged,
            self.model.currentIndexChanged,
            self.degree.valueChanged,
        ):
            signal.connect(lambda *_: self.refresh())

    def fill_columns(self) -> None:
        names = self.panel.number_columns()
        _refill(self.x, names)
        _refill(self.y, names)
        if self.y.currentIndex() == 0 and len(names) > 1 and self.x.currentIndex() == 0:
            self.y.setCurrentIndex(1)

    def refresh(self) -> None:
        self.degree.setEnabled(self.model.currentData() == "polynomial")
        x, y = self.x.currentText(), self.y.currentText()
        if not x or not y:
            self.report.setPlainText("Needs number columns.")
            return
        oids = self.panel.oids()
        result = fit_curve(
            self.panel.values(x, oids),
            self.panel.values(y, oids),
            self.model.currentData(),
            degree=self.degree.value(),
            x_name=x,
            y_name=y,
        )
        self.report.setPlainText(result.report())


class TestTab(QWidget):
    def __init__(self, panel: StatisticsPanel) -> None:
        super().__init__()
        self.panel = panel
        self.test = QComboBox()
        for key, (title, _columns) in TESTS.items():
            self.test.addItem(title, key)
        compact(self.test, chars=12)
        self.a = _column_combo()
        self.b = _column_combo()
        self.mean = QDoubleSpinBox()
        self.mean.setRange(-1e9, 1e9)
        self.mean.setDecimals(4)
        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        self.form = _form()
        self.form.addRow("Test:", self.test)
        self.form.addRow("Column:", self.a)
        self.form.addRow("Second column:", self.b)
        self.form.addRow("Compare with mean:", self.mean)
        layout = QVBoxLayout(self)
        layout.addLayout(self.form)
        layout.addWidget(self.report, 1)
        for signal in (
            self.test.currentIndexChanged,
            self.a.currentIndexChanged,
            self.b.currentIndexChanged,
            self.mean.valueChanged,
        ):
            signal.connect(lambda *_: self.refresh())

    def fill_columns(self) -> None:
        names = self.panel.number_columns()
        _refill(self.a, names)
        _refill(self.b, names)
        if self.a.currentIndex() == self.b.currentIndex() and len(names) > 1:
            self.b.setCurrentIndex(1)

    def refresh(self) -> None:
        test = self.test.currentData()
        two = TESTS[test][1] == 2
        self.form.setRowVisible(2, two)
        self.form.setRowVisible(3, test == "one_sample_t")
        a, b = self.a.currentText(), self.b.currentText()
        if not a or (two and not b):
            self.report.setPlainText("Needs number columns.")
            return
        oids = self.panel.oids()
        result = run_test(
            test,
            self.panel.values(a, oids),
            self.panel.values(b, oids) if two else None,
            a_name=a,
            b_name=b,
            mean=self.mean.value(),
        )
        self.report.setPlainText(result.report())
