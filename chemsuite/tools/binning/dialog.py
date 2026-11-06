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
"""The Bin Column dialog, with a live preview of the bins and their counts."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from chemsuite.analysis.binning import (
    METHODS,
    Binned,
    BinSpec,
    bin_values,
    parse_names,
    parse_numbers,
)
from chemsuite.app.dialogs import selected_only
from chemsuite.app.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers


def _spin(value: float = 0.0) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(-1e12, 1e12)
    spin.setDecimals(4)
    spin.setValue(value)
    return spin


class BinningDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Bin Column")
        self.ctx = ctx
        self.result: Binned | None = None

        self.column = QComboBox()
        self.column.addItems(number_columns(ctx.table))
        current = ctx.current_column()
        if current and self.column.findText(current) >= 0:
            self.column.setCurrentText(current)
        self.method = QComboBox()
        for key, title in METHODS.items():
            self.method.addItem(title, key)
        self.count = QSpinBox()
        self.count.setRange(2, 256)
        self.count.setValue(3)
        self.use_width = QCheckBox("Fixed width:")
        self.width = _spin(1.0)
        self.width.setMinimum(1e-9)
        self.use_min = QCheckBox("From:")
        self.minimum = _spin()
        self.use_max = QCheckBox("To:")
        self.maximum = _spin()
        self.edges = QLineEdit()
        self.edges.setPlaceholderText("e.g. 5, 6, 7")
        self.open_ends = QCheckBox("Add bins below the first and above the last edge")
        self.open_ends.setChecked(True)
        self.names = QLineEdit()
        self.names.setPlaceholderText("Low, Medium, High…  (one per bin, low to high)")
        self.range_name = QLineEdit()
        self.class_name = QLineEdit()
        self.write_range = QCheckBox("Range column:")
        self.write_range.setChecked(True)
        self.write_class = QCheckBox("Class column:")
        self.write_class.setChecked(True)
        self.selected_only = selected_only(ctx)
        self.preview = QTableWidget(0, 3)
        self.preview.setHorizontalHeaderLabels(["Range", "Class", "Rows"])
        self.preview.verticalHeader().hide()
        self.preview.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.preview.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Add Columns")

        width_row = QHBoxLayout()
        width_row.addWidget(self.use_width)
        width_row.addWidget(self.width, 1)
        limits = QHBoxLayout()
        for box, spin in ((self.use_min, self.minimum), (self.use_max, self.maximum)):
            limits.addWidget(box)
            limits.addWidget(spin, 1)
        self.form = QFormLayout()
        self.form.addRow("Column:", self.column)
        self.form.addRow("Method:", self.method)
        self.form.addRow("Bins:", self.count)
        self.form.addRow("", width_row)
        self.form.addRow("Limits:", limits)
        self.form.addRow("Edges:", self.edges)
        self.form.addRow("", self.open_ends)
        self.form.addRow("Class names:", self.names)
        self.form.addRow(self.write_range, self.range_name)
        self.form.addRow(self.write_class, self.class_name)
        layout = QVBoxLayout(self)
        layout.addLayout(self.form)
        layout.addWidget(self.selected_only)
        layout.addWidget(self.preview, 1)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(520, 600)

        self.column.currentIndexChanged.connect(self._column_changed)
        for signal in (
            self.method.currentIndexChanged,
            self.count.valueChanged,
            self.use_width.toggled,
            self.width.valueChanged,
            self.use_min.toggled,
            self.minimum.valueChanged,
            self.use_max.toggled,
            self.maximum.valueChanged,
            self.edges.textChanged,
            self.open_ends.toggled,
            self.names.textChanged,
            self.selected_only.toggled,
            self.write_range.toggled,
            self.write_class.toggled,
        ):
            signal.connect(lambda *_: self.update_preview())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._column_changed()

    def _column_changed(self, *_args) -> None:
        name = self.column.currentText()
        self.range_name.setText(f"{name} Range" if name else "")
        self.class_name.setText(f"{name} Class" if name else "")
        values = self._values()
        finite = values[np.isfinite(values)]
        if len(finite):
            for spin, value in ((self.minimum, finite.min()), (self.maximum, finite.max())):
                spin.blockSignals(True)
                spin.setValue(float(value))
                spin.blockSignals(False)
        self.update_preview()

    def _values(self) -> np.ndarray:
        name = self.column.currentText()
        if not name:
            return np.empty(0)
        oids = self.ctx.target_oids(self.selected_only.isChecked())
        return numbers(self.ctx.table, name, oids)

    def spec(self) -> BinSpec:
        method = self.method.currentData()
        names = parse_names(self.names.text())
        if method == "custom":
            return BinSpec(
                "custom",
                edges=parse_numbers(self.edges.text()),
                open_ends=self.open_ends.isChecked(),
                names=names,
            )
        if method == "equal_count":
            return BinSpec("equal_count", self.count.value(), names=names)
        return BinSpec(
            "equal_width",
            self.count.value(),
            width=self.width.value() if self.use_width.isChecked() else None,
            minimum=self.minimum.value() if self.use_min.isChecked() else None,
            maximum=self.maximum.value() if self.use_max.isChecked() else None,
            names=names,
        )

    def update_preview(self) -> None:
        method = self.method.currentData()
        width = method == "equal_width"
        self.form.setRowVisible(2, method != "custom")
        self.form.setRowVisible(3, width)
        self.form.setRowVisible(4, width)
        self.form.setRowVisible(5, method == "custom")
        self.form.setRowVisible(6, method == "custom")
        self.count.setEnabled(not (width and self.use_width.isChecked()))
        self.width.setEnabled(self.use_width.isChecked())
        self.minimum.setEnabled(self.use_min.isChecked())
        self.maximum.setEnabled(self.use_max.isChecked())
        self.range_name.setEnabled(self.write_range.isChecked())
        self.class_name.setEnabled(self.write_class.isChecked())
        self.result = None
        self.preview.setRowCount(0)
        ok = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        try:
            result = bin_values(self._values(), self.spec())
        except ValueError as err:
            self.status.setText(str(err))
            ok.setEnabled(False)
            return
        self.result = result
        counts = result.counts()
        self.preview.setRowCount(len(result.intervals))
        for row, (rng, cls, count) in enumerate(
            zip(result.range_names, result.class_names, counts.tolist(), strict=True)
        ):
            for col, text in enumerate((rng, cls, f"{count:,}")):
                self.preview.setItem(row, col, QTableWidgetItem(text))
        unbinned = int((result.index < 0).sum())
        self.status.setText(
            f"{int(counts.sum()):,} row(s) in {len(counts)} bin(s)"
            + (f"; {unbinned:,} blank or outside every bin" if unbinned else "")
        )
        ok.setEnabled(self.write_range.isChecked() or self.write_class.isChecked())

    def accept(self) -> None:
        if self.result is None:
            return
        self.apply()
        super().accept()

    def apply(self) -> dict[str, str]:
        result = self.result
        assert result is not None
        oids = self.ctx.target_oids(self.selected_only.isChecked())
        source = self.column.currentText()
        results = {}
        if self.write_range.isChecked():
            name = self.range_name.text().strip() or f"{source} Range"
            results[name] = (oids, result.labels(result.range_names))
        if self.write_class.isChecked():
            name = self.class_name.text().strip() or f"{source} Class"
            results[name] = (oids, result.labels(result.class_names))
        names = self.ctx.apply_results(results, label=f"Bin {source}", overwrite=True)
        self.ctx.show_message(f"Binned {source} into {len(result.intervals)} bin(s)")
        return names
