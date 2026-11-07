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
"""The MPO dialog: one row per property, each with a desirability shape and weight."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.mpo import COMBINE, SHAPES, Desirability, MpoProfile
from chemsuite.app.dialogs import selected_only
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import number_columns, numbers

HEADERS = ("Column", "Shape", "Low", "High", "Target", "Width", "Weight", "Step")
COLUMN, SHAPE, LOW, HIGH, TARGET, WIDTH, WEIGHT, HARD = range(len(HEADERS))
DEFAULT_OUTPUT = "MPO Score"
FILE_FILTER = "MPO profiles (*.json)"

_last: MpoProfile | None = None  # remembered for the session
_last_output = DEFAULT_OUTPUT


def _spin(value: float = 0.0, *, minimum: float = -1e9) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(minimum, 1e9)
    spin.setDecimals(3)
    spin.setValue(value)
    spin.setButtonSymbols(QDoubleSpinBox.ButtonSymbols.NoButtons)
    return spin


class MpoDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("MPO Score")
        self.ctx = ctx
        self.columns = number_columns(ctx.table)

        self.grid = QTableWidget(0, len(HEADERS))
        self.grid.setHorizontalHeaderLabels(HEADERS)
        self.grid.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.grid.verticalHeader().hide()
        header = self.grid.horizontalHeader()
        # Cell widgets do not size their columns: give each a fixed width.
        widths = {SHAPE: 230, LOW: 84, HIGH: 84, TARGET: 84, WIDTH: 84, WEIGHT: 70, HARD: 44}
        for col, width in widths.items():
            header.setSectionResizeMode(col, QHeaderView.ResizeMode.Fixed)
            self.grid.setColumnWidth(col, width)
        header.setSectionResizeMode(COLUMN, QHeaderView.ResizeMode.Stretch)
        tips = {
            LOW: "Maximise: 0 here. Minimise: 1 here. Range: the lower limit.",
            HIGH: "Maximise: 1 here. Minimise: 0 here. Range: the upper limit.",
            TARGET: "Range: the best value (blank: the middle). Bell curve: the centre.",
            WIDTH: "Bell curve: the spread (one standard deviation).",
            WEIGHT: "How much this property counts (0 ignores it).",
            HARD: "A step instead of a slope: 1 on the good side of the limit, 0 elsewhere.",
        }
        for column, tip in tips.items():
            self.grid.horizontalHeaderItem(column).setToolTip(tip)
        add = QPushButton("Add Property")
        remove = QPushButton("Remove")
        load = QPushButton("Load…")
        save = QPushButton("Save…")
        row_buttons = QHBoxLayout()
        for button in (add, remove):
            row_buttons.addWidget(button)
        row_buttons.addStretch(1)
        for button in (load, save):
            row_buttons.addWidget(button)

        self.method = QComboBox()
        for key, title in COMBINE.items():
            self.method.addItem(title, key)
        self.method.setToolTip(
            "Geometric: one property scoring 0 makes the whole score 0.\n"
            "Arithmetic: a good property can make up for a poor one."
        )
        self.skip_missing = QCheckBox("Score rows missing a property on the others")
        self.output = QLineEdit(_last_output)
        self.per_property = QCheckBox("Also add a desirability column per property")
        self.selected_only = selected_only(ctx)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Score")

        form = QFormLayout()
        form.addRow("Combine with:", self.method)
        form.addRow("", self.skip_missing)
        form.addRow("Result column:", self.output)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Properties:"))
        layout.addWidget(self.grid, 1)
        layout.addLayout(row_buttons)
        layout.addLayout(form)
        layout.addWidget(self.per_property)
        layout.addWidget(self.selected_only)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(960, 520)

        add.clicked.connect(lambda: self.add_property())
        remove.clicked.connect(self.remove_selected)
        load.clicked.connect(self.load_dialog)
        save.clicked.connect(self.save_dialog)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        if _last is not None:
            self.set_profile(_last)
        self._update()

    # -- rows ---------------------------------------------------------------------------------
    def add_property(self, spec: Desirability | None = None) -> int:
        row = self.grid.rowCount()
        self.grid.insertRow(row)
        column = QComboBox()
        column.addItems(self.columns)
        shape = QComboBox()
        for key, title in SHAPES.items():
            shape.addItem(title, key)
        target = QLineEdit()
        target.setPlaceholderText("middle")
        hard = QCheckBox()
        holder = QWidget()  # centre the check box in its cell
        box = QHBoxLayout(holder)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(hard, alignment=Qt.AlignmentFlag.AlignCenter)
        widgets = {
            COLUMN: column,
            SHAPE: shape,
            LOW: _spin(),
            HIGH: _spin(1.0),
            TARGET: target,
            WIDTH: _spin(1.0, minimum=0.0),
            WEIGHT: _spin(1.0, minimum=0.0),
            HARD: holder,
        }
        for col, widget in widgets.items():
            self.grid.setCellWidget(row, col, widget)
        hard.toggled.connect(lambda _on: self._update())
        for col in (LOW, HIGH, WIDTH, WEIGHT):
            widgets[col].valueChanged.connect(lambda _v: self._update())
        target.textChanged.connect(lambda _t: self._update())
        shape.currentIndexChanged.connect(lambda _i: self._update())
        if spec is None:
            used = {self._widget(r, COLUMN).currentText() for r in range(row)}
            name = next((c for c in self.columns if c not in used), "")
            if name:
                column.setCurrentText(name)
                self._default_range(row)
            column.currentIndexChanged.connect(
                lambda _i, w=column: self._default_range(self._row_of(w))
            )
        else:
            self._set_row(row, spec)
            column.currentIndexChanged.connect(
                lambda _i, w=column: self._default_range(self._row_of(w))
            )
        self._update()
        return row

    def _row_of(self, widget: QWidget) -> int:
        for row in range(self.grid.rowCount()):
            if self.grid.cellWidget(row, COLUMN) is widget:
                return row
        return -1

    def _widget(self, row: int, col: int):
        widget = self.grid.cellWidget(row, col)
        return widget.findChild(QCheckBox) if col == HARD else widget

    def _default_range(self, row: int) -> None:
        """Start from the column's 10th to 90th percentile, a sensible first guess."""
        if row < 0:
            return
        name = self._widget(row, COLUMN).currentText()
        values = numbers(self.ctx.table, name, self.ctx.table.oids)
        values = values[np.isfinite(values)]
        if not len(values):
            return
        low, high = np.percentile(values, [10, 90])
        if high == low:
            high = low + 1
        self._widget(row, LOW).setValue(float(low))
        self._widget(row, HIGH).setValue(float(high))
        self._widget(row, WIDTH).setValue(float(max(np.std(values), 1e-3)))

    def _set_row(self, row: int, spec: Desirability) -> None:
        column = self._widget(row, COLUMN)
        if column.findText(spec.column) < 0:
            column.addItem(spec.column)  # a profile from another table: shown, flagged below
        column.setCurrentText(spec.column)
        shape = self._widget(row, SHAPE)
        shape.setCurrentIndex(max(0, shape.findData(spec.shape)))
        self._widget(row, LOW).setValue(spec.low)
        self._widget(row, HIGH).setValue(spec.high)
        self._widget(row, TARGET).setText("" if spec.target is None else f"{spec.target:g}")
        self._widget(row, WIDTH).setValue(spec.width)
        self._widget(row, WEIGHT).setValue(spec.weight)
        self._widget(row, HARD).setChecked(spec.hard)

    def remove_selected(self) -> None:
        rows = sorted({i.row() for i in self.grid.selectedIndexes()}, reverse=True)
        if not rows and self.grid.rowCount():
            rows = [self.grid.rowCount() - 1]
        for row in rows:
            self.grid.removeRow(row)
        self._update()

    # -- the profile --------------------------------------------------------------------------
    def profile(self) -> MpoProfile:
        """The scheme as entered; raises ValueError with a message for the user."""
        specs = []
        for row in range(self.grid.rowCount()):
            target_text = self._widget(row, TARGET).text().strip()
            try:
                target = float(target_text) if target_text else None
            except ValueError:
                raise ValueError(f"Row {row + 1}: the target must be a number.") from None
            specs.append(
                Desirability(
                    column=self._widget(row, COLUMN).currentText(),
                    shape=self._widget(row, SHAPE).currentData(),
                    low=self._widget(row, LOW).value(),
                    high=self._widget(row, HIGH).value(),
                    target=target,
                    width=self._widget(row, WIDTH).value(),
                    weight=self._widget(row, WEIGHT).value(),
                    hard=self._widget(row, HARD).isChecked(),
                )
            )
        if not specs:
            raise ValueError("Add at least one property.")
        if not any(s.weight > 0 for s in specs):
            raise ValueError("At least one property needs a weight above 0.")
        missing = sorted({s.column for s in specs} - set(self.columns))
        if missing:
            raise ValueError(f"Not number columns in this table: {', '.join(missing)}.")
        return MpoProfile(tuple(specs), self.method.currentData(), self.skip_missing.isChecked())

    def set_profile(self, profile: MpoProfile) -> None:
        self.grid.setRowCount(0)
        for spec in profile.properties:
            self.add_property(spec)
        self.method.setCurrentIndex(max(0, self.method.findData(profile.method)))
        self.skip_missing.setChecked(profile.skip_missing)
        self._update()

    def _update(self) -> None:
        for row in range(self.grid.rowCount()):
            shape = self._widget(row, SHAPE).currentData()
            bell = shape == "gaussian"
            for col in (LOW, HIGH, HARD):
                self._widget(row, col).setEnabled(not bell)
            self._widget(row, TARGET).setEnabled(shape in ("target", "gaussian"))
            self._widget(row, WIDTH).setEnabled(bell)
        try:
            self.profile()
        except ValueError as err:
            self.status.setText(str(err))
            self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(False)
            return
        self.status.setText("")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(True)

    def load_dialog(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Load MPO Profile", "", FILE_FILTER)
        if path:
            try:
                self.set_profile(MpoProfile.from_json(Path(path).read_text(encoding="utf-8")))
            except (OSError, ValueError) as err:
                QMessageBox.warning(self, "Load MPO Profile", str(err))

    def save_dialog(self) -> None:
        try:
            profile = self.profile()
        except ValueError as err:
            QMessageBox.warning(self, "Save MPO Profile", str(err))
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save MPO Profile", "", FILE_FILTER)
        if path:
            Path(path if path.endswith(".json") else path + ".json").write_text(
                profile.to_json(), encoding="utf-8"
            )

    # -- scoring ------------------------------------------------------------------------------
    def accept(self) -> None:
        try:
            profile = self.profile()
        except ValueError as err:
            self.status.setText(str(err))
            return
        name = self.output.text().strip() or DEFAULT_OUTPUT
        global _last, _last_output
        _last, _last_output = profile, name
        self.score(profile, name)
        super().accept()

    def score(self, profile: MpoProfile, name: str = DEFAULT_OUTPUT) -> dict[str, str]:
        ctx = self.ctx
        oids = ctx.target_oids(self.selected_only.isChecked())
        values = {s.column: numbers(ctx.table, s.column, oids) for s in profile.properties}
        overall, parts = profile.score(values)
        results = {name: (oids, NumberColumn.from_floats(overall, 3))}
        if self.per_property.isChecked():
            for spec, part in zip(profile.properties, parts, strict=True):
                results[f"{name}: {spec.column}"] = (oids, NumberColumn.from_floats(part, 3))
        names = ctx.apply_results(results, label="MPO score")
        scored = int(np.isfinite(overall).sum())
        ctx.show_message(f"Scored {scored:,} of {len(oids):,} row(s) in “{names[name]}”")
        return names
