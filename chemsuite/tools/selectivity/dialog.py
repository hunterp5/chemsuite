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
"""The Selectivity dialog: a primary assay, its off-targets, and the narrow-window flag."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
)

from chemsuite.analysis.selectivity import (
    MIN_NAME,
    WINDOW_NAME,
    looks_like_assay,
    selectivity,
    selectivity_name,
)
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import number_columns, numbers

UNITS = {
    False: "p-values: pIC50, pKi… (higher is more potent)",
    True: "Concentrations: IC50, Ki… in one unit (lower is more potent)",
}


class SelectivityDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Selectivity")
        self.ctx = ctx
        self.columns = number_columns(ctx.table)
        self.primary = QComboBox()
        self.primary.addItems(self.columns)
        first = next((c for c in self.columns if looks_like_assay(c)), None)
        if first:
            self.primary.setCurrentText(first)
        self.off_targets = QListWidget()
        self.units = QComboBox()
        for key, title in UNITS.items():
            self.units.addItem(title, key)
        self.flag = QCheckBox("Flag a narrow window")
        self.flag.setChecked(True)
        self.cutoff = QDoubleSpinBox()
        self.cutoff.setRange(-1e9, 1e9)
        self.cutoff.setDecimals(3)
        self.cutoff.setValue(7.0)
        self.window = QDoubleSpinBox()
        self.window.setRange(0.0, 10.0)
        self.window.setDecimals(2)
        self.window.setValue(1.0)
        self.window.setSuffix(" log units")
        self.selected_only = QCheckBox(f"Selected rows only ({len(ctx.selection):,})")
        self.selected_only.setEnabled(len(ctx.selection) > 0)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Add Columns")

        window_row = QHBoxLayout()
        window_row.addWidget(QLabel("potent rows (primary at least"))
        window_row.addWidget(self.cutoff)
        window_row.addWidget(QLabel(") with selectivity below"))
        window_row.addWidget(self.window)
        form = QFormLayout()
        form.addRow("Primary assay:", self.primary)
        form.addRow("Off-targets:", self.off_targets)
        form.addRow("Values are:", self.units)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.flag)
        layout.addLayout(window_row)
        layout.addWidget(self.selected_only)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(620, 480)

        self.primary.currentIndexChanged.connect(lambda _i: self._fill_off_targets())
        self.off_targets.itemChanged.connect(lambda _item: self._update())
        self.units.currentIndexChanged.connect(self._units_changed)
        self.flag.toggled.connect(lambda _on: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._fill_off_targets()

    def _fill_off_targets(self) -> None:
        checked = set(self.chosen_off_targets())
        primary = self.primary.currentText()
        self.off_targets.blockSignals(True)
        self.off_targets.clear()
        for name in self.columns:
            if name == primary:
                continue
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            on = name in checked or (not checked and looks_like_assay(name))
            item.setCheckState(Qt.CheckState.Checked if on else Qt.CheckState.Unchecked)
            self.off_targets.addItem(item)
        self.off_targets.blockSignals(False)
        self._update()

    def _units_changed(self, *_args) -> None:
        concentrations = self.units.currentData()
        self.cutoff.setValue(100.0 if concentrations else 7.0)
        self.cutoff.setToolTip("At most this concentration" if concentrations else "")
        self._update()

    def set_off_targets(self, names: list[str]) -> None:
        wanted = set(names)
        for i in range(self.off_targets.count()):
            item = self.off_targets.item(i)
            item.setCheckState(
                Qt.CheckState.Checked if item.text() in wanted else Qt.CheckState.Unchecked
            )

    def chosen_off_targets(self) -> list[str]:
        items = (self.off_targets.item(i) for i in range(self.off_targets.count()))
        return [i.text() for i in items if i.checkState() == Qt.CheckState.Checked]

    def _update(self) -> None:
        for widget in (self.cutoff, self.window):
            widget.setEnabled(self.flag.isChecked())
        targets = self.chosen_off_targets()
        ok = bool(self.primary.currentText()) and bool(targets)
        if not self.columns:
            self.status.setText("Needs number columns with assay results.")
        elif not targets:
            self.status.setText("Tick at least one off-target.")
        else:
            self.status.setText("")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(ok)

    def accept(self) -> None:
        if not self.chosen_off_targets():
            return
        self.apply()
        super().accept()

    def apply(self) -> dict[str, str]:
        ctx = self.ctx
        table = ctx.table
        oids = ctx.target_oids(self.selected_only.isChecked())
        primary = self.primary.currentText()
        targets = self.chosen_off_targets()
        concentrations = bool(self.units.currentData())
        result = selectivity(
            numbers(table, primary, oids),
            {t: numbers(table, t, oids) for t in targets},
            concentrations=concentrations,
            potency_cutoff=self.cutoff.value(),
            window=self.window.value(),
        )
        results: dict = {
            selectivity_name(t): (oids, NumberColumn.from_floats(result.by_target[t], 2))
            for t in targets
        }
        if len(targets) > 1:
            results[MIN_NAME] = (oids, NumberColumn.from_floats(result.minimum, 2))
        if self.flag.isChecked():
            results[WINDOW_NAME] = (oids, result.narrow)
        names = ctx.apply_results(results, label="Selectivity")
        narrow = sum(1 for flag in result.narrow if flag == "yes")
        measured = int(np.isfinite(result.minimum).sum())
        message = f"Selectivity for {measured:,} row(s)"
        if self.flag.isChecked():
            message += f"; {narrow:,} in a narrow window"
        ctx.show_message(message)
        return names
