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
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QListWidget, QListWidgetItem

from chemsuite.analysis.selectivity import (
    MIN_NAME,
    WINDOW_NAME,
    looks_like_assay,
    selectivity,
    selectivity_name,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural
from chemsuite.plotting.data import number_columns, numbers


class SelectivityDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Selectivity", guide="selectivity")
        self.ctx = ctx
        self.columns = number_columns(ctx.table)

        inputs = self.section("Input")
        self.assay = QComboBox()
        self.assay.addItems(self.columns)
        first = next((c for c in self.columns if looks_like_assay(c)), None)
        if first:
            self.assay.setCurrentText(first)
        inputs.row("Primary assay", self.assay)
        self.off_targets = QListWidget()
        self.off_targets.setMaximumHeight(120)
        inputs.row("Off-targets", self.off_targets, help="Tick the off-target assays.")
        self.units = Segmented((("p", "p-Values"), ("concentration", "Concentrations")))
        inputs.row(
            "Values are",
            self.units,
            help="p-Values: pIC50, pKi… (higher is more potent). Concentrations: IC50, Ki… "
            "in one unit (lower is more potent).",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        narrow = self.section("Narrow Window")
        self.flag = QCheckBox("Flag potent rows with a narrow window")
        self.flag.setChecked(True)
        narrow.row("Flag", self.flag, help=f"Writes a “{WINDOW_NAME}” column: yes or no.")
        self.cutoff = QDoubleSpinBox()
        self.cutoff.setRange(-1e9, 1e9)
        self.cutoff.setDecimals(3)
        self.cutoff.setValue(7.0)
        narrow.row("Potent at", self.cutoff)
        self.window = QDoubleSpinBox()
        self.window.setRange(0.0, 10.0)
        self.window.setDecimals(2)
        self.window.setValue(1.0)
        self.window.setSuffix(" log units")
        narrow.row(
            "Window below",
            self.window,
            help="A potent row is flagged when its selectivity over any off-target is below this.",
        )
        self._narrow = narrow

        self.assay.currentIndexChanged.connect(lambda _i: self._fill_off_targets())
        self.off_targets.itemChanged.connect(lambda _item: self._update())
        self.units.changed.connect(self._units_changed)
        self.flag.toggled.connect(lambda _on: self._update())
        self.rows.changed.connect(lambda *_: self._update())
        self._units_changed()
        self._fill_off_targets()

    def concentrations(self) -> bool:
        return self.units.value() == "concentration"

    def _fill_off_targets(self) -> None:
        checked = set(self.chosen_off_targets())
        primary = self.assay.currentText()
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
        concentrations = self.concentrations()
        self.cutoff.setValue(100.0 if concentrations else 7.0)
        self._narrow.set_help(
            self.cutoff,
            "Rows at most this concentration count as potent."
            if concentrations
            else "Rows with the primary assay at least this count as potent.",
        )
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
        self.refresh()

    def problem(self) -> str:
        if not self.columns:
            return "Needs number columns with assay results."
        if not self.chosen_off_targets():
            return "Tick at least one off-target."
        return ""

    def verb(self) -> str:
        return f"Calculate {len(self.rows.oids()):,} Rows"

    def run(self):
        self.apply()

    def apply(self) -> dict[str, str]:
        ctx = self.ctx
        table = ctx.table
        oids = self.rows.oids()
        primary = self.assay.currentText()
        targets = self.chosen_off_targets()
        concentrations = self.concentrations()
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
        message = f"Calculated selectivity for {plural(measured, 'row')}"
        if self.flag.isChecked():
            message += f"; {narrow:,} in a narrow window"
        ctx.show_message(message)
        return names
