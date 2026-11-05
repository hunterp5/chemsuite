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
"""Pieces shared by the similarity tools: choosing features and rows, cached fingerprints."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QFormLayout,
    QListWidget,
    QListWidgetItem,
    QRadioButton,
    QWidget,
)

from chemsuite.analysis.fingerprints import (
    DEFAULT_KIND,
    KINDS,
    METRICS,
    Fingerprints,
)
from chemsuite.analysis.fingerprints import fingerprint_task as cached_fingerprint_task
from chemsuite.app.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers, shown_oids

SCOPES = {"all": "All rows", "shown": "Shown rows", "selected": "Selected rows"}
LABEL_COLUMNS = ("ID", "Name", "Title", "Compound", "Compound ID")


def scope_combo(ctx: ToolContext, default: str = "all") -> QComboBox:
    combo = QComboBox()
    for key, title in SCOPES.items():
        if key == "selected":
            title += f" ({len(ctx.selection):,})"
        combo.addItem(title, key)
    combo.setCurrentIndex(max(0, combo.findData(default)))
    return combo


def scope_oids(ctx: ToolContext, scope: str) -> np.ndarray:
    if scope == "selected":
        return ctx.selection.as_array()
    if scope == "shown":
        return shown_oids(ctx.table)
    return ctx.table.oids.copy()


def label_column(ctx: ToolContext) -> str | None:
    """A column that names compounds (ID, Name…), to say which compound is meant."""
    table = ctx.table
    return next((n for n in LABEL_COLUMNS if n in table and table.kind(n) != "structure"), None)


def labels(ctx: ToolContext, oids: np.ndarray) -> list[str]:
    table = ctx.table
    name = label_column(ctx)
    if name is not None:
        return table.column(name).displays(table.slots_of(oids))
    position = np.full(len(table), -1, dtype=np.int64)
    view = table.view()
    position[view] = np.arange(len(view))
    return [f"Row {p + 1}" if p >= 0 else "(hidden row)" for p in position[table.slots_of(oids)]]


def fingerprint_task(
    ctx: ToolContext, column: str, key: str, oids: np.ndarray
) -> tuple[Callable[[object], Fingerprints], Callable[[Fingerprints], None]]:
    """``(compute, remember)``: *compute(job)* runs in the job (cache hit: instant);
    *remember(fps)* stores the result on the GUI thread once the job is done."""
    return cached_fingerprint_task(ctx.table, column, key, oids, ctx.process_pool)


class FeatureChooser(QWidget):
    """Fingerprints of a structure column, or a set of number columns."""

    changed = Signal()

    def __init__(self, ctx: ToolContext, *, metric: bool = True, numbers_allowed: bool = True):
        super().__init__()
        self.ctx = ctx
        self.use_fingerprints = QRadioButton("Fingerprints")
        self.use_numbers = QRadioButton("Number columns")
        group = QButtonGroup(self)
        group.addButton(self.use_fingerprints)
        group.addButton(self.use_numbers)
        self.use_fingerprints.setChecked(True)
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        self.kind = QComboBox()
        for kind in (k for k in KINDS.values() if k.listed):
            self.kind.addItem(kind.title, kind.key)
            if kind.help:
                self.kind.setItemData(self.kind.count() - 1, kind.help, Qt.ItemDataRole.ToolTipRole)
        self.kind.setCurrentIndex(self.kind.findData(DEFAULT_KIND))
        self.metric = QComboBox()
        for key, title in METRICS.items():
            self.metric.addItem(title, key)
        self.columns = QListWidget()
        self.columns.setMaximumHeight(140)
        for name in number_columns(ctx.table):
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked)
            self.columns.addItem(item)

        self.form = QFormLayout(self)
        self.form.setContentsMargins(0, 0, 0, 0)
        if numbers_allowed:
            self.form.addRow("Compare by:", self.use_fingerprints)
            self.form.addRow("", self.use_numbers)
        self.form.addRow("Structures:", self.structure)
        self.form.addRow("Fingerprint:", self.kind)
        self.form.addRow("Similarity:", self.metric)
        self.form.addRow("Columns:", self.columns)
        self._show_metric = metric
        if not ctx.structure_columns():
            self.use_numbers.setChecked(True)
            self.use_fingerprints.setEnabled(False)
        for signal in (
            self.use_fingerprints.toggled,
            self.structure.currentIndexChanged,
            self.kind.currentIndexChanged,
            self.metric.currentIndexChanged,
            self.columns.itemChanged,
        ):
            signal.connect(lambda *_: self._changed())
        self._changed()

    def _changed(self) -> None:
        fingerprints = self.fingerprints()
        self.form.setRowVisible(self.structure, fingerprints)
        self.form.setRowVisible(self.kind, fingerprints)
        self.form.setRowVisible(self.metric, fingerprints and self._show_metric)
        self.form.setRowVisible(self.columns, not fingerprints)
        self.changed.emit()

    def fingerprints(self) -> bool:
        return self.use_fingerprints.isChecked()

    def set_number_columns(self, names) -> None:
        wanted = set(names)
        for i in range(self.columns.count()):
            item = self.columns.item(i)
            item.setCheckState(
                Qt.CheckState.Checked if item.text() in wanted else Qt.CheckState.Unchecked
            )

    def number_columns(self) -> list[str]:
        items = (self.columns.item(i) for i in range(self.columns.count()))
        return [i.text() for i in items if i.checkState() == Qt.CheckState.Checked]

    def problem(self) -> str:
        """Why the features cannot be used yet ("" when they can)."""
        if self.fingerprints():
            return "" if self.structure.currentText() else "There is no structure column."
        return "" if self.number_columns() else "Tick at least one number column."

    def describe(self) -> str:
        if self.fingerprints():
            kind = self.kind.currentText().split(",")[0]
            return f"{kind} fingerprints of {self.structure.currentText()}"
        return ", ".join(self.number_columns())

    def number_values(self, oids: np.ndarray) -> list[np.ndarray]:
        return [numbers(self.ctx.table, n, oids) for n in self.number_columns()]
