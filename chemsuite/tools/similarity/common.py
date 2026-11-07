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
from chemsuite.app.dialogs import follow_rows
from chemsuite.app.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers, shown_oids

SCOPES = {"all": "All rows", "shown": "Shown rows", "selected": "Selected rows"}
LABEL_COLUMNS = ("ID", "Name", "Title", "Compound", "Compound ID")


def scope_combo(ctx: ToolContext, default: str = "all") -> QComboBox:
    combo = QComboBox()
    for key, title in SCOPES.items():
        combo.addItem(title, key)
    combo.setCurrentIndex(max(0, combo.findData(default)))
    selected = combo.findData("selected")
    follow_rows(
        combo,
        ctx,
        lambda: combo.setItemText(selected, f"{SCOPES['selected']} ({len(ctx.selection):,})"),
    )
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
    """Fingerprints of a structure column, a set of number columns, or (with *both*) the two
    together.

    Given a dialog *section* (see ``chemsuite.app.dialogs``), its rows go into that section,
    sharing the dialog's label column, and "Compare by" is a segmented choice. Without one it
    lays its rows out in its own form (dialogs not yet on the dialog kit).

    With *both* there is no "Compare by": the Fingerprint field (which has a "None" choice) and
    the number columns are always shown, and any combination of the two is used.
    """

    changed = Signal()

    def __init__(
        self,
        ctx: ToolContext,
        *,
        metric: bool = True,
        numbers_allowed: bool = True,
        section=None,
        both: bool = False,
    ):
        super().__init__()
        self.ctx = ctx
        self.section = section
        self._free = both and section is not None  # fingerprints and columns chosen freely
        if self._free:
            self.compare = None
            self.use_both = None
            # never shown: kept so callers can still ask which is chosen
            self.use_fingerprints = QRadioButton()
            self.use_numbers = QRadioButton()
        elif section is not None:
            from chemsuite.app.dialogs import Segmented

            choices = [
                ("fingerprints", "Fingerprints"),
                ("numbers", "Numbers" if both else "Number columns"),
            ]
            if both:
                choices.append(("both", "Both"))
            self.compare = Segmented(tuple(choices))
            self.use_fingerprints = self.compare.buttons["fingerprints"]
            self.use_numbers = self.compare.buttons["numbers"]
            self.use_both = self.compare.buttons.get("both")
        else:
            self.compare = None
            self.use_both = None
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
        if self._free:
            self.kind.addItem("None", None)
            self.kind.setItemData(
                self.kind.count() - 1,
                "No fingerprints: only the ticked number columns.",
                Qt.ItemDataRole.ToolTipRole,
            )
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

        if section is not None:
            self.form = None
            if numbers_allowed and not self._free:
                section.row("Compare by", self.compare)
            section.row("Structures", self.structure)
            section.row("Fingerprint", self.kind)
            section.row("Similarity", self.metric)
            section.row(
                "Columns",
                self.columns,
                help="Tick the number columns to use, with or without fingerprints."
                if self._free
                else "Tick the number columns to compare by.",
            )
            self.columns.setMaximumHeight(120)
        else:
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
            if self._free:
                self.kind.setCurrentIndex(self.kind.findData(None))
                self.kind.setEnabled(False)
            self.use_numbers.setChecked(True)
            self.use_fingerprints.setEnabled(False)
            if self.use_both is not None:
                self.use_both.setEnabled(False)
        toggles = [self.use_fingerprints.toggled, self.use_numbers.toggled]
        if self.use_both is not None:
            toggles.append(self.use_both.toggled)
        for signal in (
            *toggles,
            self.structure.currentIndexChanged,
            self.kind.currentIndexChanged,
            self.metric.currentIndexChanged,
            self.columns.itemChanged,
        ):
            signal.connect(lambda *_: self._changed())
        self._changed()

    def _changed(self) -> None:
        if self._free:
            self._show(self.structure, self.kind.currentData() is not None)
            self._show(self.metric, False)
        else:
            bits = self.fingerprints() or self.both()
            self._show(self.structure, bits)
            self._show(self.kind, bits)
            self._show(self.metric, self.fingerprints() and self._show_metric)
            self._show(self.columns, not self.fingerprints())
        if self.section is not None:
            kind = KINDS.get(self.kind.currentData())
            self.section.set_help(self.kind, kind.help if kind is not None else "")
        self.changed.emit()

    def _show(self, field: QWidget, visible: bool) -> None:
        if self.section is not None:
            self.section.set_row_visible(field, visible)
        else:
            self.form.setRowVisible(field, visible)

    def fingerprints(self) -> bool:
        """Fingerprints alone (not with number columns: see ``both``)."""
        if self._free:
            return self.kind.currentData() is not None and not self.number_columns()
        return self.use_fingerprints.isChecked()

    def both(self) -> bool:
        if self._free:
            return self.kind.currentData() is not None and bool(self.number_columns())
        return self.use_both is not None and self.use_both.isChecked()

    def set_fingerprint(self, key: str | None) -> None:
        """Choose the fingerprint kind (None: no fingerprints; only with *both*)."""
        self.kind.setCurrentIndex(self.kind.findData(key))

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
        if self._free and self.kind.currentData() is None and not self.number_columns():
            return "Choose a fingerprint or tick at least one number column."
        if self.fingerprints() or self.both():
            if not self.structure.currentText():
                return "There is no structure column."
            if self.fingerprints():
                return ""
        return "" if self.number_columns() else "Tick at least one number column."

    def describe(self) -> str:
        kind = self.kind.currentText().split(",")[0]
        bits = f"{kind} fingerprints of {self.structure.currentText()}"
        if self.fingerprints():
            return bits
        if self.both():
            return f"{bits} and {', '.join(self.number_columns())}"
        return ", ".join(self.number_columns())

    def number_values(self, oids: np.ndarray) -> list[np.ndarray]:
        return [numbers(self.ctx.table, n, oids) for n in self.number_columns()]
