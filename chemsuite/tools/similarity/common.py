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
    QCheckBox,
    QComboBox,
    QListWidget,
    QListWidgetItem,
    QRadioButton,
    QSpinBox,
    QWidget,
)

from chemsuite.analysis.clustering import SCALINGS
from chemsuite.analysis.fingerprints import (
    DEFAULT_KIND,
    KINDS,
    METRICS,
    RADII,
    RADIUS,
    SIZED,
    SIZES,
    Fingerprints,
    variant,
)
from chemsuite.analysis.fingerprints import fingerprint_task as cached_fingerprint_task
from chemsuite.app.document.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers

LABEL_COLUMNS = ("ID", "Name", "Title", "Compound", "Compound ID")


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

    Its rows go into a dialog *section* (see ``chemsuite.app.widgets.dialogs``), sharing the
    dialog's label column, and "Compare by" is a segmented choice.

    With *both* there is no "Compare by": the Fingerprint field (which has a "None" choice) and
    the number columns are always shown, and any combination of the two is used.

    A *more* section (with *both*) adds how the features are made: the fingerprint's size,
    radius and counting, the fingerprints' share when columns are used too, and how the
    columns are scaled.
    """

    changed = Signal()

    def __init__(
        self,
        ctx: ToolContext,
        *,
        metric: bool = True,
        numbers_allowed: bool = True,
        section,
        both: bool = False,
        more=None,
    ):
        super().__init__()
        self.ctx = ctx
        self.section = section
        self._free = both  # fingerprints and columns chosen freely
        if self._free:
            self.compare = None
            # never shown: kept so callers can still ask which is chosen
            self.use_fingerprints = QRadioButton()
            self.use_numbers = QRadioButton()
        else:
            from chemsuite.app.widgets.dialogs import Segmented

            self.compare = Segmented(
                (("fingerprints", "Fingerprints"), ("numbers", "Number columns"))
            )
            self.use_fingerprints = self.compare.buttons["fingerprints"]
            self.use_numbers = self.compare.buttons["numbers"]
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
        self._show_metric = metric
        self.more = more
        self._build_more(more)
        if not ctx.structure_columns():
            if self._free:
                self.kind.setCurrentIndex(self.kind.findData(None))
                self.kind.setEnabled(False)
            self.use_numbers.setChecked(True)
            self.use_fingerprints.setEnabled(False)
        toggles = [self.use_fingerprints.toggled, self.use_numbers.toggled]
        for signal in (
            *toggles,
            self.structure.currentIndexChanged,
            self.kind.currentIndexChanged,
            self.metric.currentIndexChanged,
            self.columns.itemChanged,
        ):
            signal.connect(lambda *_: self._changed())
        for signal in (
            self.size.currentIndexChanged,
            self.radius.valueChanged,
            self.count.toggled,
            self.weight.valueChanged,
            self.scaling.currentIndexChanged,
        ):
            signal.connect(lambda *_: self.changed.emit())
        self._changed()

    def _build_more(self, more) -> None:
        """How the features are made (see the class notes); the rows go in *more*."""
        self.size = QComboBox()
        for bits in SIZES:
            self.size.addItem(f"{bits:,} bits", bits)
        self.radius = QSpinBox()
        self.radius.setRange(RADII[0], RADII[-1])
        self.count = QCheckBox("Count each feature")
        self.weight = QSpinBox()
        self.weight.setRange(0, 100)
        self.weight.setSuffix(" %")
        self.weight.setValue(50)
        self.scaling = QComboBox()
        for key, title in SCALINGS.items():
            self.scaling.addItem(title, key)
        self._kind_shown = None  # the kind whose size and radius the boxes show
        if more is not None:
            more.row(
                "Fingerprint size",
                self.size,
                help="Bits the features are folded onto: more bits, fewer collisions.",
            )
            more.row(
                "Radius", self.radius, help="Bonds out from each atom: 2 is ECFP4, 3 is ECFP6."
            )
            more.row(
                "Counts",
                self.count,
                help="How often each feature occurs, not only whether it does. Counts are "
                "compared by Euclidean-type distances, not Jaccard.",
            )
            more.row(
                "Fingerprint weight",
                self.weight,
                help="The fingerprints' share when number columns are ticked too; the "
                "columns get the rest.",
            )
            more.row("Column scaling", self.scaling,
                     help="How the number columns are put on one scale. Z-score: mean 0, "
                     "standard deviation 1. Min-max: 0 to 1. Robust: median 0, interquartile "
                     "range 1, so outliers weigh less. None: as they are.")  # fmt: skip

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
        kind = KINDS.get(self.kind.currentData())
        self.section.set_help(self.kind, kind.help if kind is not None else "")
        if self.more is not None:
            self._show_more(kind)
        self.changed.emit()

    def _show_more(self, kind) -> None:
        key = kind.key if kind is not None else None
        if key != self._kind_shown:  # a new kind starts from its own size and radius
            self._kind_shown = key
            if kind is not None:
                self.size.setCurrentIndex(max(0, self.size.findData(kind.bits)))
                self.radius.setValue(RADIUS.get(key, 2))
        more = self.more
        bits = self.uses_fingerprints()
        more.set_row_visible(self.size, bits and key in SIZED)
        more.set_row_visible(self.radius, bits and key in RADIUS)
        more.set_row_visible(self.count, bits and key in SIZED)
        more.set_row_visible(self.weight, self.both())
        more.set_row_visible(self.scaling, self.uses_columns() and bool(self.number_columns()))

    def uses_fingerprints(self) -> bool:
        """Fingerprints are part of the features (alone or with columns)."""
        return self.fingerprints() or self.both()

    def uses_columns(self) -> bool:
        """Number columns are part of the features (alone or with fingerprints)."""
        return self.both() or not self.fingerprints()

    def _show(self, field: QWidget, visible: bool) -> None:
        self.section.set_row_visible(field, visible)

    def fingerprints(self) -> bool:
        """Fingerprints alone (not with number columns: see ``both``)."""
        if self._free:
            return self.kind.currentData() is not None and not self.number_columns()
        return self.use_fingerprints.isChecked()

    def both(self) -> bool:
        if self._free:
            return self.kind.currentData() is not None and bool(self.number_columns())
        return False

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

    def fingerprint_key(self) -> str | None:
        """The fingerprint with its chosen size and radius (never counts: see
        :meth:`counts`); None for none."""
        base = self.kind.currentData()
        if base is None or not self.uses_fingerprints():
            return None
        if self.more is None or base not in SIZED:
            return base
        radius = self.radius.value() if base in RADIUS else None
        return variant(base, self.size.currentData(), radius)

    def counts(self) -> bool:
        """Count fingerprints rather than bits."""
        return self.more is not None and self.count.isChecked() and self.kind.currentData() in SIZED

    def bits(self) -> bool:
        """Fingerprint bits alone: what Jaccard-type distances need."""
        return self.fingerprints() and not self.counts()

    def fingerprint_weight(self) -> float:
        return self.weight.value() / 100.0

    def column_scaling(self) -> str:
        return self.scaling.currentData()

    def remembered(self) -> dict[str, QWidget]:
        """The fields a dialog remembers (see ``ToolDialog.remember``): not the structure
        column or the ticked columns, which belong to one table."""
        fields = {"fingerprint": self.kind}
        if self.more is not None:
            fields |= {"size": self.size, "radius": self.radius, "counts": self.count,
                       "weight": self.weight, "scaling": self.scaling}  # fmt: skip
        return fields

    def number_values(self, oids: np.ndarray) -> list[np.ndarray]:
        return [numbers(self.ctx.table, n, oids) for n in self.number_columns()]
