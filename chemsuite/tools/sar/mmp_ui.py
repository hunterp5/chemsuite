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
"""Matched Molecular Pairs…: single-change pairs grouped by transformation."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.selectivity import looks_like_assay
from chemsuite.app.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers
from chemsuite.sar.mmp import Pair, Transform, find_pairs, readable, transforms
from chemsuite.tools.similarity.common import labels, scope_combo, scope_oids

from .common import ResultTable, activity_combo, structure_cells

WARN_ROWS = 50_000
POOL_MIN_ROWS = 1000  # fragmenting is slow per molecule: worth the worker start-up early
NO_ACTIVITY = "(none)"


class MmpDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Matched Molecular Pairs")
        self.ctx = ctx
        self.job = None
        self.panel: MmpPanel | None = None
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        self.activity = QComboBox()
        self.activity.addItem(NO_ACTIVITY)
        self.activity.addItems(number_columns(ctx.table))

        guess = activity_combo(ctx).currentText()
        if guess and looks_like_assay(guess):
            self.activity.setCurrentText(guess)
        self.rows = scope_combo(ctx, "shown")
        self.max_variable = QSpinBox()
        self.max_variable.setRange(1, 40)
        self.max_variable.setValue(10)
        self.max_variable.setToolTip("Largest changing part, in heavy atoms.")
        self.min_pairs = QSpinBox()
        self.min_pairs.setRange(1, 1000)
        self.min_pairs.setValue(2)
        self.min_pairs.setToolTip("Show transformations seen in at least this many pairs.")
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Find Pairs")
        form = QFormLayout(self)
        form.addRow("Structures:", self.structure)
        form.addRow("Activity:", self.activity)
        form.addRow("Rows:", self.rows)
        form.addRow("Largest change (heavy atoms):", self.max_variable)
        form.addRow("Minimum pairs:", self.min_pairs)
        form.addRow(self.status)
        form.addRow(self.buttons)
        self.resize(480, 0)
        self.rows.currentIndexChanged.connect(lambda *_: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    def _update(self) -> None:
        count = len(scope_oids(self.ctx, self.rows.currentData()))
        problem = "" if self.structure.currentText() else "There is no structure column."
        note = problem
        if not problem and count > WARN_ROWS:
            note = f"{count:,} rows: fragmenting may take several minutes (it can be cancelled)."
        self.status.setText(note)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if not self.structure.currentText():
            return
        self.job = self.find()
        super().accept()

    def find(self):
        ctx = self.ctx
        oids = scope_oids(ctx, self.rows.currentData())
        cells = structure_cells(ctx, self.structure.currentText(), oids)
        name = self.activity.currentText()
        activity = (
            numbers(ctx.table, name, oids) if name != NO_ACTIVITY else np.full(len(oids), np.nan)
        )
        max_variable, min_pairs = self.max_variable.value(), self.min_pairs.value()
        executor = ctx.process_pool() if len(cells) >= POOL_MIN_ROWS else None

        def work(job):
            pairs = find_pairs(
                cells,
                activity,
                max_variable=max_variable,
                executor=executor,
                progress=job.report,
                cancel=job.cancel_token,
            )
            return pairs, transforms(pairs, min_pairs=min_pairs)

        def done(result) -> None:
            pairs, found = result
            self.finish(oids, activity, name, pairs, found)

        return ctx.run("Finding matched pairs…", work, done)

    def finish(self, oids, activity, name, pairs: list[Pair], found: list[Transform]):
        ctx = self.ctx
        self.panel = MmpPanel(ctx, oids, activity, name, found)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"{len(pairs):,} matched pair(s); {len(found):,} transformation(s) shown", 10_000
        )
        return self.panel


class MmpPanel(QWidget):
    """Transformations, and the pairs behind the chosen one."""

    titleChanged = Signal(str)
    tool_id = "mmp"

    def __init__(self, ctx: ToolContext, oids, activity, name: str, found: list[Transform]):
        super().__init__()
        self.ctx = ctx
        self.oids = np.asarray(oids)
        self.activity = np.asarray(activity)
        self.name = name
        self.found = found
        self.names = labels(ctx, self.oids)
        self.transforms = ResultTable(["From", "To", "Pairs", "Mean Δ", "SD", "Up", "Down"])
        self.transforms.fill(
            [
                [
                    readable(t.from_group),
                    readable(t.to_group),
                    t.count,
                    t.mean,
                    t.sd,
                    t.ups,
                    t.downs,
                ]
                for t in found
            ],
            found,
        )
        self.transforms.setToolTip(
            f"Δ = {name} of the “to” compound minus the “from” compound."
            if name != NO_ACTIVITY
            else "No activity column was chosen."
        )
        self.pairs = ResultTable(["From compound", "To compound", "From", "To", "Δ"])
        self.transforms.itemSelectionChanged.connect(self.show_pairs)
        self.transforms.sortByColumn(2, Qt.SortOrder.DescendingOrder)
        select = QPushButton("Select Compounds")
        select.setToolTip("Select the compounds in the chosen transformations or pairs.")
        select.clicked.connect(lambda: self.select_rows())
        copy = QPushButton("Copy Transformations")
        copy.clicked.connect(lambda: self.transforms.copy_text())
        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.transforms)
        split.addWidget(self.pairs)
        buttons = QHBoxLayout()
        buttons.addWidget(select)
        buttons.addStretch(1)
        buttons.addWidget(copy)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(split, 1)
        layout.addLayout(buttons)

    def title(self) -> str:
        return "Matched Pairs" + (f": {self.name}" if self.name != NO_ACTIVITY else "")

    def show_pairs(self) -> None:
        chosen = [t for t in self.transforms.selected_data() if t is not None]
        pairs = [p for t in chosen for p in t.pairs]
        self.pairs.fill(
            [
                [self.names[p.a], self.names[p.b], self.activity[p.a], self.activity[p.b], p.delta]
                for p in pairs
            ],
            pairs,
        )
        heads = ["From compound", "To compound", f"From {self.name}", f"To {self.name}", "Δ"]
        self.pairs.setHorizontalHeaderLabels(heads)

    def select_rows(self) -> list[int]:
        pairs = [p for p in self.pairs.selected_data() if p is not None]
        if pairs:
            rows = {r for p in pairs for r in (p.a, p.b)}
        else:
            chosen = [t for t in self.transforms.selected_data() if t is not None]
            rows = set().union(*(t.rows() for t in chosen)) if chosen else set()
        oids = self.oids[sorted(rows)].tolist()
        if oids:
            self.ctx.selection.set(oids, "mmp")
        return oids

    def state(self) -> dict:
        return {}

    def set_state(self, state: dict) -> None:
        pass

    def close_panel(self) -> None:
        pass
