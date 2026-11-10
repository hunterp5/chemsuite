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
    QHBoxLayout,
    QPushButton,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.selectivity import looks_like_assay
from chemsuite.app.dialogs import RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers
from chemsuite.sar.mmp import Pair, Transform, find_pairs, readable, transforms
from chemsuite.tools.similarity.common import labels

from .common import RerunPanel, ResultTable, activity_combo, saved_rows, structure_cells

POOL_MIN_ROWS = 1000  # fragmenting is slow per molecule: worth the worker start-up early
NO_ACTIVITY = "(none)"


class MmpDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Matched Molecular Pairs", guide="matched-molecular-pairs")
        self.ctx = ctx
        self.panel: MmpPanel | None = None
        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.activity = QComboBox()
        self.activity.addItem(NO_ACTIVITY)
        self.activity.addItems(number_columns(ctx.table))
        guess = activity_combo(ctx).currentText()
        if guess and looks_like_assay(guess):
            self.activity.setCurrentText(guess)
        inputs.row(
            "Activity",
            self.activity,
            help="Transformations are summarised with the activity shift across their pairs.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx, "shown"))

        method = self.section("Method")
        self.max_variable = QSpinBox()
        self.max_variable.setRange(1, 40)
        self.max_variable.setValue(10)
        self.max_variable.setSuffix(" heavy atoms")
        method.row("Largest change", self.max_variable, help="The largest changing part.")
        self.min_pairs = QSpinBox()
        self.min_pairs.setRange(1, 1000)
        self.min_pairs.setValue(2)
        method.row(
            "Min. pairs",
            self.min_pairs,
            help="List transformations seen in at least this many pairs.",
        )

        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        return "" if self.structure.currentText() else "There is no structure column."

    def verb(self) -> str:
        return f"Find Pairs in {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.find()

    def recipe(self) -> dict:
        """Everything the analysis needs, saved with its panel to rerun it on reopening."""
        return {
            "oids": self.rows.oids().tolist(),
            "structure": self.structure.currentText(),
            "activity": self.activity.currentText(),
            "max_variable": self.max_variable.value(),
            "min_pairs": self.min_pairs.value(),
        }

    def find(self):
        recipe = self.recipe()
        return _compute(
            self.ctx,
            recipe,
            lambda oids, activity, pairs, found: self.finish(
                oids, activity, recipe["activity"], pairs, found, recipe
            ),
        )

    def finish(self, oids, activity, name, pairs: list[Pair], found: list[Transform], recipe=None):
        ctx = self.ctx
        self.panel = MmpPanel(ctx, oids, activity, name, found, recipe)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"{len(pairs):,} matched pair(s); {len(found):,} transformation(s) shown", 10_000
        )
        return self.panel


class MmpPanel(QWidget):
    """Transformations, and the pairs behind the chosen one."""

    titleChanged = Signal(str)
    tool_id = "mmp"

    def __init__(
        self, ctx: ToolContext, oids, activity, name: str, found: list[Transform], recipe=None
    ):
        super().__init__()
        self.ctx = ctx
        self.recipe = recipe or {}
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
        return _title(self.name)

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
        return dict(self.recipe)

    def close_panel(self) -> None:
        pass


def _compute(ctx: ToolContext, recipe: dict, done):
    """Find the pairs *recipe* describes in the background; then
    ``done(oids, activity, pairs, found)``."""
    oids = np.asarray(recipe["oids"], dtype=np.int64)
    cells = structure_cells(ctx, recipe["structure"], oids)
    name = recipe["activity"]
    activity = numbers(ctx.table, name, oids) if name != NO_ACTIVITY else np.full(len(oids), np.nan)
    max_variable, min_pairs = int(recipe["max_variable"]), int(recipe["min_pairs"])
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

    return ctx.run("Finding matched pairs…", work, lambda result: done(oids, activity, *result))


def restore(ctx: ToolContext, state: dict) -> RerunPanel:
    """Reopen a saved panel: rerun the analysis on the rows it covered."""
    name = state["activity"]
    columns = [state["structure"]] + ([name] if name != NO_ACTIVITY else [])
    recipe = {**state, "oids": saved_rows(ctx, state, *columns)}
    holder = RerunPanel(MmpPanel.tool_id, _title(name), recipe)
    ctx.add_panel(holder, holder.title())

    def done(oids, activity, _pairs, found) -> None:
        holder.fill(MmpPanel(ctx, oids, activity, name, found, recipe))

    holder.job = _compute(ctx, recipe, done)
    return holder


def _title(activity: str) -> str:
    return "Matched Pairs" + (f": {activity}" if activity != NO_ACTIVITY else "")
