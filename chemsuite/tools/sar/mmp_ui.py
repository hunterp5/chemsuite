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
from PySide6.QtWidgets import QComboBox, QSpinBox

from chemsuite.analysis.selectivity import looks_like_assay
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.wording import plural
from chemsuite.plotting.data import number_columns, numbers
from chemsuite.sar.mmp import Pair, Transform, find_pairs, transforms

from .common import RerunPanel, activity_combo, saved_rows, structure_cells
from .mmp_panel import NO_ACTIVITY, MmpPanel, panel_title, transform_rows

POOL_MIN_ROWS = 1000  # fragmenting is slow per molecule: worth the worker start-up early


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
            help="Transformations are summarized with the activity shift across their pairs.",
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
            lambda oids, activity, pairs, found, rows, cells: self.finish(
                oids, activity, recipe["activity"], pairs, found, recipe, rows, cells
            ),
        )

    def finish(
        self,
        oids,
        activity,
        name,
        pairs: list[Pair],
        found: list[Transform],
        recipe=None,
        rows: list[list] | None = None,
        cells: list | None = None,
    ):
        ctx = self.ctx
        self.panel = MmpPanel(ctx, oids, activity, name, found, recipe, rows, cells)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"Found {plural(len(pairs), 'matched pair')}; "
            f"{plural(len(found), 'transformation')} shown"
        )
        return self.panel


def _compute(ctx: ToolContext, recipe: dict, done):
    """Find the pairs *recipe* describes in the background; then
    ``done(oids, activity, pairs, found, rows, cells)`` (*cells*: the structures searched)."""
    oids = np.asarray(recipe["oids"], dtype=np.int64)
    cells = structure_cells(ctx, recipe["structure"], oids)
    name = recipe["activity"]
    activity = numbers(ctx.table, name, oids) if name != NO_ACTIVITY else np.full(len(oids), np.nan)
    max_variable, min_pairs = int(recipe["max_variable"]), int(recipe["min_pairs"])
    executor = ctx.workers_for(len(cells), POOL_MIN_ROWS)

    def work(job):
        pairs = find_pairs(
            cells,
            activity,
            max_variable=max_variable,
            executor=executor,
            progress=job.report,
            cancel=job.cancel_token,
        )
        found = transforms(pairs, min_pairs=min_pairs)
        return pairs, found, transform_rows(found)

    return ctx.run(
        "Finding matched pairs…", work, lambda result: done(oids, activity, *result, cells)
    )


def restore(ctx: ToolContext, state: dict) -> RerunPanel:
    """Reopen a saved panel: rerun the analysis on the rows it covered."""
    name = state["activity"]
    columns = [state["structure"]] + ([name] if name != NO_ACTIVITY else [])
    recipe = {**state, "oids": saved_rows(ctx, state, *columns)}
    holder = RerunPanel(MmpPanel.tool_id, panel_title(name), recipe)
    ctx.add_panel(holder, holder.title())

    def done(oids, activity, _pairs, found, rows, cells) -> None:
        holder.fill(MmpPanel(ctx, oids, activity, name, found, recipe, rows, cells))

    holder.job = _compute(ctx, recipe, done)
    return holder
