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
"""Diverse Subset: pick compounds that cover the chemistry, by MaxMin."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QLineEdit, QSpinBox

from chemsuite.analysis.diversity import Picks, maxmin_fingerprints, maxmin_numeric
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog, selected_only
from chemsuite.core.jobs import Cancelled
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural

from .common import FeatureChooser, fingerprint_task

RANK_COLUMN = "Diversity Rank"


class DiverseSubsetDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Diverse Subset", guide="diverse-subset")
        self.ctx = ctx

        inputs = self.section("Input")
        self.features = FeatureChooser(ctx, section=inputs)
        self.rows = inputs.row(
            "Pick from", RowsChoice(ctx, "shown"), help="The rows the picks are taken from."
        )

        method = self.section("Method")
        self.count = QSpinBox()
        self.count.setRange(1, 1_000_000)
        self.count.setValue(min(50, max(1, len(ctx.table))))
        method.row(
            "How many",
            self.count,
            help="MaxMin: each pick is the compound farthest from all the earlier picks.",
        )
        self.seeded = selected_only(ctx, "Add to the selected rows")
        method.row(
            "Start from",
            self.seeded,
            help="Pick the compounds the selected rows do not already cover.",
        )

        output = self.section("Output")
        self.output = QLineEdit(RANK_COLUMN)
        output.row("Column name", self.output, help="Each pick's rank: 1 is the first pick.")
        self.select = QCheckBox("The picked rows")
        self.select.setChecked(True)
        output.row("Select", self.select)

        for signal in (self.features.changed, self.rows.changed, self.count.valueChanged):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        if not len(self.rows.oids()):
            return "There are no rows to pick from."
        return ""

    def verb(self) -> str:
        return f"Pick {self.count.value():,} of {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.pick()

    def pick(self):
        ctx = self.ctx
        pool = self.rows.oids()
        seeds = ctx.selection.as_array() if self.seeded.isChecked() else np.empty(0, np.int64)
        everything = np.union1d(pool, seeds)
        seed_rows = np.searchsorted(everything, seeds)
        in_pool = np.isin(everything, pool)
        count = self.count.value()
        name = self.output.text().strip() or RANK_COLUMN

        if self.features.fingerprints():
            metric = self.features.metric.currentData()
            compute, remember = fingerprint_task(
                ctx,
                self.features.structure.currentText(),
                self.features.kind.currentData(),
                everything,
            )

            def work(job):
                fps = compute(job)
                try:
                    return fps, maxmin_fingerprints(
                        _only_pool(fps, in_pool, seed_rows),
                        count,
                        seeds=seed_rows,
                        metric=metric,
                        progress=job.report,
                        cancel=job.cancel_token,
                    )
                except Cancelled as stopped:  # keep the picks made so far
                    raise Cancelled(stopped.partial and (fps, stopped.partial)) from None
        else:
            values = self.features.number_values(everything)
            remember = None

            def work(job):
                masked = [
                    np.where(in_pool | _mask(len(everything), seed_rows), v, np.nan) for v in values
                ]
                try:
                    return None, maxmin_numeric(
                        masked, count, seeds=seed_rows, progress=job.report, cancel=job.cancel_token
                    )
                except Cancelled as stopped:
                    raise Cancelled(stopped.partial and (None, stopped.partial)) from None

        def done(result) -> None:
            fps, picks = result
            if remember is not None:
                remember(fps)
            self.finish(everything, picks, name, len(seeds))

        return ctx.run("Picking diverse compounds…", work, done, keep_partial=True)

    def finish(self, everything: np.ndarray, picks: Picks, name: str, seeds: int) -> None:
        ctx = self.ctx
        rank = np.full(len(everything), np.nan)
        rank[picks.rows] = np.arange(1, len(picks) + 1)
        names = ctx.apply_results(
            {name: (everything, NumberColumn.from_floats(rank, 0))},
            label="Diverse subset",
        )
        picked = everything[picks.rows]
        if self.select.isChecked():
            ctx.selection.set(picked.tolist(), "diverse")
        message = f"Picked {plural(len(picks), 'compound')} in “{names[name]}”"
        if len(picks) < self.count.value() and not self.job_cancelled():
            message += " (the rest duplicate earlier picks)"
        if seeds:
            message += f", adding to {seeds:,} selected"
        ctx.show_message(message)

    def job_cancelled(self) -> bool:
        """Whether these picks are what a cancelled job had made so far."""
        return self.job is not None and self.job.cancel_token.cancelled


def _mask(n: int, rows: np.ndarray) -> np.ndarray:
    out = np.zeros(n, dtype=np.bool_)
    out[rows] = True
    return out


def _only_pool(fps, in_pool: np.ndarray, seed_rows: np.ndarray):
    """Rows outside the pool (seeds only) can anchor distances but never be picked."""
    from chemsuite.analysis.fingerprints import Fingerprints

    pickable = fps.valid & (in_pool | _mask(len(fps), seed_rows))
    return Fingerprints(fps.kind, fps.packed, pickable)
