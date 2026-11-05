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
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
)

from chemsuite.analysis.diversity import Picks, maxmin_fingerprints, maxmin_numeric
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn

from .common import FeatureChooser, fingerprint_task, scope_combo, scope_oids

RANK_COLUMN = "Diversity Rank"


class DiverseSubsetDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Diverse Subset")
        self.ctx = ctx
        self.job = None
        self.features = FeatureChooser(ctx)
        self.rows = scope_combo(ctx, "shown")
        self.count = QSpinBox()
        self.count.setRange(1, 1_000_000)
        self.count.setValue(min(50, max(1, len(ctx.table))))
        self.seeded = QCheckBox(
            f"Add to the selected rows ({len(ctx.selection):,}): pick what they do not cover"
        )
        self.seeded.setEnabled(len(ctx.selection) > 0)
        self.output = QLineEdit(RANK_COLUMN)
        self.select = QCheckBox("Select the picked rows")
        self.select.setChecked(True)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Pick")
        form = QFormLayout()
        form.addRow("Pick from:", self.rows)
        form.addRow("How many:", self.count)
        form.addRow("Result column:", self.output)
        layout = QVBoxLayout(self)
        layout.addWidget(self.features)
        layout.addLayout(form)
        layout.addWidget(self.seeded)
        layout.addWidget(self.select)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(520, 0)
        self.features.changed.connect(self._update)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    def _update(self) -> None:
        problem = self.features.problem()
        self.status.setText(
            problem
            or "Each pick is the compound least like everything picked so far (MaxMin). "
            "Rank 1 is picked first."
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if self.features.problem():
            return
        self.job = self.pick()
        super().accept()

    def pick(self):
        ctx = self.ctx
        pool = scope_oids(ctx, self.rows.currentData())
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
                return fps, maxmin_fingerprints(
                    _only_pool(fps, in_pool, seed_rows),
                    count,
                    seeds=seed_rows,
                    metric=metric,
                    progress=job.report,
                    cancel=job.cancel_token,
                )
        else:
            values = self.features.number_values(everything)
            remember = None

            def work(job):
                masked = [
                    np.where(in_pool | _mask(len(everything), seed_rows), v, np.nan) for v in values
                ]
                return None, maxmin_numeric(
                    masked, count, seeds=seed_rows, progress=job.report, cancel=job.cancel_token
                )

        def done(result) -> None:
            fps, picks = result
            if remember is not None:
                remember(fps)
            self.finish(everything, picks, name, len(seeds))

        return ctx.run(f"Picking {count:,} diverse compounds…", work, done)

    def finish(self, everything: np.ndarray, picks: Picks, name: str, seeds: int) -> None:
        ctx = self.ctx
        rank = np.full(len(everything), np.nan)
        rank[picks.rows] = np.arange(1, len(picks) + 1)
        names = ctx.apply_results(
            {name: (everything, NumberColumn.from_floats(rank, 0))},
            label="Diverse subset",
            overwrite=True,
        )
        picked = everything[picks.rows]
        if self.select.isChecked():
            ctx.selection.set(picked.tolist(), "diverse")
        message = f"Picked {len(picks):,} compound(s) in “{names[name]}”"
        if len(picks) < self.count.value():
            message += " (the rest duplicate earlier picks)"
        if seeds:
            message += f", adding to {seeds:,} selected"
        ctx.show_message(message, 8000)


def _mask(n: int, rows: np.ndarray) -> np.ndarray:
    out = np.zeros(n, dtype=np.bool_)
    out[rows] = True
    return out


def _only_pool(fps, in_pool: np.ndarray, seed_rows: np.ndarray):
    """Rows outside the pool (seeds only) can anchor distances but never be picked."""
    from chemsuite.analysis.fingerprints import Fingerprints

    pickable = fps.valid & (in_pool | _mask(len(fps), seed_rows))
    return Fingerprints(fps.kind, fps.packed, pickable)
