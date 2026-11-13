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
"""Activity Cliffs (SALI)…: similar compounds with very different activity."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
)

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural
from chemsuite.plotting.data import numbers
from chemsuite.sar.sali import Cliffs, sali_pairs
from chemsuite.tools.similarity.common import (
    FeatureChooser,
    fingerprint_task,
    labels,
)

from .common import RerunPanel, ResultTable, activity_combo, saved_rows

SALI_COLUMN = "Max SALI"
PARTNER_COLUMN = "Cliff Partner"
SHOWN_PAIRS = 2000


class SaliDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Activity Cliffs (SALI)", guide="activity-cliffs-sali")
        self.ctx = ctx
        self.panel: SaliPanel | None = None
        inputs = self.section("Input")
        self.features = FeatureChooser(ctx, numbers_allowed=False, section=inputs)
        self.activity = inputs.row("Activity", activity_combo(ctx))
        self.rows = inputs.row("Rows", RowsChoice(ctx, "shown"))

        method = self.section("Method")
        self.min_similarity = QDoubleSpinBox()
        self.min_similarity.setRange(0.3, 1.0)
        self.min_similarity.setSingleStep(0.05)
        self.min_similarity.setValue(0.6)
        method.row(
            "Min. similarity",
            self.min_similarity,
            help="Only pairs at least this similar are compared. "
            "SALI = |Δ activity| / (1 − similarity).",
        )
        self.skip_identical = QCheckBox("Skip identical fingerprints")
        self.skip_identical.setChecked(True)
        method.row(
            "Identical",
            self.skip_identical,
            help="Salts and stereoisomers share a fingerprint and would divide by zero.",
        )

        for signal in (self.features.changed, self.activity.currentIndexChanged, self.rows.changed):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.activity.currentText():
            return "There is no number column for activity."
        return self.features.problem()

    def verb(self) -> str:
        return f"Find Cliffs in {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.find()

    def recipe(self) -> dict:
        """Everything the analysis needs, saved with its panel to rerun it on reopening."""
        return {
            "oids": self.rows.oids().tolist(),
            "activity": self.activity.currentText(),
            "structure": self.features.structure.currentText(),
            "fingerprint": self.features.kind.currentData(),
            "metric": self.features.metric.currentData(),
            "min_similarity": self.min_similarity.value(),
            "skip_identical": self.skip_identical.isChecked(),
        }

    def find(self):
        recipe = self.recipe()
        return _compute(
            self.ctx,
            recipe,
            lambda oids, activity, cliffs: self.finish(
                oids, activity, recipe["activity"], cliffs, recipe
            ),
        )

    def finish(self, oids, activity, name: str, cliffs: Cliffs, recipe: dict | None = None):
        ctx = self.ctx
        best, partner = cliffs.per_row(len(oids))
        names = labels(ctx, oids)
        partners = [names[p] if p >= 0 else "" for p in partner.tolist()]
        written = ctx.apply_results(
            {
                SALI_COLUMN: (oids, NumberColumn.from_floats(best, 2)),
                PARTNER_COLUMN: (oids, partners),
            },
            label="Activity cliffs",
        )
        self.panel = SaliPanel(ctx, oids, activity, name, cliffs, names, recipe)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"Scored {plural(len(cliffs), 'similar pair')}; highest SALI per compound in "
            f"“{written[SALI_COLUMN]}”"
        )
        return self.panel


class SaliPanel(ToolPanel):
    tool_id = "sali"

    def __init__(
        self, ctx: ToolContext, oids, activity, name, cliffs: Cliffs, names, recipe=None
    ) -> None:
        super().__init__()
        self.ctx = ctx
        self.recipe = recipe or {}
        self.oids = np.asarray(oids)
        self.name = name
        shown = min(len(cliffs), SHOWN_PAIRS)
        pairs = list(zip(cliffs.a[:shown].tolist(), cliffs.b[:shown].tolist(), strict=True))
        self.table = ResultTable(
            ["Compound A", "Compound B", "Similarity", f"{name} A", f"{name} B", "SALI"]
        )
        self.table.fill(
            [
                [names[a], names[b], cliffs.similarity[i], activity[a], activity[b], cliffs.sali[i]]
                for i, (a, b) in enumerate(pairs)
            ],
            pairs,
        )
        self.table.itemDoubleClicked.connect(lambda *_: self.select_pairs())
        which = (
            f"The {shown:,} steepest of {len(cliffs):,}"
            if shown < len(cliffs)
            else f"All {len(cliffs):,}"
        )
        note = QLabel(f"{which} similar pair(s). Double-click a pair to select both compounds.")
        note.setWordWrap(True)
        select = QPushButton("Select Compounds")
        select.clicked.connect(lambda: self.select_pairs())
        copy = QPushButton("Copy")
        copy.clicked.connect(lambda: self.table.copy_text())
        buttons = QHBoxLayout()
        buttons.addWidget(select)
        buttons.addStretch(1)
        buttons.addWidget(copy)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(note)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)

    def title(self) -> str:
        return f"Activity Cliffs: {self.name}"

    def select_pairs(self) -> list[int]:
        rows = sorted({r for pair in self.table.selected_data() if pair for r in pair})
        oids = self.oids[rows].tolist()
        if oids:
            self.ctx.selection.set(oids, "sali")
        return oids

    def state(self) -> dict:
        return dict(self.recipe)

    def close_panel(self) -> None:
        pass


def _compute(ctx: ToolContext, recipe: dict, done):
    """Score the pairs *recipe* describes in the background; then
    ``done(oids, activity, cliffs)``."""
    oids = np.asarray(recipe["oids"], dtype=np.int64)
    activity = numbers(ctx.table, recipe["activity"], oids)
    threshold, metric = float(recipe["min_similarity"]), recipe["metric"]
    skip_identical = bool(recipe["skip_identical"])
    compute, remember = fingerprint_task(ctx, recipe["structure"], recipe["fingerprint"], oids)

    def work(job):
        fps = compute(job)
        cliffs = sali_pairs(
            fps,
            activity,
            min_similarity=threshold,
            metric=metric,
            skip_identical=skip_identical,
            progress=job.report,
            cancel=job.cancel_token,
        )
        return fps, cliffs

    def finished(result) -> None:
        fps, cliffs = result
        remember(fps)
        done(oids, activity, cliffs)

    return ctx.run("Finding activity cliffs…", work, finished)


def restore(ctx: ToolContext, state: dict) -> RerunPanel:
    """Reopen a saved panel: rerun the analysis on the rows it covered (columns stay as they
    were saved)."""
    recipe = {**state, "oids": saved_rows(ctx, state, state["activity"], state["structure"])}
    holder = RerunPanel(SaliPanel.tool_id, f"Activity Cliffs: {recipe['activity']}", recipe)
    ctx.add_panel(holder, holder.title())

    def done(oids, activity, cliffs) -> None:
        names = labels(ctx, oids)
        holder.fill(SaliPanel(ctx, oids, activity, recipe["activity"], cliffs, names, recipe))

    holder.job = _compute(ctx, recipe, done)
    return holder
