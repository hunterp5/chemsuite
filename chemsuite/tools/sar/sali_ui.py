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
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import numbers
from chemsuite.sar.sali import Cliffs, sali_pairs
from chemsuite.tools.similarity.common import (
    FeatureChooser,
    fingerprint_task,
    labels,
    scope_combo,
    scope_oids,
)

from .common import ResultTable, activity_combo

SALI_COLUMN = "Max SALI"
PARTNER_COLUMN = "Cliff Partner"
SHOWN_PAIRS = 2000


class SaliDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Activity Cliffs (SALI)")
        self.ctx = ctx
        self.job = None
        self.panel: SaliPanel | None = None
        self.features = FeatureChooser(ctx, numbers_allowed=False)
        self.activity = activity_combo(ctx)
        self.rows = scope_combo(ctx, "shown")
        self.min_similarity = QDoubleSpinBox()
        self.min_similarity.setRange(0.3, 1.0)
        self.min_similarity.setSingleStep(0.05)
        self.min_similarity.setValue(0.6)
        self.min_similarity.setToolTip("Only pairs at least this similar are compared.")
        self.skip_identical = QCheckBox("Skip identical fingerprints (salts, stereoisomers)")
        self.skip_identical.setChecked(True)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Find Cliffs")
        form = QFormLayout()
        form.addRow("Activity:", self.activity)
        form.addRow("Rows:", self.rows)
        form.addRow("Minimum similarity:", self.min_similarity)
        form.addRow("", self.skip_identical)
        layout = QVBoxLayout(self)
        layout.addWidget(self.features)
        layout.addLayout(form)
        layout.addWidget(
            QLabel("SALI = |Δ activity| / (1 − similarity). Writes Max SALI and Cliff Partner.")
        )
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(500, 0)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        problem = self.problem()
        self.status.setText(problem)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def problem(self) -> str:
        if not self.activity.currentText():
            return "There is no number column for activity."
        return self.features.problem()

    def accept(self) -> None:
        if self.problem():
            return
        self.job = self.find()
        super().accept()

    def find(self):
        ctx = self.ctx
        oids = scope_oids(ctx, self.rows.currentData())
        name = self.activity.currentText()
        activity = numbers(ctx.table, name, oids)
        metric = self.features.metric.currentData()
        threshold = self.min_similarity.value()
        skip_identical = self.skip_identical.isChecked()
        compute, remember = fingerprint_task(
            ctx, self.features.structure.currentText(), self.features.kind.currentData(), oids
        )

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

        def done(result) -> None:
            fps, cliffs = result
            remember(fps)
            self.finish(oids, activity, name, cliffs)

        return ctx.run("Finding activity cliffs…", work, done)

    def finish(self, oids, activity, name: str, cliffs: Cliffs):
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
            overwrite=True,
        )
        self.panel = SaliPanel(ctx, oids, activity, name, cliffs, names)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"{len(cliffs):,} similar pair(s) scored; highest SALI per compound in "
            f"“{written[SALI_COLUMN]}”",
            10_000,
        )
        return self.panel


class SaliPanel(QWidget):
    titleChanged = Signal(str)
    tool_id = "sali"

    def __init__(self, ctx: ToolContext, oids, activity, name, cliffs: Cliffs, names) -> None:
        super().__init__()
        self.ctx = ctx
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
        return {}

    def set_state(self, state: dict) -> None:
        pass

    def close_panel(self) -> None:
        pass
