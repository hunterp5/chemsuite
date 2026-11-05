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
"""Chemical Space Map…: PCA, UMAP or t-SNE coordinates as columns, shown as a scatter plot."""

from __future__ import annotations

import importlib.util

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
)

from chemsuite.analysis import embedding as em
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import colorable_columns
from chemsuite.tools.sar.common import activity_combo
from chemsuite.tools.similarity.common import FeatureChooser, scope_combo, scope_oids

from .common import features_task

NONE = "(none)"
HELP = {
    "pca": "Fast and exact; distances along the axes are meaningful, but chemistry often "
    "needs more than two components, so groups can overlap.",
    "umap": "Keeps neighbourhoods together, so series form islands; distances between islands "
    "mean little. Above 50,000 rows it is fitted on a sample.",
    "tsne": "Like UMAP, good at showing local groups; slower, and only for the rows given.",
}


class SpaceMapDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Chemical Space Map")
        self.ctx = ctx
        self.job = None
        self.plot = None
        self.features = FeatureChooser(ctx, metric=False)
        self.method = QComboBox()
        has_umap = importlib.util.find_spec("umap") is not None
        for key, title in em.METHODS.items():
            self.method.addItem(title, key)
            if key == "umap" and not has_umap:
                item = self.method.model().item(self.method.count() - 1)
                item.setEnabled(False)
                item.setToolTip("UMAP needs the umap-learn package.")
        self.method.setCurrentIndex(self.method.findData("umap" if has_umap else "pca"))
        self.about = QLabel()
        self.about.setWordWrap(True)
        self.neighbours = QSpinBox()
        self.neighbours.setRange(2, 200)
        self.neighbours.setValue(15)
        self.neighbours.setToolTip("Small: tight local groups. Large: more of the global layout.")
        self.min_dist = QDoubleSpinBox()
        self.min_dist.setRange(0.0, 0.99)
        self.min_dist.setSingleStep(0.05)
        self.min_dist.setValue(0.1)
        self.min_dist.setToolTip("How tightly points may pack together.")
        self.perplexity = QDoubleSpinBox()
        self.perplexity.setRange(2, 200)
        self.perplexity.setValue(30)
        self.perplexity.setToolTip("Roughly the number of neighbours each point attends to.")
        self.repeatable = QCheckBox("Same map every run (slower: one processor core)")
        self.repeatable.setChecked(True)
        self.rows = scope_combo(ctx)
        self.color = QComboBox()
        self.color.addItems([NONE, *colorable_columns(ctx.table)])
        guess = activity_combo(ctx).currentText()
        if guess:
            self.color.setCurrentText(guess)
        self.output = QLineEdit()
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Make Map")
        self.form = QFormLayout()
        self.form.addRow("Method:", self.method)
        self.form.addRow("Neighbours:", self.neighbours)
        self.form.addRow("Minimum distance:", self.min_dist)
        self.form.addRow("Perplexity:", self.perplexity)
        self.form.addRow("", self.repeatable)
        self.form.addRow("Rows:", self.rows)
        self.form.addRow("Colour by:", self.color)
        self.form.addRow("Column names:", self.output)
        layout = QVBoxLayout(self)
        layout.addWidget(self.features)
        layout.addLayout(self.form)
        layout.addWidget(self.about)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(540, 0)
        self._auto_name = ""
        self.method.currentIndexChanged.connect(lambda *_: self._update())
        self.rows.currentIndexChanged.connect(lambda *_: self._update())
        self.features.changed.connect(self._update)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    def _count(self) -> int:
        return len(scope_oids(self.ctx, self.rows.currentData()))

    def problem(self) -> str:
        return self.features.problem()

    def _update(self) -> None:
        method = self.method.currentData()
        self.about.setText(HELP[method])
        self.form.setRowVisible(self.neighbours, method == "umap")
        self.form.setRowVisible(self.min_dist, method == "umap")
        self.form.setRowVisible(self.perplexity, method == "tsne")
        self.form.setRowVisible(self.repeatable, method != "pca")
        name = em.METHODS[method]
        if self.output.text() in ("", self._auto_name):
            self.output.setText(name)
        self._auto_name = name
        problem = self.problem()
        note = problem
        rows = self._count()
        if not problem and method == "tsne" and rows > em.TSNE_WARN_ROWS:
            note = f"t-SNE on {rows:,} rows may take many minutes; UMAP is quicker at this size."
        elif not problem and method == "umap":
            note = "The first UMAP run in a session takes longer while it compiles."
        self.status.setText(note)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if self.problem():
            return
        self.job = self.make_map()
        super().accept()

    def options(self) -> dict:
        return {
            "neighbours": self.neighbours.value(),
            "min_dist": self.min_dist.value(),
            "perplexity": self.perplexity.value(),
            "seed": 42 if self.repeatable.isChecked() else None,
        }

    def make_map(self):
        ctx = self.ctx
        oids = scope_oids(ctx, self.rows.currentData())
        method = self.method.currentData()
        options = self.options()
        build, remember = features_task(ctx, self.features, oids)
        prefix = self.output.text().strip() or em.METHODS[method]
        color = self.color.currentText()

        def work(job):
            features = build(job)
            return em.embed(
                features, method, progress=job.report, cancel=job.cancel_token, **options
            )

        def done(result: em.Embedding) -> None:
            remember()
            self.finish(oids, result, prefix, None if color == NONE else color)

        return ctx.run(f"Calculating {em.METHODS[method]}…", work, done)

    def finish(self, oids: np.ndarray, result: em.Embedding, prefix: str, color: str | None):
        ctx = self.ctx
        names = ctx.apply_results(
            {
                f"{prefix} 1": (oids, NumberColumn.from_floats(result.coords[:, 0], 3)),
                f"{prefix} 2": (oids, NumberColumn.from_floats(result.coords[:, 1], 3)),
            },
            label=f"{em.METHODS[result.method]} map",
            overwrite=True,
        )
        self.plot = ctx.open_tool("scatter")
        self.plot.set_spec(x=names[f"{prefix} 1"], y=names[f"{prefix} 2"], color=color)
        message = f"{result.placed:,} of {len(oids):,} row(s) placed"
        if result.notes:
            message += "; " + "; ".join(result.notes)
        ctx.show_message(message, 10_000)
        return self.plot
