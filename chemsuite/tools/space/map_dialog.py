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
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox

from chemsuite.analysis import embedding as em
from chemsuite.app.dialogs import OutputChoice, RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import colorable_columns
from chemsuite.tools.sar.common import activity_combo
from chemsuite.tools.similarity.common import FeatureChooser

from .common import features_task

NONE = "(none)"
HELP = {
    "pca": "Fast and exact; distances along the axes mean something, but two components "
    "often cannot separate every group.",
    "umap": "Keeps neighbourhoods together, so series form islands; distances between "
    "islands mean little. Above 50,000 rows it is fitted on a sample.",
    "tsne": "Like UMAP, good at local groups; slower, and only for the rows given.",
}


class SpaceMapDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(
            ctx.parent,
            "Chemical Space Map",
            guide="chemical-space-map",
        )
        self.ctx = ctx
        self.plot = None

        inputs = self.section("Input")
        self.features = FeatureChooser(ctx, metric=False, section=inputs)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.method = QComboBox()
        has_umap = importlib.util.find_spec("umap") is not None
        for key, title in em.METHODS.items():
            self.method.addItem(title, key)
            if key == "umap" and not has_umap:
                item = self.method.model().item(self.method.count() - 1)
                item.setEnabled(False)
                item.setToolTip("UMAP needs the umap-learn package.")
        self.method.setCurrentIndex(self.method.findData("umap" if has_umap else "pca"))
        method.row("Method", self.method)
        self.neighbours = QSpinBox()
        self.neighbours.setRange(2, 200)
        self.neighbours.setValue(15)
        method.row(
            "Neighbours",
            self.neighbours,
            help="Fewer keeps local detail; more shows the overall layout.",
        )
        self.min_dist = QDoubleSpinBox()
        self.min_dist.setRange(0.0, 0.99)
        self.min_dist.setSingleStep(0.05)
        self.min_dist.setValue(0.1)
        method.row("Min. distance", self.min_dist, help="How tightly points may pack together.")
        self.perplexity = QDoubleSpinBox()
        self.perplexity.setRange(2, 200)
        self.perplexity.setValue(30)
        method.row(
            "Perplexity",
            self.perplexity,
            help="Roughly how many neighbours each point attends to.",
        )
        self.repeatable = QCheckBox("Same map every run (slower: one processor core)")
        self.repeatable.setChecked(True)
        method.row("", self.repeatable)
        self._method_section = method

        output = self.section("Output")
        self.output = QLineEdit()
        output.row("Column names", self.output)
        self.if_exists = output.row("If they exist", OutputChoice())
        self.color = QComboBox()
        self.color.addItems([NONE, *colorable_columns(ctx.table)])
        guess = activity_combo(ctx).currentText()
        if guess:
            self.color.setCurrentText(guess)
        output.row("Colour plot by", self.color)
        self._output_section = output

        self._auto_name = ""
        for signal in (
            self.method.currentIndexChanged,
            self.rows.changed,
            self.features.changed,
            self.output.textChanged,
            self.color.currentIndexChanged,
        ):
            signal.connect(lambda *_: self._update())
        self._update()

    def _count(self) -> int:
        return len(self.rows.oids())

    def _prefix(self) -> str:
        return self.output.text().strip() or em.METHODS[self.method.currentData()]

    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        if self._count() < 3:
            return "A map needs at least 3 rows."
        return ""

    def note(self) -> str:
        method = self.method.currentData()
        rows = self._count()
        if method == "tsne" and rows > em.TSNE_WARN_ROWS:
            return f"t-SNE on {rows:,} rows may take many minutes; UMAP is quicker at this size."
        if method == "umap":
            return "The first UMAP run in a session takes longer while it compiles."
        return ""

    def summary(self) -> str:
        method = em.METHODS[self.method.currentData()]
        prefix = self._prefix()
        color = self.color.currentText()
        plot = "a scatter plot" + (f" coloured by {color}" if color != NONE else "")
        return f"{method} of {self.features.describe()} → {prefix} 1, {prefix} 2, and {plot}"

    def verb(self) -> str:
        return f"Map {self._count():,} Rows"

    def _update(self) -> None:
        method = self.method.currentData()
        section = self._method_section
        section.set_help(self.method, HELP[method])
        section.set_row_visible(self.neighbours, method == "umap")
        section.set_row_visible(self.min_dist, method == "umap")
        section.set_row_visible(self.perplexity, method == "tsne")
        section.set_row_visible(self.repeatable, method != "pca")
        name = em.METHODS[method]
        if self.output.text() in ("", self._auto_name):
            self.output.blockSignals(True)
            self.output.setText(name)
            self.output.blockSignals(False)
        self._auto_name = name
        self.refresh()

    def run(self):
        return self.make_map()

    def options(self) -> dict:
        return {
            "neighbours": self.neighbours.value(),
            "min_dist": self.min_dist.value(),
            "perplexity": self.perplexity.value(),
            "seed": 42 if self.repeatable.isChecked() else None,
        }

    def make_map(self):
        ctx = self.ctx
        oids = self.rows.oids()
        method = self.method.currentData()
        overwrite = self.if_exists.overwrite()
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
            self.finish(oids, result, prefix, None if color == NONE else color, overwrite)

        return ctx.run(f"Calculating {em.METHODS[method]}…", work, done)

    def finish(
        self,
        oids: np.ndarray,
        result: em.Embedding,
        prefix: str,
        color: str | None,
        overwrite: bool = True,
    ):
        ctx = self.ctx
        names = ctx.apply_results(
            {
                f"{prefix} 1": (oids, NumberColumn.from_floats(result.coords[:, 0], 3)),
                f"{prefix} 2": (oids, NumberColumn.from_floats(result.coords[:, 1], 3)),
            },
            label=f"{em.METHODS[result.method]} map",
            overwrite=overwrite,
        )
        self.plot = ctx.open_tool("scatter")
        self.plot.set_spec(x=names[f"{prefix} 1"], y=names[f"{prefix} 2"], color=color)
        message = f"{result.placed:,} of {len(oids):,} row(s) placed"
        if result.notes:
            message += "; " + "; ".join(result.notes)
        ctx.show_message(message, 10_000)
        return self.plot
