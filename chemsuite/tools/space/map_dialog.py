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
from PySide6.QtWidgets import QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox

from chemsuite.analysis import embedding as em
from chemsuite.app.dialogs import OutputChoice, RowsChoice, Segmented, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import colorable_columns
from chemsuite.tools.sar.common import activity_combo
from chemsuite.tools.similarity.common import FeatureChooser

from .common import features_task

NONE = "(none)"
#: the PCA axes offered: (first, second) component, 1-based
COMPONENT_PAIRS = ((1, 2), (1, 3), (2, 3), (3, 4))
HELP = {
    "pca": "Fast and exact; distances along the axes mean something, but two components "
    "often cannot separate every group.",
    "umap": "Keeps neighbourhoods together, so series form islands; distances between "
    "islands mean little. Above 50,000 rows it is fitted on a sample.",
    "tsne": "Like UMAP, good at local groups; slower, and only for the rows given.",
    "pacmap": "Keeps local groups and the overall layout together, so distances between "
    "islands mean more than in UMAP or t-SNE; fast on large sets.",
}
#: default neighbours per method (the box is shared, so switching method updates it)
NEIGHBOURS = {"umap": 15, "pacmap": 10}


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
        missing = {
            "umap": None if has_umap else "UMAP needs the umap-learn package.",
            "pacmap": None
            if importlib.util.find_spec("pacmap") is not None
            else "PaCMAP needs the pacmap package.",
        }
        for key, title in em.METHODS.items():
            self.method.addItem(title, key)
            if missing.get(key):
                item = self.method.model().item(self.method.count() - 1)
                item.setEnabled(False)
                item.setToolTip(missing[key])
        self.method.setCurrentIndex(self.method.findData("umap" if has_umap else "pca"))
        method.row("Method", self.method)
        self.dims = Segmented((("2", "2D"), ("3", "3D")), "2")
        method.row(
            "Dimensions",
            self.dims,
            help="3D writes three columns and opens a 3D Scatter you can turn.",
        )
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
        self.components = QComboBox()
        for first, second in COMPONENT_PAIRS:
            self.components.addItem(f"PC{first} and PC{second}")
        method.row(
            "Components",
            self.components,
            help="Which principal components make the two axes. Later components show "
            "differences the first two miss.",
        )
        self._method_section = method
        self._method_key = self.method.currentData()

        more = self.section("More Options", collapsed=True)
        self.distance = QComboBox()
        self.distance.addItem("Automatic", None)
        for key, title in em.UMAP_METRICS.items():
            self.distance.addItem(title, key)
        more.row(
            "Distance",
            self.distance,
            help="How unlike two compounds are. Automatic: Jaccard (Tanimoto) on fingerprint "
            "bits, Euclidean on z-scored columns. Jaccard, Dice and Rogers-Tanimoto need "
            "fingerprints.",
        )
        self.spread = QDoubleSpinBox()
        self.spread.setRange(0.1, 10.0)
        self.spread.setSingleStep(0.1)
        self.spread.setValue(1.0)
        more.row(
            "Spread",
            self.spread,
            help="The overall scale of the map; with Min. distance, how clumped the islands are.",
        )
        self.epochs = QSpinBox()
        self.epochs.setRange(0, 5000)
        self.epochs.setSingleStep(50)
        self.epochs.setSpecialValueText("Automatic")
        more.row(
            "Epochs",
            self.epochs,
            help="Optimisation rounds. Automatic: 500 for small sets, 200 for large ones.",
        )
        self.fit_max = QSpinBox()
        self.fit_max.setRange(1_000, 1_000_000)
        self.fit_max.setSingleStep(10_000)
        self.fit_max.setGroupSeparatorShown(True)
        self.fit_max.setValue(em.UMAP_FIT_MAX)
        more.row(
            "Fit on up to",
            self.fit_max,
            help="Larger sets are fitted on a random sample of this many rows and the rest "
            "placed onto that map: faster, and nearly the same picture.",
        )
        self.learning_rate = QDoubleSpinBox()
        self.learning_rate.setRange(0, 5000)
        self.learning_rate.setSingleStep(50)
        self.learning_rate.setSpecialValueText("Automatic")
        more.row(
            "Learning rate",
            self.learning_rate,
            help="Step size. Automatic scales it with the number of rows; raise it if the map "
            "looks like a compressed ball.",
        )
        self.iterations = QSpinBox()
        self.iterations.setRange(250, 10_000)
        self.iterations.setSingleStep(250)
        self.iterations.setValue(1000)
        more.row("Iterations", self.iterations, help="More can settle a large map further.")
        self.exaggeration = QDoubleSpinBox()
        self.exaggeration.setRange(1.0, 100.0)
        self.exaggeration.setValue(12.0)
        more.row(
            "Exaggeration",
            self.exaggeration,
            help="Early exaggeration: higher pushes the groups further apart at the start.",
        )
        self.reduce_to = QSpinBox()
        self.reduce_to.setRange(2, 500)
        self.reduce_to.setValue(em.REDUCED_DIMS)
        more.row(
            "PCA dimensions",
            self.reduce_to,
            help="t-SNE runs on this many principal components of the features.",
        )
        self.init = QComboBox()
        self.init.addItem("PCA", "pca")
        self.init.addItem("Random", "random")
        more.row(
            "Start from",
            self.init,
            help="PCA keeps the overall layout more stable; Random is the classic t-SNE.",
        )
        self.seed = QSpinBox()
        self.seed.setRange(-1, 2_147_483_647)
        self.seed.setSpecialValueText("Random")
        self.seed.setValue(42)
        more.row(
            "Seed",
            self.seed,
            help="The same seed gives the same map every run (UMAP then uses one processor "
            "core). Random is faster but differs each time.",
        )
        self.mn_ratio = QDoubleSpinBox()
        self.mn_ratio.setRange(0.0, 10.0)
        self.mn_ratio.setSingleStep(0.1)
        self.mn_ratio.setValue(0.5)
        more.row(
            "Mid-near pairs",
            self.mn_ratio,
            help="Mid-near pairs per neighbour pair; more keeps the overall layout (global "
            "structure) better.",
        )
        self.fp_ratio = QDoubleSpinBox()
        self.fp_ratio.setRange(0.0, 20.0)
        self.fp_ratio.setSingleStep(0.5)
        self.fp_ratio.setValue(2.0)
        more.row(
            "Further pairs",
            self.fp_ratio,
            help="Further pairs per neighbour pair; more pushes separate groups further apart.",
        )
        self.pacmap_distance = QComboBox()
        for key, title in em.PACMAP_METRICS.items():
            self.pacmap_distance.addItem(title, key)
        more.row(
            "PaCMAP distance",
            self.pacmap_distance,
            help="Measured on the PCA dimensions of the features. Angular is the cosine distance.",
        )
        self._more_section = more

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
            self.dims.changed,
            self.distance.currentIndexChanged,
            self.min_dist.valueChanged,
            self.spread.valueChanged,
        ):
            signal.connect(lambda *_: self._update())
        self._update()

    def _count(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        if self._count() < 3:
            return "A map needs at least 3 rows."
        if self.method.currentData() == "umap":
            metric = self.distance.currentData()
            if metric in em.BIT_METRICS and not self.features.fingerprints():
                return f"{em.UMAP_METRICS[metric]} distance needs fingerprints."
            if self.min_dist.value() > self.spread.value():
                return "Min. distance must not be larger than the spread."
        return ""

    def verb(self) -> str:
        return f"Map {self._count():,} Rows"

    def _update(self) -> None:
        method = self.method.currentData()
        section = self._method_section
        section.set_help(self.method, HELP[method])
        if method != self._method_key:  # each method's own default number of neighbours
            if self.neighbours.value() == NEIGHBOURS.get(self._method_key):
                self.neighbours.setValue(NEIGHBOURS.get(method, self.neighbours.value()))
            self._method_key = method
        three = self.dims.value() == "3"
        section.set_row_visible(self.neighbours, method in ("umap", "pacmap"))
        section.set_row_visible(self.min_dist, method == "umap")
        section.set_row_visible(self.perplexity, method == "tsne")
        section.set_row_visible(self.components, method == "pca" and not three)
        more = self._more_section
        for widget in (self.distance, self.spread, self.epochs, self.fit_max):
            more.set_row_visible(widget, method == "umap")
        for widget in (
            self.learning_rate,
            self.iterations,
            self.exaggeration,
            self.init,
        ):
            more.set_row_visible(widget, method == "tsne")
        more.set_row_visible(self.reduce_to, method in ("tsne", "pacmap"))
        for widget in (self.mn_ratio, self.fp_ratio, self.pacmap_distance):
            more.set_row_visible(widget, method == "pacmap")
        more.setVisible(method != "pca")
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
            "components": COMPONENT_PAIRS[self.components.currentIndex()],
            "neighbours": self.neighbours.value(),
            "min_dist": self.min_dist.value(),
            "spread": self.spread.value(),
            "metric": self.distance.currentData(),
            "epochs": self.epochs.value() or None,
            "fit_max": self.fit_max.value(),
            "perplexity": self.perplexity.value(),
            "learning_rate": self.learning_rate.value() or None,
            "iterations": self.iterations.value(),
            "exaggeration": self.exaggeration.value(),
            "reduce_to": self.reduce_to.value(),
            "init": self.init.currentData(),
            "seed": None if self.seed.value() < 0 else self.seed.value(),
            "mn_ratio": self.mn_ratio.value(),
            "fp_ratio": self.fp_ratio.value(),
            "distance": self.pacmap_distance.currentData(),
            "dims": int(self.dims.value()),
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
        dims = result.coords.shape[1]
        axes = [f"{prefix} {i + 1}" for i in range(dims)]
        names = ctx.apply_results(
            {
                axis: (oids, NumberColumn.from_floats(result.coords[:, i], 3))
                for i, axis in enumerate(axes)
            },
            label=f"{em.METHODS[result.method]} map",
            overwrite=overwrite,
        )
        columns = [names[a] for a in axes]
        if dims == 3:
            self.plot = ctx.open_tool("scatter3d")
            self.plot.set_spec(x=columns[0], y=columns[1], z=columns[2], color=color)
        else:
            self.plot = ctx.open_tool("scatter")
            self.plot.set_spec(x=columns[0], y=columns[1], color=color)
        message = f"{result.placed:,} of {len(oids):,} row(s) placed"
        if result.notes:
            message += "; " + "; ".join(result.notes)
        ctx.show_message(message, 10_000)
        return self.plot
