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
"""Data ▸ Dimensionality Reduction ▸ Generative Topographic Map…: train a map
(``analysis.gtm``), write each compound's place and log-likelihood, keep the map in the
document as an asset, and open the GTM Landscape plot on it."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QDoubleSpinBox, QLabel, QSpinBox

from chemsuite.analysis.gtm import (
    ASSET_KIND,
    ITERATIONS,
    MAX_SIDE,
    RBF_WIDTH,
    REGULARIZATION,
    GtmMap,
    default_basis,
    default_side,
    gtm,
)
from chemsuite.analysis.selectivity import looks_like_assay
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import FitChoice, RowsChoice, ToolDialog, inline
from chemsuite.core.assets import Asset, PutAsset
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural
from chemsuite.plotting.data import number_columns
from chemsuite.tools.similarity.common import FeatureChooser

from .common import features_task

X_COLUMN = "GTM X"
Y_COLUMN = "GTM Y"
LOGLIK_COLUMN = "GTM Log-likelihood"


class GtmDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(
            ctx.parent, "Generative Topographic Map", guide="generative-topographic-map"
        )
        self.ctx = ctx
        self.plot = None
        self.map: GtmMap | None = None
        inputs = self.section("Input")
        method = self.section("Method")
        more = self.section("More Options", collapsed=True)
        self.features = FeatureChooser(ctx, metric=False, section=inputs, more=more)
        self.rows = inputs.row("Rows", RowsChoice(ctx))
        self.fit = inputs.row("Fit on", FitChoice(ctx),
                              help="The rows the map is trained on; every row is then placed "
                              "on it.")  # fmt: skip

        self.side = QSpinBox()
        self.side.setRange(2, MAX_SIDE)
        self.auto = QCheckBox("Automatic")
        self.auto.setChecked(True)
        method.row("Nodes", inline(self.side, QLabel("per side"), self.auto),
                   help="The map is a square grid of nodes. Automatic: 10 per side below 300 "
                   "compounds, 16 below 3,000, else 25.")  # fmt: skip
        self.basis = QSpinBox()
        self.basis.setRange(2, MAX_SIDE)
        method.row("Basis functions", inline(self.basis, QLabel("per side")),
                   help="How freely the map can bend: fewer gives a smoother map. Set from the "
                   "nodes (a quarter as many, at least 3) while Nodes is Automatic.")  # fmt: skip
        self.iterations = QSpinBox()
        self.iterations.setRange(2, 1000)
        self.iterations.setValue(ITERATIONS)
        method.row("Iterations", self.iterations,
                   help="The most training passes; training stops sooner once the fit no "
                   "longer improves.")  # fmt: skip

        self.width = QDoubleSpinBox()
        self.width.setRange(0.1, 10.0)
        self.width.setSingleStep(0.1)
        self.width.setValue(RBF_WIDTH)
        more.row("Basis width", self.width,
                 help="Each basis function's width, in spacings between them: wider gives a "
                 "smoother map.")  # fmt: skip
        self.regularization = QDoubleSpinBox()
        self.regularization.setDecimals(3)
        self.regularization.setRange(0.0, 100.0)
        self.regularization.setSingleStep(0.05)
        self.regularization.setValue(REGULARIZATION)
        more.row("Regularization", self.regularization,
                 help="How strongly the map is kept from folding into the data: higher gives "
                 "a smoother, more general map.")  # fmt: skip

        for signal in (self.auto.toggled, self.rows.changed, self.fit.currentIndexChanged):
            signal.connect(lambda *_: self._update())
        self.features.changed.connect(self._update)
        self._update()
        self.remember(ctx, "gtm", {
            "auto": self.auto, "side": self.side, "basis": self.basis,
            "iterations": self.iterations, "width": self.width,
            "regularization": self.regularization, "fit": self.fit,
            **self.features.remembered(),
        })  # fmt: skip

    def _update(self) -> None:
        auto = self.auto.isChecked()
        if auto:
            oids = self.rows.oids()
            fit = self.fit.rows(oids)
            side = default_side(len(oids) if fit is None else len(fit))
            self.side.setValue(side)
            self.basis.setValue(default_basis(side))
        for spin in (self.side, self.basis):
            spin.setEnabled(not auto)
        self.refresh()

    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        fit = self.fit.rows(self.rows.oids())
        if fit is not None and len(fit) < 2:
            return "Train on at least 2 of the rows mapped."
        return ""

    def verb(self) -> str:
        return f"Train on {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.train()

    def recipe(self) -> dict:
        """How the features were made, kept with the map for reference."""
        chooser = self.features
        columns = chooser.number_columns() if chooser.uses_columns() else []
        return {
            "structure": chooser.structure.currentText() if chooser.fingerprint_key() else "",
            "fingerprint": chooser.fingerprint_key() or "",
            "columns": list(columns),
            "scaling": chooser.column_scaling(),
            "weight": chooser.fingerprint_weight(),
        }

    def train(self):
        ctx = self.ctx
        oids = self.rows.oids()
        build, remember = features_task(ctx, self.features, oids)
        options = {
            "side": self.side.value(),
            "per_side": self.basis.value(),
            "iterations": self.iterations.value(),
            "width": self.width.value(),
            "regularization": self.regularization.value(),
            "fit": self.fit.rows(oids),
        }
        recipe = self.recipe()

        def work(job) -> GtmMap:
            features = build(job)
            return gtm(features, progress=job.report, cancel=job.cancel_token, **options)

        def done(result: GtmMap) -> None:
            remember()
            self.finish(result.with_oids(oids[result.oids], recipe))

        return ctx.run("Calculating Generative Topographic Map…", work, done)

    def finish(self, gmap: GtmMap):
        ctx = self.ctx
        self.map = gmap
        asset_id = ctx.table.unique_asset_id("gtm")
        title = f"GTM {gmap.side} × {gmap.side}, {plural(len(gmap), 'compound')}"
        asset = Asset(asset_id, ASSET_KIND, title, gmap, {"side": str(gmap.side)})
        oids = gmap.oids
        value = next((n for n in number_columns(ctx.table) if looks_like_assay(n)), None)
        ctx.apply_results(
            {
                X_COLUMN: (oids, NumberColumn.from_floats(gmap.position[:, 0], 3)),
                Y_COLUMN: (oids, NumberColumn.from_floats(gmap.position[:, 1], 3)),
                LOGLIK_COLUMN: (oids, NumberColumn.from_floats(gmap.loglik, 3)),
            },
            label="Generative topographic map",
            extra_ops=[PutAsset(asset_id, None, asset)],
        )
        self.plot = ctx.open_tool("gtm")
        self.plot.set_spec(gtm_map=asset_id, gtm_value=value)
        low = float(np.percentile(gmap.loglik, 5)) if len(gmap) else float("nan")
        ctx.show_message(
            f"Mapped {plural(len(gmap), 'compound')} on a {gmap.side} × {gmap.side} map in "
            f"{plural(gmap.iterations, 'iteration')}; 5% have a log-likelihood below {low:.3g}"
        )
        return self.plot
