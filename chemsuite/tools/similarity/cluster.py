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
"""Cluster: Butina, k-means or hierarchical clustering by fingerprints or number columns."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox

from chemsuite.analysis import clustering as cl
from chemsuite.app.dialogs import OutputChoice, RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn

from .common import FeatureChooser, fingerprint_task

CLUSTER_COLUMN = "Cluster"
METHOD_HELP = {
    "butina": "Any number of clusters: compounds at least as similar as the cutoff group "
    "around a centre.",
    "kmeans": "A fixed number of clusters of similar size.",
    "hierarchical": "A fixed number of clusters, merged step by step; choose how distances "
    "between groups are measured.",
}


class ClusterDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(
            ctx.parent,
            "Cluster",
            guide="cluster",
        )
        self.ctx = ctx

        inputs = self.section("Input")
        self.features = FeatureChooser(ctx, section=inputs)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.method = QComboBox()
        for key, title in cl.METHODS.items():
            self.method.addItem(title, key)
        method.row("Method", self.method)
        self.cutoff = QDoubleSpinBox()
        self.cutoff.setRange(0.05, 1.0)
        self.cutoff.setSingleStep(0.05)
        self.cutoff.setValue(0.6)
        method.row(
            "Cutoff",
            self.cutoff,
            help="Compounds at least this similar are neighbours. Higher gives smaller, "
            "tighter clusters.",
        )
        self.k = QSpinBox()
        self.k.setRange(2, 10_000)
        self.k.setValue(10)
        method.row("Clusters", self.k)
        self.linkage = QComboBox()
        for key, title in cl.LINKAGES.items():
            self.linkage.addItem(title, key)
        method.row("Linkage", self.linkage)
        self._method_section = method

        output = self.section("Output")
        self.output = QLineEdit(CLUSTER_COLUMN)
        output.row("Column name", self.output)
        self.if_exists = output.row("If it exists", OutputChoice())
        self._output_section = output

        for signal in (
            self.method.currentIndexChanged,
            self.linkage.currentIndexChanged,
            self.cutoff.valueChanged,
            self.k.valueChanged,
            self.rows.changed,
            self.features.changed,
            self.output.textChanged,
        ):
            signal.connect(lambda *_: self._update())
        self._update()

    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        problem = self.features.problem()
        method = self.method.currentData()
        fingerprints = self.features.fingerprints()
        rows = self._rows()
        if problem:
            return problem
        if not rows:
            return "There are no rows to cluster."
        if method == "butina" and not fingerprints:
            return "Butina works on fingerprints; choose k-means or hierarchical for numbers."
        if method == "hierarchical" and rows > cl.MAX_HIERARCHICAL:
            return (
                f"Hierarchical clustering handles up to {cl.MAX_HIERARCHICAL:,} rows "
                f"(here {rows:,}). Use Butina or k-means, or fewer rows."
            )
        if method == "hierarchical" and fingerprints and self.linkage.currentData() == "ward":
            return "Ward linkage needs number columns: choose average or complete."
        return ""

    def verb(self) -> str:
        return f"Cluster {self._rows():,} Rows"

    def _update(self) -> None:
        method = self.method.currentData()
        section = self._method_section
        section.set_row_visible(self.cutoff, method == "butina")
        section.set_row_visible(self.k, method != "butina")
        section.set_row_visible(self.linkage, method == "hierarchical")
        section.set_help(self.method, METHOD_HELP[method])
        self.refresh()

    def run(self):
        return self.cluster()

    def cluster(self):
        ctx = self.ctx
        oids = self.rows.oids()
        method = self.method.currentData()
        overwrite = self.if_exists.overwrite()
        k, cutoff, linkage = self.k.value(), self.cutoff.value(), self.linkage.currentData()
        name = self.output.text().strip() or CLUSTER_COLUMN
        if self.features.fingerprints():
            metric = self.features.metric.currentData()
            compute, remember = fingerprint_task(
                ctx,
                self.features.structure.currentText(),
                self.features.kind.currentData(),
                oids,
            )

            def work(job):
                fps = compute(job)
                options = {"progress": job.report, "cancel": job.cancel_token}
                if method == "butina":
                    return fps, cl.butina(fps, cutoff, metric, **options)
                if method == "kmeans":
                    return fps, cl.kmeans_fingerprints(fps, k, **options)
                return fps, cl.hierarchical_fingerprints(
                    fps, k, linkage, metric, cancel=job.cancel_token
                )
        else:
            values = self.features.number_values(oids)
            remember = None

            def work(job):
                if method == "kmeans":
                    return None, cl.kmeans_numeric(values, k)
                return None, cl.hierarchical_numeric(values, k, linkage)

        def done(result) -> None:
            fps, clusters = result
            if remember is not None:
                remember(fps)
            self.finish(oids, clusters, name, overwrite)

        return ctx.run("Clustering…", work, done)

    def finish(
        self, oids: np.ndarray, clusters: cl.Clusters, name: str, overwrite: bool = True
    ) -> dict[str, str]:
        ctx = self.ctx
        member = clusters.labels >= 0
        number = np.where(member, clusters.labels + 1, np.nan)
        size = np.where(member, clusters.sizes(), np.nan)
        centre = ["yes" if c else "" for c in clusters.centres.tolist()]
        names = ctx.apply_results(
            {
                name: (oids, NumberColumn.from_floats(number, 0)),
                f"{name} Size": (oids, NumberColumn.from_floats(size, 0)),
                f"{name} Centre": (oids, centre),
            },
            label="Cluster",
            overwrite=overwrite,
        )
        left_out = int((~member).sum())
        message = (
            f"{clusters.count:,} cluster(s) in “{names[name]}”, "
            f"{clusters.singletons:,} singleton(s); cluster 1 is the largest"
        )
        if left_out:
            message += f"; {left_out:,} row(s) without a structure or value"
        ctx.show_message(message, 10_000)
        return names
