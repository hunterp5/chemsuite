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
from PySide6.QtWidgets import (
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

from chemsuite.analysis import clustering as cl
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn

from .common import FeatureChooser, fingerprint_task, scope_combo, scope_oids

CLUSTER_COLUMN = "Cluster"


class ClusterDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Cluster")
        self.ctx = ctx
        self.job = None
        self.features = FeatureChooser(ctx)
        self.method = QComboBox()
        for key, title in cl.METHODS.items():
            self.method.addItem(title, key)
        self.cutoff = QDoubleSpinBox()
        self.cutoff.setRange(0.05, 1.0)
        self.cutoff.setSingleStep(0.05)
        self.cutoff.setValue(0.6)
        self.cutoff.setToolTip(
            "Compounds at least this similar are neighbours. Higher: smaller, tighter clusters."
        )
        self.k = QSpinBox()
        self.k.setRange(2, 10_000)
        self.k.setValue(10)
        self.linkage = QComboBox()
        for key, title in cl.LINKAGES.items():
            self.linkage.addItem(title, key)
        self.rows = scope_combo(ctx)
        self.output = QLineEdit(CLUSTER_COLUMN)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Cluster")
        self.form = QFormLayout()
        self.form.addRow("Method:", self.method)
        self.form.addRow("Similarity cutoff:", self.cutoff)
        self.form.addRow("Clusters:", self.k)
        self.form.addRow("Linkage:", self.linkage)
        self.form.addRow("Rows:", self.rows)
        self.form.addRow("Result column:", self.output)
        layout = QVBoxLayout(self)
        layout.addWidget(self.features)
        layout.addLayout(self.form)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(520, 0)
        for signal in (
            self.method.currentIndexChanged,
            self.rows.currentIndexChanged,
            self.linkage.currentIndexChanged,
            self.features.changed,
        ):
            signal.connect(lambda *_: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    def _rows(self) -> int:
        return len(scope_oids(self.ctx, self.rows.currentData()))

    def problem(self) -> str:
        problem = self.features.problem()
        method = self.method.currentData()
        fingerprints = self.features.fingerprints()
        if problem:
            return problem
        if method == "butina" and not fingerprints:
            return "Butina works on fingerprints; choose k-means or hierarchical for numbers."
        if method == "hierarchical" and self._rows() > cl.MAX_HIERARCHICAL:
            return (
                f"Hierarchical clustering handles up to {cl.MAX_HIERARCHICAL:,} rows "
                f"(here {self._rows():,}). Use Butina or k-means, or fewer rows."
            )
        if method == "hierarchical" and fingerprints and self.linkage.currentData() == "ward":
            return "Ward linkage needs number columns: choose average or complete."
        return ""

    def _update(self) -> None:
        method = self.method.currentData()
        self.form.setRowVisible(self.cutoff, method == "butina")
        self.form.setRowVisible(self.k, method != "butina")
        self.form.setRowVisible(self.linkage, method == "hierarchical")
        problem = self.problem()
        note = problem
        if not problem and method == "butina" and self._rows() > cl.BUTINA_WARN_ROWS:
            note = (
                f"Butina compares every pair: {self._rows():,} rows may take minutes "
                "(it can be cancelled)."
            )
        self.status.setText(note)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if self.problem():
            return
        self.job = self.cluster()
        super().accept()

    def cluster(self):
        ctx = self.ctx
        oids = scope_oids(ctx, self.rows.currentData())
        method = self.method.currentData()
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
            self.finish(oids, clusters, name)

        return ctx.run("Clustering…", work, done)

    def finish(self, oids: np.ndarray, clusters: cl.Clusters, name: str) -> dict[str, str]:
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
            overwrite=True,
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
