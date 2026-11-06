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
"""Fingerprint Similarity: every row's similarity to one structure, or to the selected rows."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QVBoxLayout,
)

from chemsuite.analysis.fingerprints import max_similarity, query_fingerprint, similarity
from chemsuite.app.dialogs import follow_rows
from chemsuite.app.tools import ToolContext
from chemsuite.core.jobs import Cancelled
from chemsuite.core.structure_input import read_structure
from chemsuite.core.table.columns import NumberColumn

from .common import FeatureChooser, fingerprint_task, labels, scope_combo, scope_oids

NEAREST_COLUMN = "Nearest Selected"


class SimilarityDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Fingerprint Similarity")
        self.ctx = ctx
        self.job = None
        self.features = FeatureChooser(ctx, numbers_allowed=False)
        self.to_structure = QRadioButton("A structure:")
        self.to_selection = QRadioButton()
        group = QButtonGroup(self)
        group.addButton(self.to_structure)
        group.addButton(self.to_selection)
        self.to_structure.setChecked(True)

        def follow() -> None:
            count = len(ctx.selection)
            self.to_selection.setText(
                f"The selected rows ({count:,}): similarity to the nearest one"
            )
            self.to_selection.setEnabled(count > 0)
            if not count and self.to_selection.isChecked():
                self.to_structure.setChecked(True)

        follow_rows(self.to_selection, ctx, follow)
        self.query = QLineEdit()
        self.query.setPlaceholderText("SMILES, or Draw… / Current Row")
        draw = QPushButton("Draw…")
        current = QPushButton("Current Row")
        self.query_info = QLabel()
        self.rows = scope_combo(ctx)
        self.output = QLineEdit("Similarity")
        self.sort = QCheckBox("Sort the table, most similar first")
        self.sort.setChecked(True)
        self.select = QCheckBox("Select rows with similarity of at least")
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.0, 1.0)
        self.threshold.setSingleStep(0.05)
        self.threshold.setValue(0.7)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Calculate")

        query_row = QHBoxLayout()
        query_row.addWidget(self.query, 1)
        query_row.addWidget(draw)
        query_row.addWidget(current)
        select_row = QHBoxLayout()
        select_row.addWidget(self.select)
        select_row.addWidget(self.threshold)
        select_row.addStretch(1)
        form = QFormLayout()
        form.addRow("Rows:", self.rows)
        form.addRow("Result column:", self.output)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Compare every row with:"))
        layout.addWidget(self.to_structure)
        layout.addLayout(query_row)
        layout.addWidget(self.query_info)
        layout.addWidget(self.to_selection)
        layout.addWidget(self.features)
        layout.addLayout(form)
        layout.addWidget(self.sort)
        layout.addLayout(select_row)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(560, 0)

        draw.clicked.connect(self.draw)
        current.clicked.connect(self.use_current_row)
        self.query.textChanged.connect(lambda _t: self._update())
        self.to_structure.toggled.connect(self._mode_changed)
        self.features.changed.connect(self._update)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        if ctx.current_oid() is not None:
            self.use_current_row()
        self._update()

    # -- the query ----------------------------------------------------------------------------
    def use_current_row(self) -> None:
        oid = self.ctx.current_oid()
        column = self.features.structure.currentText()
        if oid is None or not column:
            return
        table = self.ctx.table
        mol = table.column(column).mol(table.slots_of(np.array([oid]))[0])
        if mol is not None:
            from chemsuite.core.mols import smiles_of

            self.query.setText(smiles_of(mol))

    def draw(self) -> None:
        from chemsuite.app.structure_editor import StructureEditor

        editor = StructureEditor(self.query.text(), "Query Structure", self)
        if editor.exec() and editor.result_input.mol is not None:
            from chemsuite.core.mols import smiles_of

            self.query.setText(smiles_of(editor.result_input.mol))

    def _mode_changed(self, *_args) -> None:
        to_structure = self.to_structure.isChecked()
        self.query.setEnabled(to_structure)
        current = self.output.text()
        if current in ("Similarity", "Max Similarity to Selection"):
            self.output.setText("Similarity" if to_structure else "Max Similarity to Selection")
        self._update()

    def _query_mol(self):
        found = read_structure(self.query.text())
        return found.mol, ("" if found.mol is not None else found.error or "Enter a structure.")

    def _update(self) -> None:
        problem = self.features.problem()
        if not problem and self.to_structure.isChecked():
            _mol, problem = self._query_mol()
        self.query_info.setText(problem if self.to_structure.isChecked() else "")
        self.status.setText(problem if not self.to_structure.isChecked() else "")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    # -- running ------------------------------------------------------------------------------
    def accept(self) -> None:
        self.job = self.calculate()
        if self.job is not None:
            super().accept()

    def calculate(self):
        ctx = self.ctx
        column = self.features.structure.currentText()
        key = self.features.kind.currentData()
        metric = self.features.metric.currentData()
        targets = scope_oids(ctx, self.rows.currentData())
        name = self.output.text().strip() or "Similarity"
        if self.to_structure.isChecked():
            mol, problem = self._query_mol()
            if mol is None:
                self.query_info.setText(problem)
                return None
            query = query_fingerprint(mol, key)
            refs = np.empty(0, dtype=np.int64)
        else:
            query = None
            refs = ctx.selection.as_array()
        everything = np.union1d(targets, refs)
        compute, remember = fingerprint_task(ctx, column, key, everything)
        target_rows = np.searchsorted(everything, targets)
        ref_rows = np.searchsorted(everything, refs)

        def work(job):
            fps = compute(job)
            mine = fps.take(target_rows)
            if query is not None:
                return fps, similarity(mine, query, metric), None
            own = np.full(len(everything), -1, dtype=np.int64)
            own[ref_rows] = np.arange(len(ref_rows))  # a selected row is not its own neighbour
            try:
                best, which = max_similarity(
                    mine,
                    fps.take(ref_rows),
                    metric,
                    exclude_self=own[target_rows],
                    progress=job.report,
                    cancel=job.cancel_token,
                )
            except Cancelled as stopped:  # keep the first rows' similarities
                raise Cancelled(stopped.partial and (fps, *stopped.partial)) from None
            return fps, best, which

        def done(result) -> None:
            fps, sims, which = result
            remember(fps)
            self.finish(targets[: len(sims)], sims, which, refs, name)

        return ctx.run("Calculating similarity…", work, done, keep_partial=True)

    def finish(self, targets, sims, which, refs, name) -> None:
        ctx = self.ctx
        results = {name: (targets, NumberColumn.from_floats(sims, 3))}
        if which is not None:
            ref_labels = labels(ctx, refs)
            nearest = [ref_labels[w] if w >= 0 else "" for w in which.tolist()]
            results[NEAREST_COLUMN] = (targets, nearest)
        names = ctx.apply_results(results, label="Fingerprint similarity", overwrite=True)
        column = names[name]
        if self.sort.isChecked():
            ctx.table.set_sort([(column, False)])
        message = f"Similarity for {int(np.isfinite(sims).sum()):,} row(s) in “{column}”"
        if self.select.isChecked():
            hits = targets[np.nan_to_num(sims, nan=-1) >= self.threshold.value()]
            ctx.selection.set(hits.tolist(), "similarity")
            message += f"; selected {len(hits):,} at ≥ {self.threshold.value():.2f}"
        ctx.show_message(message, 8000)
