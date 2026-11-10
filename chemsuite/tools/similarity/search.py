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
    QCheckBox,
    QDoubleSpinBox,
    QLineEdit,
    QPushButton,
    QToolButton,
)

from chemsuite.analysis.fingerprints import max_similarity, query_fingerprint, similarity
from chemsuite.app.dialogs import RowsChoice, Segmented, ToolDialog, follow_rows, inline
from chemsuite.app.glyphs import glyph
from chemsuite.app.tools import ToolContext
from chemsuite.core.jobs import Cancelled
from chemsuite.core.structure_input import read_structure
from chemsuite.core.table.columns import NumberColumn

from .common import FeatureChooser, fingerprint_task, labels

NEAREST_COLUMN = "Nearest Selected"
STRUCTURE_COLUMN = "Similarity"
SELECTION_COLUMN = "Max Similarity to Selection"


class SimilarityDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Similarity Search", guide="similarity-search")
        self.ctx = ctx

        inputs = self.section("Input")
        self.features = FeatureChooser(ctx, numbers_allowed=False, section=inputs)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        reference = self.section("Compare With")
        self.reference = Segmented((("structure", "A Structure"), ("selection", "Selected Rows")))
        self.to_structure = self.reference.buttons["structure"]
        self.to_selection = self.reference.buttons["selection"]
        reference.row(
            "Reference",
            self.reference,
            help="A Structure: every row's similarity to it. Selected Rows: every row's "
            "similarity to the nearest selected compound, and a column naming that compound.",
        )
        follow_rows(
            self.to_selection,
            ctx,
            lambda: self.reference.set_text("selection", f"Selected Rows ({len(ctx.selection):,})"),
        )
        self.query = QLineEdit()
        self.query.setPlaceholderText("SMILES")
        draw = QToolButton()
        draw.setAutoRaise(True)
        draw.setFixedSize(22, 22)
        draw.setIcon(glyph("pencil", draw.palette().text().color()))
        draw.setToolTip("Draw the structure in the sketcher")
        current = QPushButton("Current Row")
        current.setAutoDefault(False)
        current.setToolTip("Use the structure in the table's current row")
        self.query_field = inline(self.query, draw, current, stretch=0)
        reference.row("Structure", self.query_field, help="Type or paste SMILES, or draw it.")
        self._reference = reference

        output = self.section("Output")
        self.output = QLineEdit(STRUCTURE_COLUMN)
        output.row("Column name", self.output)
        self.sort = QCheckBox("Most similar first")
        self.sort.setChecked(True)
        output.row("Sort", self.sort, help="Sort the table by the new column.")
        self.select = QCheckBox("At least")
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.0, 1.0)
        self.threshold.setSingleStep(0.05)
        self.threshold.setValue(0.7)
        self.threshold.setEnabled(False)
        self.select.toggled.connect(self.threshold.setEnabled)
        output.row(
            "Select",
            inline(self.select, self.threshold),
            help="Select the rows at least this similar.",
        )

        draw.clicked.connect(self.draw)
        current.clicked.connect(self.use_current_row)
        self.query.textChanged.connect(lambda _t: self.refresh())
        self.reference.changed.connect(self._mode_changed)
        self.features.changed.connect(self.refresh)
        self.rows.changed.connect(lambda *_: self.refresh())
        if ctx.current_oid() is not None:
            self.use_current_row()
        self._mode_changed()

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
        from chemsuite.app.sketcher import SketcherDialog
        from chemsuite.core.mols import smiles_of

        dialog = SketcherDialog(self.query.text(), "Query Structure", self)
        if not dialog.exec():
            return
        found = read_structure(dialog.molblock) if dialog.molblock else None
        if found is not None and found.mol is not None:
            self.query.setText(smiles_of(found.mol))

    def _mode_changed(self, *_args) -> None:
        to_structure = self.to_structure.isChecked()
        self._reference.set_row_visible(self.query_field, to_structure)
        current = self.output.text()
        if current in (STRUCTURE_COLUMN, SELECTION_COLUMN):
            self.output.setText(STRUCTURE_COLUMN if to_structure else SELECTION_COLUMN)
        self.refresh()

    def _query_mol(self):
        found = read_structure(self.query.text())
        return found.mol, ("" if found.mol is not None else found.error or "Enter a structure.")

    # -- what will happen ---------------------------------------------------------------------
    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        if self.to_structure.isChecked():
            return self._query_mol()[1]
        if not len(self.ctx.selection):
            return "No rows are selected. Select the compounds to compare with."
        return ""

    def verb(self) -> str:
        return f"Compare {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.calculate()

    # -- running ------------------------------------------------------------------------------
    def calculate(self):
        ctx = self.ctx
        column = self.features.structure.currentText()
        key = self.features.kind.currentData()
        metric = self.features.metric.currentData()
        targets = self.rows.oids()
        name = self.output.text().strip() or STRUCTURE_COLUMN
        if self.to_structure.isChecked():
            mol, _problem = self._query_mol()
            if mol is None:
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
        names = ctx.apply_results(results, label="Fingerprint similarity")
        column = names[name]
        if self.sort.isChecked():
            ctx.table.set_sort([(column, False)])
        message = f"Similarity for {int(np.isfinite(sims).sum()):,} row(s) in “{column}”"
        if self.select.isChecked():
            hits = targets[np.nan_to_num(sims, nan=-1) >= self.threshold.value()]
            ctx.selection.set(hits.tolist(), "similarity")
            message += f"; selected {len(hits):,} at ≥ {self.threshold.value():.2f}"
        ctx.show_message(message, 8000)
