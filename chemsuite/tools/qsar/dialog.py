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
"""Data ▸ QSAR ▸ Auto QSAR…: pick a target, and let the fit choose the rest."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs.kit import RowsChoice, Segmented, ToolDialog, inline
from chemsuite.plotting.data import number_columns, numbers
from chemsuite.qsar.recipe import (
    BUDGETS,
    MIN_ROWS,
    Recipe,
    default_descriptors,
    default_fingerprints,
    detect_task,
)

from .columns import write_fit
from .work import fit_in_job

#: a target with more distinct values than this is a number to predict, not a class
TASKS = (("regression", "Number"), ("classification", "Class"))


class AutoQsarDialog(ToolDialog):
    """Fit the best model the budget allows for one target column.

    The dialog asks for three things — the structures, the target, and how long it may take —
    and detects the rest. Everything a chemist might want to overrule is in More Options,
    folded, because the defaults are what should be used nine times in ten.
    """

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Auto QSAR", guide="auto-qsar")
        self.ctx = ctx

        data = self.section("Input")
        self.rows = data.row("Rows", RowsChoice(ctx), "The rows the model learns from.")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        if ctx.default_structure_column():
            self.structure.setCurrentText(ctx.default_structure_column())
        data.row("Structures", self.structure, "The structures the features are computed from.")
        self.target = QComboBox()
        self.target.addItems(self._targets())
        data.row(
            "Target",
            self.target,
            "The column to predict: a measured property or activity, or a class label.",
        )

        method = self.method = self.section("Method")
        self.task = method.row(
            "Predict",
            Segmented(TASKS),
            "A number to predict, or a class to assign. Detected from the target column; "
            "change it to treat a numeric code as a class.",
        )
        self.budget = method.row(
            "Effort",
            Segmented(tuple((b.key, b.title) for b in BUDGETS.values()), "standard"),
            "How many models to try. Quick takes under a minute on a few hundred rows; "
            "Thorough searches ten times as hard for the last little gain.",
        )
        self.split = method.row(
            "Test on",
            Segmented(
                (
                    ("scaffold", "Scaffolds"),
                    ("cluster", "Clusters"),
                    ("random", "Random"),
                    ("column", "A column"),
                ),
                "scaffold",
            ),
            "How the held-out compounds are chosen. Holding out whole scaffolds or clusters "
            "asks the honest question — can it predict a series it has not seen? A random "
            "split flatters a model that has merely memorized close analogs.",
        )
        self.split_column = method.row(
            "Ordered by",
            QComboBox(),
            "The column to split on, oldest first: the model learns from the early rows and "
            "is judged on the later ones, as it will be used.",
        )

        more = self.section("More Options", collapsed=True)
        self.descriptors = QCheckBox("Descriptors")
        self.descriptors.setChecked(True)
        self.fingerprints = QCheckBox("Fingerprints")
        self.fingerprints.setChecked(True)
        self.use_columns = QCheckBox("Numeric columns")
        more.row(
            "Features from",
            inline(self.descriptors, self.fingerprints, width=None),
            "What goes into the pool the fit chooses from. The fit keeps the few features "
            "that earn their place, so a wide pool costs time, not accuracy.",
        )
        more.row("", inline(self.use_columns, width=None), more.help(self.descriptors))
        self.confidence = QDoubleSpinBox()
        self.confidence.setRange(0.5, 0.99)
        self.confidence.setSingleStep(0.05)
        self.confidence.setDecimals(2)
        self.confidence.setValue(0.9)
        more.row(
            "Confidence",
            self.confidence,
            "How often the prediction interval, or the set of classes, should contain the "
            "truth. Checked on held-out compounds, not assumed.",
        )
        self.max_features = QSpinBox()
        self.max_features.setRange(2, 2000)
        self.max_features.setValue(200)
        more.row(
            "Most features",
            self.max_features,
            "The most features the model may end up using. Fewer is usually better: a model "
            "on 20 features is one a chemist can argue with.",
        )
        self.seed = QSpinBox()
        self.seed.setRange(0, 9999)
        more.row(
            "Seed",
            self.seed,
            "The same seed, rows and settings give the same model. Change it to see how much "
            "the result depends on where the folds fell.",
        )

        output = self.section("Output")
        self.write_prediction = QCheckBox("Prediction")
        self.write_prediction.setChecked(True)
        self.write_interval = QCheckBox("Interval")
        self.write_interval.setChecked(True)
        self.write_residual = QCheckBox("Residual")
        self.write_residual.setChecked(True)
        self.write_domain = QCheckBox("Domain")
        self.write_domain.setChecked(True)
        self.write_fold = QCheckBox("Fold")
        output.row(
            "Columns",
            inline(self.write_prediction, self.write_interval, width=None),
            "The columns written for the rows the model was fitted on. Every number is from "
            "a fold that did not include that row, so the table shows an honest estimate.",
        )
        output.row("", inline(self.write_residual, self.write_domain, width=None), "")
        output.row("", inline(self.write_fold, width=None), "")

        for widget in (self.structure, self.target, self.split_column):
            widget.currentTextChanged.connect(self._changed)
        for segment in (self.task, self.budget, self.split):
            segment.changed.connect(self._changed)
        for box in (
            self.descriptors,
            self.fingerprints,
            self.use_columns,
            self.write_prediction,
            self.write_interval,
            self.write_residual,
            self.write_domain,
            self.write_fold,
        ):
            box.toggled.connect(self._changed)
        self.rows.changed.connect(self._changed)
        self.remember(
            ctx,
            "auto_qsar",
            {
                "rows": self.rows,
                "structure": self.structure,
                "budget": self.budget,
                "split": self.split,
                "descriptors": self.descriptors,
                "fingerprints": self.fingerprints,
                "columns": self.use_columns,
                "confidence": self.confidence,
                "max_features": self.max_features,
                "seed": self.seed,
                "prediction": self.write_prediction,
                "interval": self.write_interval,
                "residual": self.write_residual,
                "domain": self.write_domain,
                "fold": self.write_fold,
            },
        )
        self._changed()

    # -- the fields ------------------------------------------------------------------------
    def _targets(self) -> list[str]:
        """Columns worth predicting: the numbers, and the short text columns that are labels.

        A structure column is never a target, and neither is a text column with a different
        value in every row — that is a name, not a class.
        """
        table = self.ctx.table
        out = list(number_columns(table))
        rows = max(len(table), 1)
        for name in table.column_names:
            if name in out or name in self.ctx.structure_columns():
                continue
            distinct = len(self.ctx.distinct_values(name, table.oids))
            if 2 <= distinct <= min(10, max(2, rows // 5)):
                out.append(name)
        return out

    def _changed(self, *_a) -> None:
        target = self.target.currentText()
        if target:
            numeric = target in number_columns(self.ctx.table)
            distinct = len(self.ctx.distinct_values(target, self.rows.oids()))
            self.task.set_value(detect_task(numeric, distinct))
            self.task.set_enabled(
                "regression", numeric, "This target's values are labels, not numbers."
            )
        orderings = self._orderings()
        self.split.set_enabled(
            "column", bool(orderings), "There is no column to order the rows by."
        )
        by_column = self.split.value() == "column"
        if (
            by_column
            and [self.split_column.itemText(i) for i in range(self.split_column.count())]
            != orderings
        ):
            self.split_column.clear()
            self.split_column.addItems(orderings)
        self.method.set_row_visible(self.split_column, by_column)
        self.refresh()

    def _orderings(self) -> list[str]:
        """Columns a time split can order by: the numbers and the dates."""
        return [n for n in number_columns(self.ctx.table) if n != self.target.currentText()]

    def chosen(self) -> list[str]:
        boxes = {
            "prediction": self.write_prediction,
            "interval": self.write_interval,
            "residual": self.write_residual,
            "domain": self.write_domain,
            "fold": self.write_fold,
        }
        return [key for key, box in boxes.items() if box.isChecked()]

    def recipe(self) -> Recipe:
        task = self.task.value()
        return Recipe(
            task=task,
            target=self.target.currentText(),
            structures=self.structure.currentText(),
            descriptors=default_descriptors() if self.descriptors.isChecked() else (),
            fingerprints=default_fingerprints() if self.fingerprints.isChecked() else (),
            columns=self._extra_columns() if self.use_columns.isChecked() else (),
            split=self.split.value(),
            split_column=self.split_column.currentText() if self.split.value() == "column" else "",
            budget=self.budget.value(),
            confidence=round(self.confidence.value(), 2),
            seed=self.seed.value(),
            max_features=self.max_features.value(),
        )

    def _extra_columns(self) -> tuple[str, ...]:
        """The table's own numeric columns, minus the target and anything derived from it."""
        target = self.target.currentText()
        return tuple(
            name
            for name in number_columns(self.ctx.table)
            if name != target and not name.startswith(f"{target} ")
        )

    # -- what the kit asks -----------------------------------------------------------------
    def problem(self) -> str:
        if not self.ctx.structure_columns():
            return "There is no structure column to compute features from."
        target = self.target.currentText()
        if not target:
            return "There is no column to predict. Add a column of measured values."
        if target not in self.ctx.table:
            # the dialog does not block the window, so the column can go while it is open
            return f"The {target} column has gone. Choose another column to predict."
        oids = self.rows.oids()
        if len(oids) < MIN_ROWS:
            return (
                f"{len(oids)} rows cannot be modeled honestly; {MIN_ROWS} is the fewest. "
                "A model needs enough compounds to be held out from."
            )
        if self.split.value() == "column" and not self.split_column.currentText():
            return "Choose the column to order the rows by."
        if not (
            self.descriptors.isChecked()
            or self.fingerprints.isChecked()
            or (self.use_columns.isChecked() and self._extra_columns())
        ):
            return "Choose at least one kind of feature to learn from."
        if not self.chosen():
            return "Choose at least one column to write."
        y = self._values(oids)
        if y is not None and np.isfinite(y).sum() < MIN_ROWS:
            return (
                f"Only {int(np.isfinite(y).sum())} of these rows have a value for "
                f"{self.target.currentText()}; {MIN_ROWS} is the fewest."
            )
        return ""

    def _values(self, oids: np.ndarray) -> np.ndarray | None:
        target = self.target.currentText()
        if not target or target not in number_columns(self.ctx.table):
            return None
        return numbers(self.ctx.table, target, oids)

    def verb(self) -> str:
        return f"Fit on {len(self.rows.oids()):,} Rows"

    def run(self):
        ctx = self.ctx
        return fit_in_job(
            ctx,
            self.recipe(),
            self.rows.oids(),
            chosen=self.chosen(),
            on_done=lambda fitted, oids: write_fit(ctx, fitted, oids, self.chosen()),
        )
