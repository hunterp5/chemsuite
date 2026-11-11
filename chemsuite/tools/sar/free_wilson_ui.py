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
"""Free-Wilson…: additive substituent contributions and the best unmade combinations."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.mols import mol_to_blob, smiles_of
from chemsuite.core.table.columns import NumberColumn, StructureColumn
from chemsuite.plotting.data import numbers
from chemsuite.sar.free_wilson import FreeWilson, free_wilson
from chemsuite.sar.rgroups import HYDROGEN, Decomposition, assemble

from .common import (
    CoreChooser,
    RerunPanel,
    ResultTable,
    activity_combo,
    decompose_task,
    number_text,
    saved_rows,
)

PREDICTED_COLUMN = "FW Predicted"
UNMADE_LIMIT = 200


def fit_summary(fit, skipped) -> str:
    text = (
        f"n = {fit.n}, parameters = {fit.parameters}, R² = {fit.r2:.3f}, "
        f"q² (LOO) = {number_text(fit.q2)}, s = {number_text(fit.s)}"
    )
    left_out = ", ".join(f"{count:,} {why}" for why, count in skipped.items() if count)
    if left_out:
        text += f"\nLeft out: {left_out}"
    return text


class FreeWilsonDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Free-Wilson Analysis", guide="free-wilson")
        self.ctx = ctx
        self.panel: FreeWilsonPanel | None = None
        inputs = self.section("Input")
        self.core = CoreChooser(ctx, inputs)
        self.activity = inputs.row("Activity", activity_combo(ctx))
        self.rows = inputs.row("Rows", RowsChoice(ctx, "shown"))

        method = self.section("Method")
        self.min_count = QSpinBox()
        self.min_count.setRange(1, 100)
        self.min_count.setValue(2)
        method.row(
            "Min. occurrences",
            self.min_count,
            help="Leave out compounds carrying a substituent seen fewer times than this: a "
            "substituent seen once is fitted exactly and its contribution means nothing.",
        )
        self.higher_is_better = QCheckBox("Higher is better")
        self.higher_is_better.setChecked(True)
        method.row(
            "Activity",
            self.higher_is_better,
            help="Which way the Best Unmade tab ranks combinations.",
        )

        output = self.section("Output")
        self.output = QLineEdit(PREDICTED_COLUMN)
        output.row("Prediction column", self.output)

        for signal in (
            self.core.core.textChanged,
            self.core.structure.currentIndexChanged,
            self.activity.currentIndexChanged,
            self.rows.changed,
        ):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.core.core.text().strip():
            return "Enter or suggest a core with the varying positions open."
        if not self.activity.currentText():
            return "There is no number column for activity."
        return self.core.problem()

    def verb(self) -> str:
        return f"Analyse {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.analyse()

    def recipe(self) -> dict:
        """Everything the analysis needs, saved with its panel to rerun it on reopening."""
        return {
            **self.core.recipe(),
            "oids": self.rows.oids().tolist(),
            "activity": self.activity.currentText(),
            "min_count": self.min_count.value(),
            "higher_is_better": self.higher_is_better.isChecked(),
        }

    def analyse(self):
        recipe = self.recipe()
        name = self.output.text().strip() or PREDICTED_COLUMN
        return _compute(
            self.ctx, recipe, lambda oids, d, model: self.finish(oids, d, model, name, recipe)
        )

    def finish(self, oids, d: Decomposition, model: FreeWilson, name: str, options: dict):
        ctx = self.ctx
        names = ctx.apply_results(
            {name: (oids, NumberColumn.from_floats(model.predicted, 2))},
            label="Free-Wilson",
        )
        self.panel = FreeWilsonPanel(ctx, oids, d, model, options)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"Free-Wilson: R² {model.fit.r2:.2f} on {model.fit.n:,} compound(s); "
            f"predictions in “{names[name]}”",
            10_000,
        )
        return self.panel


class FreeWilsonPanel(ToolPanel):
    """Results of one Free-Wilson analysis. *options* is the analysis recipe (see
    :meth:`FreeWilsonDialog.recipe`); it is saved with the panel to rerun it on reopening."""

    tool_id = "free_wilson"

    def __init__(self, ctx: ToolContext, oids, d: Decomposition, model: FreeWilson, options):
        super().__init__()
        self.ctx = ctx
        self.decomposition = d
        self.model = model
        self.recipe = options
        self.oids = np.asarray(oids)
        summary = QLabel(fit_summary(model.fit, model.skipped))
        summary.setWordWrap(True)
        reference = ", ".join(f"{p} = {r}" for p, r in model.reference.items())
        summary.setToolTip(
            "Fujita-Ban: the change from the reference substituent at each position "
            f"({reference}).\nFree-Wilson: the contribution relative to the series mean."
        )
        self.terms = ResultTable(
            ["Position", "Substituent", "N", "Fujita-Ban", "± SE", "Free-Wilson"]
        )
        self.terms.fill(
            [
                [
                    t.position,
                    t.label + (" (ref)" if t.reference else ""),
                    t.count,
                    t.fujita_ban,
                    t.error,
                    t.free_wilson,
                ]
                for t in model.terms
            ],
            model.terms,
        )
        better = options.get("higher_is_better", True)
        self.unmade = model.best_unmade(UNMADE_LIMIT, higher_is_better=better)
        self.best = ResultTable([*model.positions, "Predicted"])
        self.best.fill([[*labels_, value] for labels_, value in self.unmade], self.unmade)
        add = QPushButton("Add Selected to Table")
        add.setToolTip("Build the chosen combinations and add them as new rows.")
        add.clicked.connect(lambda: self.add_to_table())
        select = QPushButton("Select Compounds")
        select.setToolTip("Select the compounds carrying the chosen substituents.")
        select.clicked.connect(lambda: self.select_carriers())
        copy = QPushButton("Copy")
        copy.clicked.connect(lambda: self.tabs.currentWidget().findChild(ResultTable).copy_text())

        self.tabs = QTabWidget()
        terms_tab = QWidget()
        box = QVBoxLayout(terms_tab)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(self.terms)
        row = QHBoxLayout()
        row.addWidget(select)
        row.addStretch(1)
        box.addLayout(row)
        best_tab = QWidget()
        box = QVBoxLayout(best_tab)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(
            QLabel(
                f"Best of {model.combinations():,} possible combinations not in the data "
                f"({'higher' if better else 'lower'} is better)."
            )
        )
        box.addWidget(self.best)
        row = QHBoxLayout()
        row.addWidget(add)
        row.addStretch(1)
        box.addLayout(row)
        self.tabs.addTab(terms_tab, "Contributions")
        self.tabs.addTab(best_tab, "Best Unmade")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(summary)
        layout.addWidget(self.tabs, 1)
        layout.addWidget(copy)

    def title(self) -> str:
        return f"Free-Wilson: {self.recipe.get('activity', '')}"

    def select_carriers(self) -> list[int]:
        terms = [t for t in self.terms.selected_data() if t is not None]
        d = self.decomposition
        hit = np.zeros(len(d), dtype=bool)
        for t in terms:
            hit |= np.array([lab == t.label for lab in d.labels[t.position]]) & d.matched
        oids = self.oids[hit].tolist()
        if oids:
            self.ctx.selection.set(oids, "free_wilson")
        return oids

    def add_to_table(self, choices=None) -> np.ndarray:
        """Add structures for *choices* (default: the selected best-unmade rows)."""
        d = self.decomposition
        choices = choices if choices is not None else [c for c in self.best.selected_data() if c]
        blobs, texts, values, failed = [], [], [], 0
        for labels_, value in choices:
            fragments = {}
            for position, label in zip(self.model.positions, labels_, strict=True):
                fragments[position] = d.fragments.get((position, label), _hydrogen(position, label))
            for position in d.positions:  # fixed positions: their single substituent
                if position not in fragments:
                    only = d.counts(position).most_common(1)
                    label = only[0][0] if only else HYDROGEN
                    fragments[position] = d.fragments.get(
                        (position, label), _hydrogen(position, label)
                    )
            mol = assemble(d.core, fragments) if None not in fragments.values() else None
            if mol is None:
                failed += 1
                continue
            blobs.append(mol_to_blob(mol))
            texts.append(smiles_of(mol))
            values.append(value)
        if not blobs:
            self.ctx.show_message("Choose rows in Best Unmade to add (none could be built).")
            return np.empty(0, dtype=np.int64)
        structure = self.recipe.get("structure") or "Structure"
        new = self.ctx.table.add_rows(
            {
                structure: StructureColumn.from_parsed(blobs, texts),
                PREDICTED_COLUMN: NumberColumn.from_floats(np.array(values), 2),
                "Source": ["Free-Wilson (not made)"] * len(blobs),
            },
            label="Add Free-Wilson proposals",
        )
        message = f"Added {len(new):,} proposed compound(s) at the end of the table"
        if failed:
            message += f"; {failed} could not be built"
        self.ctx.show_message(message, 8000)
        return new

    def state(self) -> dict:
        return dict(self.recipe)

    def close_panel(self) -> None:
        pass


def _hydrogen(position: str, label: str) -> str | None:
    if label != HYDROGEN:
        return None
    return f"[H][*:{position.removeprefix('R')}]"


def _compute(ctx: ToolContext, recipe: dict, done):
    """Decompose and fit in the background; then ``done(oids, decomposition, model)``."""
    oids = np.asarray(recipe["oids"], dtype=np.int64)
    activity = numbers(ctx.table, recipe["activity"], oids)
    decompose = decompose_task(ctx, recipe, oids)
    min_count = int(recipe["min_count"])

    def work(job):
        d = decompose(job)
        return d, free_wilson(d, activity, min_count=min_count)

    return ctx.run("Free-Wilson analysis…", work, lambda result: done(oids, *result))


def restore(ctx: ToolContext, state: dict) -> RerunPanel:
    """Reopen a saved panel: rerun the analysis on the rows it covered (the prediction column
    stays as it was saved)."""
    recipe = {**state, "oids": saved_rows(ctx, state, state["activity"], state["structure"])}
    holder = RerunPanel(FreeWilsonPanel.tool_id, f"Free-Wilson: {recipe['activity']}", recipe)
    ctx.add_panel(holder, holder.title())
    holder.job = _compute(
        ctx,
        recipe,
        lambda oids, d, model: holder.fill(FreeWilsonPanel(ctx, oids, d, model, recipe)),
    )
    return holder
