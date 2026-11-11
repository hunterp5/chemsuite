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
"""Hansch Analysis…: regress activity on substituent constants from the user's own table."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)

from chemsuite.app.dialogs import RowsChoice, Segmented, ToolDialog, inline
from chemsuite.app.tools import ToolContext, ToolPanel
from chemsuite.core.table.columns import NumberColumn
from chemsuite.plotting.data import numbers
from chemsuite.sar.hansch import (
    ELECTRONIC,
    FIELDS,
    HEADER,
    STERIC,
    Hansch,
    HanschSettings,
    hansch,
    load_constants,
)

from .common import CoreChooser, RerunPanel, ResultTable, activity_combo, decompose_task, saved_rows
from .free_wilson_ui import fit_summary

PREDICTED_COLUMN = "Hansch Predicted"
TEMPLATE = f"""# Substituent constants for chemsuite's Hansch analysis.
# One row per substituent: SMILES with * at the attachment point (e.g. *Cl, *C(F)(F)F,
# *OC), or H. Leave a cell blank when a value is unknown. Lines starting with # are ignored.
# Columns: pi (hydrophobic), sigma_m / sigma_p (Hammett), f / r (Swain-Lupton field and
# resonance), es (Taft steric), mr (molar refractivity), upsilon (Charton), sigma_star (Taft).
# Hydrogen counts as 0 for pi, sigma, f and r when it is not listed.
# Fill in values from your own reference source (e.g. Hansch, Leo & Taft, Chem. Rev. 1991).
{HEADER}
"""
_last_path = ""  # the constants file used last in this session


class HanschDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Hansch Analysis", guide="hansch-analysis")
        self.ctx = ctx
        self.panel: HanschPanel | None = None
        inputs = self.section("Input")
        self.core = CoreChooser(ctx, inputs)
        self.activity = inputs.row("Activity", activity_combo(ctx))
        self.rows = inputs.row("Rows", RowsChoice(ctx, "shown"))

        constants = self.section("Constants")
        self.constants = QLineEdit(_last_path)
        self.constants.setPlaceholderText("CSV file")
        browse = QPushButton("Browse…")
        browse.setAutoDefault(False)
        browse.clicked.connect(lambda: self.browse())
        constants.row(
            "File",
            inline(self.constants, browse, stretch=0),
            help=f"A CSV of substituent constants with the columns {', '.join(FIELDS)}.",
        )
        template = QPushButton("Save Template…")
        template.setAutoDefault(False)
        template.clicked.connect(lambda: self.save_template())
        constants.row(
            "",
            inline(template),
            help="Save an empty constants file with the expected columns, to fill in.",
        )

        terms = self.section("Terms")
        self.pi = QCheckBox("π")
        self.pi.setChecked(True)
        self.pi.setToolTip("π: hydrophobicity")
        self.pi_squared = QCheckBox("π²")
        self.pi_squared.setToolTip("π²: an optimum lipophilicity")
        terms.row("Hydrophobic", inline(self.pi, self.pi_squared))
        self.electronic = QComboBox()
        for key, title in ELECTRONIC.items():
            self.electronic.addItem(title, key)
        terms.row("Electronic", self.electronic)
        self.sigma = QComboBox()
        self.sigma.addItem("para (σp)", "para")
        self.sigma.addItem("meta (σm)", "meta")
        terms.row("Hammett σ", self.sigma, help="The σ used for positions not listed as meta.")
        self.meta_positions = QLineEdit()
        self.meta_positions.setPlaceholderText("R2, R3")
        terms.row(
            "Meta positions",
            self.meta_positions,
            help="Positions meta to the reaction centre, which use σm.",
        )
        self.steric = QComboBox()
        for key, title in STERIC.items():
            self.steric.addItem(title, key)
        terms.row("Steric", self.steric)
        self.coefficients = Segmented((("sum", "Summed"), ("position", "Per Position")))
        self.per_position = self.coefficients.buttons["position"]
        terms.row(
            "Coefficients",
            self.coefficients,
            help="Summed: one coefficient per term, over all positions (classical Hansch). "
            "Per Position: one coefficient per term and position.",
        )
        self._terms = terms

        output = self.section("Output")
        self.output = QLineEdit(PREDICTED_COLUMN)
        output.row("Prediction column", self.output)

        for signal in (
            self.core.core.textChanged,
            self.core.structure.currentIndexChanged,
            self.activity.currentIndexChanged,
            self.rows.changed,
            self.constants.textChanged,
            self.pi.toggled,
            self.pi_squared.toggled,
            self.electronic.currentIndexChanged,
            self.steric.currentIndexChanged,
        ):
            signal.connect(lambda *_: self._update())
        self._update()

    def browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Substituent Constants", self.constants.text(), "CSV files (*.csv);;All files (*)"
        )
        if path:
            self.constants.setText(path)

    def save_template(self, path: str | None = None) -> str:
        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Constants Template", "substituent_constants.csv", "CSV files (*.csv)"
            )
        if path:
            Path(path).write_text(TEMPLATE, encoding="utf-8")
            self.constants.setText(path)
        return path or ""

    def settings(self) -> HanschSettings:
        meta = [p.strip().upper() for p in self.meta_positions.text().replace(";", ",").split(",")]
        return HanschSettings(
            pi=self.pi.isChecked(),
            pi_squared=self.pi_squared.isChecked(),
            electronic=self.electronic.currentData(),
            sigma={p: "meta" for p in meta if p},
            default_sigma=self.sigma.currentData(),
            steric=self.steric.currentData(),
            per_position=self.per_position.isChecked(),
        )

    def problem(self) -> str:
        if not self.core.core.text().strip():
            return "Enter or suggest a core with the varying positions open."
        if not self.activity.currentText():
            return "There is no number column for activity."
        path = self.constants.text().strip()
        if not path:
            return "Choose a CSV of substituent constants, or save a template to fill in."
        if not Path(path).is_file():
            return f"There is no file {path}."
        if not self.settings().terms():
            return "Choose at least one term."
        return self.core.problem()

    def verb(self) -> str:
        return f"Analyse {len(self.rows.oids()):,} Rows"

    def _update(self) -> None:
        sigma = self.electronic.currentData() == "sigma"
        self._terms.set_row_visible(self.sigma, sigma)
        self._terms.set_row_visible(self.meta_positions, sigma)
        self.pi_squared.setEnabled(self.pi.isChecked())
        self.refresh()

    def accept(self) -> None:
        if self._problem():
            self.refresh()
            return
        try:
            constants = load_constants(self.constants.text().strip())
        except (OSError, ValueError) as err:
            self.report(str(err))
            return
        global _last_path
        _last_path = self.constants.text().strip()
        self.job = self.analyse(constants)
        QDialog.accept(self)

    def recipe(self, constants: dict) -> dict:
        """Everything the analysis needs, saved with its panel to rerun it on reopening. The
        constants themselves are kept, so reopening does not need the CSV file."""
        return {
            **self.core.recipe(),
            "oids": self.rows.oids().tolist(),
            "activity": self.activity.currentText(),
            "constants": constants,
            "settings": asdict(self.settings()),
        }

    def analyse(self, constants: dict):
        recipe = self.recipe(constants)
        output = self.output.text().strip() or PREDICTED_COLUMN
        return _compute(
            self.ctx, recipe, lambda oids, model: self.finish(oids, model, output, recipe)
        )

    def finish(self, oids, model: Hansch, output: str, recipe: dict | None = None):
        ctx = self.ctx
        names = ctx.apply_results(
            {output: (oids, NumberColumn.from_floats(model.predicted, 2))},
            label="Hansch",
        )
        self.panel = HanschPanel(model, recipe)
        ctx.add_panel(self.panel, self.panel.title())
        ctx.show_message(
            f"Hansch: R² {model.fit.r2:.2f} on {model.fit.n:,} compound(s); "
            f"predictions in “{names[output]}”",
            10_000,
        )
        return self.panel


class HanschPanel(ToolPanel):
    """One Hansch equation. *recipe* (see :meth:`HanschDialog.recipe`) is saved with the panel
    to rerun the analysis on reopening."""

    tool_id = "hansch"

    def __init__(self, model: Hansch, recipe: dict | None = None) -> None:
        super().__init__()
        self.model = model
        self.recipe = recipe or {}
        fit = model.fit
        equation = QLabel(model.equation())
        equation.setWordWrap(True)
        font = equation.font()
        font.setPointSizeF(font.pointSizeF() * 1.15)
        equation.setFont(font)
        summary = QLabel(fit_summary(fit, model.skipped))
        summary.setWordWrap(True)
        self.coefficients = ResultTable(["Term", "Coefficient", "± SE"])
        self.coefficients.fill(
            [["Intercept", fit.intercept, fit.intercept_error]]
            + [[n, c, e] for n, c, e in zip(fit.names, fit.coefficients, fit.errors, strict=True)]
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(equation)
        layout.addWidget(summary)
        layout.addWidget(self.coefficients, 1)
        if model.missing:
            worst = ", ".join(f"{label} ({n})" for label, n in model.missing.most_common(12))
            missing = QLabel(f"Substituents without the constants needed: {worst}")
            missing.setWordWrap(True)
            layout.addWidget(missing)
        copy = QPushButton("Copy")
        copy.clicked.connect(lambda: self.coefficients.copy_text())
        layout.addWidget(copy)

    def title(self) -> str:
        return f"Hansch: {self.model.activity_name}"

    def state(self) -> dict:
        return dict(self.recipe)

    def close_panel(self) -> None:
        pass


def _compute(ctx: ToolContext, recipe: dict, done):
    """Decompose and fit in the background; then ``done(oids, model)``."""
    oids = np.asarray(recipe["oids"], dtype=np.int64)
    name = recipe["activity"]
    activity = numbers(ctx.table, name, oids)
    decompose = decompose_task(ctx, recipe, oids)
    constants = {
        str(group): {str(k): float(v) for k, v in values.items()}
        for group, values in recipe["constants"].items()
    }
    settings = HanschSettings(**recipe["settings"])

    def work(job):
        d = decompose(job)
        return hansch(d, activity, constants, settings, activity_name=name)

    return ctx.run("Hansch analysis…", work, lambda model: done(oids, model))


def restore(ctx: ToolContext, state: dict) -> RerunPanel:
    """Reopen a saved panel: rerun the analysis on the rows it covered (the prediction column
    stays as it was saved)."""
    recipe = {**state, "oids": saved_rows(ctx, state, state["activity"], state["structure"])}
    holder = RerunPanel(HanschPanel.tool_id, f"Hansch: {recipe['activity']}", recipe)
    ctx.add_panel(holder, holder.title())
    holder.job = _compute(ctx, recipe, lambda _oids, model: holder.fill(HanschPanel(model, recipe)))
    return holder
