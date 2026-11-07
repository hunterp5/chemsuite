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
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.tools import ToolContext
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
from chemsuite.tools.similarity.common import scope_combo, scope_oids

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


class HanschDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Hansch Analysis")
        self.ctx = ctx
        self.job = None
        self.panel: HanschPanel | None = None
        self.core = CoreChooser(ctx)
        self.activity = activity_combo(ctx)
        self.rows = scope_combo(ctx, "shown")
        self.constants = QLineEdit(_last_path)
        self.constants.setPlaceholderText("CSV of substituent constants")
        browse = QPushButton("Browse…")
        browse.clicked.connect(self.browse)
        template = QPushButton("Save Template…")
        template.setToolTip("Save an empty constants file with the expected columns.")
        template.clicked.connect(self.save_template)
        file_row = QHBoxLayout()
        file_row.addWidget(self.constants, 1)
        file_row.addWidget(browse)
        file_row.addWidget(template)
        self.pi = QCheckBox("π (hydrophobicity)")
        self.pi.setChecked(True)
        self.pi_squared = QCheckBox("π² (an optimum lipophilicity)")
        self.electronic = QComboBox()
        for key, title in ELECTRONIC.items():
            self.electronic.addItem(title, key)
        self.sigma = QComboBox()
        self.sigma.addItem("para (σp)", "para")
        self.sigma.addItem("meta (σm)", "meta")
        self.meta_positions = QLineEdit()
        self.meta_positions.setPlaceholderText("e.g. R2, R3")
        self.meta_positions.setToolTip("Positions meta to the reaction centre use σm.")
        self.steric = QComboBox()
        for key, title in STERIC.items():
            self.steric.addItem(title, key)
        self.per_position = QCheckBox("One coefficient per position (instead of sums)")
        self.output = QLineEdit(PREDICTED_COLUMN)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Analyse")
        self.form = QFormLayout()
        self.form.addRow("Activity:", self.activity)
        self.form.addRow("Rows:", self.rows)
        self.form.addRow("Constants:", file_row)
        self.form.addRow("Hydrophobic:", self.pi)
        self.form.addRow("", self.pi_squared)
        self.form.addRow("Electronic:", self.electronic)
        self.form.addRow("Hammett σ:", self.sigma)
        self.form.addRow("Meta positions:", self.meta_positions)
        self.form.addRow("Steric:", self.steric)
        self.form.addRow("", self.per_position)
        self.form.addRow("Prediction column:", self.output)
        layout = QVBoxLayout(self)
        layout.addWidget(self.core)
        layout.addLayout(self.form)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(640, 0)
        for signal in (
            self.core.core.textChanged,
            self.constants.textChanged,
            self.pi.toggled,
            self.electronic.currentIndexChanged,
            self.steric.currentIndexChanged,
        ):
            signal.connect(lambda *_: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
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
            return (
                f"Choose a CSV of substituent constants ({', '.join(FIELDS)}). "
                "Save Template makes an empty one to fill in."
            )
        if not Path(path).is_file():
            return f"There is no file {path}."
        if not self.settings().terms():
            return "Choose at least one term."
        return self.core.problem()

    def _update(self) -> None:
        sigma = self.electronic.currentData() == "sigma"
        self.form.setRowVisible(self.sigma, sigma)
        self.form.setRowVisible(self.meta_positions, sigma)
        self.pi_squared.setEnabled(self.pi.isChecked())
        problem = self.problem()
        self.status.setText(problem)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if self.problem():
            return
        try:
            constants = load_constants(self.constants.text().strip())
        except (OSError, ValueError) as err:
            self.status.setText(str(err))
            return
        global _last_path
        _last_path = self.constants.text().strip()
        self.job = self.analyse(constants)
        super().accept()

    def recipe(self, constants: dict) -> dict:
        """Everything the analysis needs, saved with its panel to rerun it on reopening. The
        constants themselves are kept, so reopening does not need the CSV file."""
        return {
            **self.core.recipe(),
            "oids": scope_oids(self.ctx, self.rows.currentData()).tolist(),
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


class HanschPanel(QWidget):
    """One Hansch equation. *recipe* (see :meth:`HanschDialog.recipe`) is saved with the panel
    to rerun the analysis on reopening."""

    titleChanged = Signal(str)
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
