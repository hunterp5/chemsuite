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
"""Data ▸ QSAR ▸ Refine Model…: the same recipe, fitted again on the data measured since."""

from __future__ import annotations

from PySide6.QtWidgets import QCheckBox, QComboBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs.kit import RowsChoice, Segmented, ToolDialog, inline
from chemsuite.qsar.model import Model
from chemsuite.qsar.recipe import BUDGETS

from .part import PART_KEY, QsarModels
from .work import refine_in_job


class RefineModelDialog(ToolDialog):
    """Fit a model again on its own compounds plus this table's.

    The model's recipe is kept — the same features, the same split, the same target — because
    a refit that changed the recipe would not be the same model and its history would mean
    nothing. Only the effort may change, so a model first fitted Quick can be refitted
    Thorough once the data is worth it.
    """

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Refine Model", guide="refine-model")
        self.ctx = ctx
        self.part: QsarModels = ctx.document_part(
            PART_KEY, lambda saved: QsarModels(ctx.table, saved)
        )

        data = self.section("Input")
        self.which = QComboBox()
        for entry in reversed(self.part.kept()):
            self.which.addItem(entry.title, entry.asset)
        data.row(
            "Model",
            self.which,
            "The model to fit again. Its recipe is kept; only the rows and the effort change.",
        )
        self.rows = data.row(
            "Rows",
            RowsChoice(ctx),
            "The rows to fit on, together with the compounds the model already knows. A "
            "repeated measurement replaces the old value; a new compound is added.",
        )

        method = self.section("Method")
        self.budget = method.row(
            "Effort",
            Segmented(tuple((b.key, b.title) for b in BUDGETS.values()), "standard"),
            "How hard to search this time. A model worth refitting is often worth more "
            "effort than the first exploratory fit got.",
        )
        self.replace = QCheckBox("Replace the old model")
        self.replace.setChecked(True)
        method.row(
            "Then",
            inline(self.replace, width=None),
            "Unticked, the refitted model is kept beside the old one, so the two scores can "
            "be compared before the old one goes.",
        )

        output = self.section("Output")
        self.write_prediction = QCheckBox("Prediction")
        self.write_prediction.setChecked(True)
        self.write_interval = QCheckBox("Interval")
        self.write_residual = QCheckBox("Residual")
        self.write_residual.setChecked(True)
        self.write_domain = QCheckBox("Domain")
        output.row(
            "Columns",
            inline(self.write_prediction, self.write_residual, width=None),
            "The refitted model's out-of-fold numbers for the rows it was fitted on.",
        )
        output.row("", inline(self.write_interval, self.write_domain, width=None), "")

        self.which.currentIndexChanged.connect(lambda _i: self.refresh())
        self.budget.changed.connect(lambda _v: self.refresh())
        self.rows.changed.connect(lambda _v: self.refresh())
        for box in (
            self.replace,
            self.write_prediction,
            self.write_interval,
            self.write_residual,
            self.write_domain,
        ):
            box.toggled.connect(lambda _on: self.refresh())
        self.refresh()

    def model(self) -> Model | None:
        asset = self.which.currentData()
        return None if asset is None else self.part.model(str(asset))

    def chosen(self) -> list[str]:
        boxes = {
            "prediction": self.write_prediction,
            "interval": self.write_interval,
            "residual": self.write_residual,
            "domain": self.write_domain,
        }
        return [key for key, box in boxes.items() if box.isChecked()]

    # -- what the kit asks -----------------------------------------------------------------
    def problem(self) -> str:
        model = self.model()
        if model is None:
            return "This document holds no model to refine. Fit one with Auto QSAR first."
        target = model.recipe.target
        if target not in self.ctx.table:
            return (
                f"This table has no {target} column, so there are no new measurements to "
                "refine with."
            )
        if not self.chosen():
            return "Choose at least one column to write."
        return ""

    def verb(self) -> str:
        model = self.model()
        if model is None:
            return "Refine"
        return f"Refine on {len(self.rows.oids()):,} Rows"

    def run(self):
        model = self.model()
        asset = self.which.currentData()
        if model is None or asset is None:
            return None
        return refine_in_job(
            self.ctx,
            self.part,
            str(asset),
            model,
            self.rows.oids(),
            budget=self.budget.value(),
            chosen=self.chosen(),
            replace_old=self.replace.isChecked(),
        )
