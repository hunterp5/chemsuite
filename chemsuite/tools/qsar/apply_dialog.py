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
"""Data ▸ QSAR ▸ Apply Model…: predict with a model this document holds, or one from a file."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QFileDialog, QPushButton

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs.kit import RowsChoice, ToolDialog, inline
from chemsuite.qsar.apply import ModelError
from chemsuite.qsar.model import SUFFIX, Model, load

from .columns import names_for, values_for
from .part import PART_KEY, QsarModels
from .work import missing_features, predict_in_job

FROM_FILE = "__file__"


class ApplyModelDialog(ToolDialog):
    """Predict for rows with an already fitted model.

    These numbers are not the fit's columns. The fit shows what the model would have said
    about a compound it had not seen; this shows what it says about a compound now, which for
    the rows it was trained on is its own homework and will look better than it is.
    """

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Apply Model", guide="apply-model")
        self.ctx = ctx
        self.part: QsarModels = ctx.document_part(
            PART_KEY, lambda saved: QsarModels(ctx.table, saved)
        )
        self._loaded: Model | None = None
        self._path: Path | None = None

        source = self.section("Input")
        self.which = QComboBox()
        self._fill_models()
        self.browse = QPushButton("Open…")
        self.browse.setAutoDefault(False)
        self.browse.clicked.connect(self.open_file)
        source.row(
            "Model",
            inline(self.which, self.browse, stretch=0),
            f"A model this document holds, or a {SUFFIX} file. chemsuite does not load "
            "pickled models: opening one would run whatever code it contains.",
        )
        self.rows = source.row("Rows", RowsChoice(ctx), "The rows to predict.")

        output = self.section("Output")
        self.write_prediction = QCheckBox("Prediction")
        self.write_prediction.setChecked(True)
        self.write_interval = QCheckBox("Interval")
        self.write_interval.setChecked(True)
        self.write_domain = QCheckBox("Domain")
        self.write_domain.setChecked(True)
        output.row(
            "Columns",
            inline(self.write_prediction, self.write_interval, width=None),
            "The interval, and whether the compound is close enough to the training set to "
            "be worth believing at all.",
        )
        output.row("", inline(self.write_domain, width=None), "")

        self.which.currentIndexChanged.connect(lambda _i: self.refresh())
        for box in (self.write_prediction, self.write_interval, self.write_domain):
            box.toggled.connect(lambda _on: self.refresh())
        self.rows.changed.connect(lambda _v: self.refresh())
        self.refresh()

    def _fill_models(self) -> None:
        self.which.clear()
        for entry in reversed(self.part.kept()):
            self.which.addItem(entry.title, entry.asset)
        if self._loaded is not None and self._path is not None:
            self.which.addItem(f"{self._path.name} (file)", FROM_FILE)
            self.which.setCurrentIndex(self.which.count() - 1)

    def open_file(self) -> None:
        chosen, _filter = QFileDialog.getOpenFileName(
            self, "Open Model", "", f"chemsuite model (*{SUFFIX})"
        )
        if not chosen:
            return
        try:
            self._loaded = load(Path(chosen))
        except (ModelError, OSError) as error:
            self._loaded = None
            self._path = None
            self.report(str(error))
            return
        self._path = Path(chosen)
        self._fill_models()
        self.refresh()

    def model(self) -> Model | None:
        asset = self.which.currentData()
        if asset == FROM_FILE:
            return self._loaded
        return None if asset is None else self.part.model(str(asset))

    def chosen(self) -> list[str]:
        boxes = {
            "prediction": self.write_prediction,
            "interval": self.write_interval,
            "domain": self.write_domain,
        }
        return [key for key, box in boxes.items() if box.isChecked()]

    # -- what the kit asks -----------------------------------------------------------------
    def problem(self) -> str:
        model = self.model()
        if model is None:
            return (
                "There is no model to apply. Fit one with Auto QSAR, or open a saved "
                f"{SUFFIX} file."
            )
        if not self.ctx.structure_columns():
            return "There is no structure column to compute features from."
        missing = missing_features(self.ctx, model)
        if missing:
            return (
                f"This model predicts from {', '.join(missing)}, which this document does "
                "not have. Add the column, or apply the model where it exists."
            )
        if not self.chosen():
            return "Choose at least one column to write."
        return ""

    def verb(self) -> str:
        return f"Predict {len(self.rows.oids()):,} Rows"

    def run(self):
        ctx = self.ctx
        model = self.model()
        if model is None:
            return None
        chosen = self.chosen()

        def done(said, oids: np.ndarray, usable: np.ndarray) -> None:
            write_prediction(ctx, model, said, oids, usable, chosen)

        return predict_in_job(ctx, model, self.rows.oids(), on_done=done)


def write_prediction(
    ctx: ToolContext,
    model: Model,
    said,
    oids: np.ndarray,
    usable: np.ndarray,
    chosen: list[str],
) -> None:
    """Write a prediction's columns for the rows whose structures could be read."""
    keep = np.flatnonzero(usable)
    names = names_for(model.recipe.target, chosen, model.recipe.task)
    values = values_for(model, _take(said, keep), chosen)
    fresh = ctx.table.unique_names(names)
    row_oids = np.asarray(oids, dtype=np.int64)[keep]
    written = {name: (row_oids, value) for name, value in zip(fresh, values, strict=True)}
    ctx.apply_results(written, label=f"Apply Model ({model.recipe.target})")
    outside = int((~np.asarray(said.inside, dtype=np.bool_)[keep]).sum())
    message = f"Predicted {model.recipe.target} for {len(keep):,} rows"
    skipped = len(oids) - len(keep)
    if skipped:
        message += f"; {skipped:,} structures could not be read"
    if outside:
        message += f"; {outside:,} sit outside the model's domain"
    ctx.show_message(message)


def _take(said, keep: np.ndarray):
    """A prediction narrowed to the rows that were usable."""
    from chemsuite.qsar.model import Prediction

    def pick(values):
        return None if values is None else np.asarray(values)[keep]

    return Prediction(
        np.asarray(said.predicted)[keep],
        low=pick(said.low),
        high=pick(said.high),
        scores=pick(said.scores),
        sets=pick(said.sets),
        distance=pick(said.distance),
        inside=pick(said.inside),
    )
