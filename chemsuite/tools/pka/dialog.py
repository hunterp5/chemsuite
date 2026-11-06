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

"""The Predict pKa dialog."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QVBoxLayout, QWidget

from chemsuite.app.dialogs import OutputChoice, RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn, StructureColumn
from chemsuite.ionisation.predict import PkaResults, predict_pka

POOL_MIN_ROWS = 200

#: output key -> (label in the dialog, tooltip)
OUTPUTS = {
    "pka": ("pKa (most acidic), pKa (most basic)", "The strongest acid and the strongest base."),
    "sites": ("pKa Values", "Every predicted site, lowest pKa first, marked acid or base."),
    "charge": ("Charge, % Neutral", "Mean net charge, and the share with no site ionised."),
    "logd": ("cLogD", "Crippen cLogP of the neutral parent, corrected for ionisation."),
    "species": ("Major Species", "The dominant protonation state, as a structure column."),
}
DEFAULT_OUTPUTS = ("pka", "charge", "logd")


def column_names(ph: float) -> dict[str, list[str]]:
    """The columns each output writes, at *ph*."""
    at = f"@ pH {ph:g}"
    return {
        "pka": ["pKa (most acidic)", "pKa (most basic)"],
        "sites": ["pKa Values"],
        "charge": [f"Charge {at}", f"% Neutral {at}"],
        "logd": [f"cLogD {ph:g}"],
        "species": [f"Major Species {at}"],
    }


class PkaDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(
            ctx.parent,
            "Predict pKa",
            guide="predict-pka",
        )
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.ph = QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setDecimals(1)
        self.ph.setSingleStep(0.1)
        self.ph.setValue(7.4)
        method.row(
            "pH",
            self.ph,
            help="Charge, % neutral, cLogD and the major species are given at this pH. "
            "The pKa values themselves do not depend on it.",
        )

        output = self.section("Output")
        boxes = QWidget()
        layout = QVBoxLayout(boxes)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        self.outputs: dict[str, QCheckBox] = {}
        for key, (label, tip) in OUTPUTS.items():
            box = QCheckBox(label)
            box.setToolTip(tip)
            box.setChecked(key in DEFAULT_OUTPUTS)
            layout.addWidget(box)
            self.outputs[key] = box
        output.row("Columns", boxes)
        self.if_exists = output.row("If they exist", OutputChoice())

        for signal in (self.ph.valueChanged, self.rows.changed):
            signal.connect(lambda *_: self.refresh())
        for box in self.outputs.values():
            box.toggled.connect(lambda *_: self.refresh())
        self.refresh()

    # -- what will happen ---------------------------------------------------------------------
    def chosen(self) -> list[str]:
        return [key for key, box in self.outputs.items() if box.isChecked()]

    def _columns(self) -> list[str]:
        names = column_names(self.ph.value())
        return [name for key in self.chosen() for name in names[key]]

    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to predict from."
        if not self._rows():
            return "There are no rows to predict."
        if not self.chosen():
            return "Choose at least one kind of column to write."
        return ""

    def summary(self) -> str:
        return f"MolGpKa at pH {self.ph.value():g} → " + ", ".join(self._columns())

    def note(self) -> str:
        if self._rows() > 50_000:
            return f"{self._rows():,} rows may take a few minutes (you can cancel)."
        return ""

    def verb(self) -> str:
        return f"Predict pKa for {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        ctx = self.ctx
        table = ctx.table
        oids = self.rows.oids()
        column = table.column(self.structure.currentText())
        slots = table.slots_of(oids)
        cells = [  # unparsed cells go as text; the workers parse them
            blob if blob is not None else (text or None)
            for blob, text in zip(column.blobs(slots), column.displays(slots), strict=True)
        ]
        ph = self.ph.value()
        chosen = self.chosen()
        overwrite = self.if_exists.overwrite()
        use_pool = len(cells) >= POOL_MIN_ROWS

        def work(job):
            return predict_pka(
                cells,
                ph,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(results: PkaResults) -> None:
            finish(ctx, oids[: len(results)], results, ph, chosen, overwrite)

        return ctx.run("Predicting pKa…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    results: PkaResults,
    ph: float,
    chosen: list[str],
    overwrite: bool,
) -> dict[str, str]:
    """Write the chosen columns for *oids* (the rows *results* covers)."""
    names = column_names(ph)
    values: dict[str, object] = {
        "pka": [
            NumberColumn.from_floats(np.array(results.most_acidic), 2),
            NumberColumn.from_floats(np.array(results.most_basic), 2),
        ],
        "sites": [results.sites],
        "charge": [
            NumberColumn.from_floats(np.array(results.charge), 2),
            NumberColumn.from_floats(np.array(results.neutral), 1),
        ],
        "logd": [NumberColumn.from_floats(np.array(results.logd), 2)],
        "species": [StructureColumn.from_parsed(results.species_blobs, results.species_texts)],
    }
    columns = {
        name: (oids, value)
        for key in chosen
        for name, value in zip(names[key], values[key], strict=True)
    }
    written = ctx.apply_results(columns, label="Predict pKa", overwrite=overwrite)
    message = f"Predicted pKa for {len(oids):,} row(s) at pH {ph:g}"
    if results.failed:
        message += f"; {results.failed:,} structure(s) could not be processed"
    ctx.show_message(message, 8000)
    return written
