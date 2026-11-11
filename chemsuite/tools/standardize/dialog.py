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

"""The Standardize Structures dialog."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QHBoxLayout, QWidget

from chemsuite.analysis.standardize import Options, StandardizeResults, standardize_cells
from chemsuite.app.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import StructureColumn

POOL_MIN_ROWS = 200
NOTES_COLUMN = "Standardization Notes"

#: option -> (label, help), in the order the steps run
STEPS = {
    "cleanup": (
        "Clean up",
        "Remove explicit hydrogens, disconnect metals bonded to N, O or F, draw functional "
        "groups one way (nitro, N-oxides, azides…) and put charges on the strongest acid.",
    ),
    "strip_salts": ("Strip salts", "Keep the largest organic fragment."),
    "neutralise": (
        "Neutralise",
        "Add or remove hydrogens so charges cancel. Quaternary nitrogens stay charged.",
    ),
    "remove_isotopes": ("Remove isotopes", "Deuterium becomes hydrogen, ¹³C carbon, …"),
    "tautomer": (
        "Canonical tautomer",
        "One consistent tautomer per molecule (RDKit), so duplicates collapse. It is not "
        "necessarily the most stable or the bioactive form. Drawn stereocentres are kept.",
    ),
}
PROTONATE = (
    "Protonate at pH",
    "Each ionisable site in its majority form at the pH, from MolGpKa's predicted pKa "
    "values (as the Major Species of Predict pKa).",
)
NOTES_HELP = (
    f"A “{NOTES_COLUMN}” column is always added: what changed in each structure, and any "
    "undefined stereo to resolve before docking."
)


def output_name(source: str) -> str:
    return f"{source} (Standardized)"


class StandardizeDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Standardize Structures", guide="standardize")
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        steps = self.section("Steps")
        defaults = Options()
        self.steps: dict[str, QCheckBox] = {}
        for key, (label, tip) in STEPS.items():
            box = QCheckBox(label)
            box.setChecked(getattr(defaults, key))
            self.steps[key] = steps.row("", box, help=tip)
        protonation = QWidget()
        line = QHBoxLayout(protonation)
        line.setContentsMargins(0, 0, 0, 0)
        line.setSpacing(8)
        self.protonate = QCheckBox(PROTONATE[0])
        self.ph = QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setDecimals(1)
        self.ph.setSingleStep(0.1)
        self.ph.setValue(defaults.ph)
        self.ph.setEnabled(False)
        line.addWidget(self.protonate)
        line.addWidget(self.ph)
        line.addStretch(1)
        steps.row("", protonation, help=PROTONATE[1])
        self.protonate.toggled.connect(self.ph.setEnabled)

        output = self.section("Output")
        self.target = output.row(
            "Write to",
            Segmented((("source", "Source"), ("new", "New Column")), "new"),
            help="Source replaces the structures (a row keeps its structure when it cannot be "
            "standardized). New Column writes “<column> (Standardized)”, replacing it if it "
            "exists. Either way, Undo restores the table. " + NOTES_HELP,
        )

        self.structure.currentTextChanged.connect(lambda *_: self.refresh())
        for signal in (self.ph.valueChanged, self.rows.changed, self.target.changed):
            signal.connect(lambda *_: self.refresh())
        for box in (*self.steps.values(), self.protonate):
            box.toggled.connect(lambda *_: self.refresh())
        self.refresh()

    # -- what will happen ---------------------------------------------------------------------
    def options(self) -> Options:
        return Options(
            **{key: box.isChecked() for key, box in self.steps.items()},
            protonate=self.protonate.isChecked(),
            ph=self.ph.value(),
        )

    def _replacing(self) -> bool:
        return self.target.value() == "source"

    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to standardize."
        if not self._rows():
            return "There are no rows to standardize."
        if not self.options().steps():
            return "Choose at least one step."
        return ""

    def verb(self) -> str:
        return f"Standardize {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        ctx = self.ctx
        source = self.structure.currentText()
        oids = self.rows.oids()
        cells = ctx.structure_inputs(source, oids)  # workers parse unparsed cells
        options = self.options()
        replace = self._replacing()

        def work(job):
            return standardize_cells(
                cells,
                options,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(results: StandardizeResults) -> None:
            finish(ctx, oids[: len(results)], results, source, replace)

        return ctx.run("Standardizing structures…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    results: StandardizeResults,
    source: str,
    replace: bool,
) -> dict[str, str]:
    """Write the results for *oids* (the rows *results* covers), with the notes column.

    Replacing writes into *source*; otherwise the results go to new columns.
    """
    blobs, texts = results.blobs, results.texts
    target_oids = oids
    if replace:  # a row that could not be standardized keeps its structure
        keep = np.array([bool(b) for b in blobs], dtype=bool)
        target_oids = oids[keep]
        blobs = [b for b, k in zip(blobs, keep, strict=True) if k]
        texts = [t for t, k in zip(texts, keep, strict=True) if k]
    columns: dict[str, tuple[np.ndarray, object]] = {
        source if replace else output_name(source): (
            target_oids,
            StructureColumn.from_parsed(blobs, texts),
        ),
        NOTES_COLUMN: (oids, results.notes),
    }
    written = ctx.apply_results(
        columns, label="Standardize Structures", overwrite={source} if replace else False
    )
    message = f"Standardized {len(oids):,} row(s); {results.changed:,} changed"
    if results.failed:
        message += f"; {results.failed:,} could not be standardized"
    ctx.show_message(message, 8000)
    return written
