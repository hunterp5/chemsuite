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
from chemsuite.app.dialogs import OutputChoice, RowsChoice, Segmented, ToolDialog
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
    "strip_salts": ("Strip salts and solvents", "Keep the largest organic fragment."),
    "neutralise": (
        "Neutralise charges",
        "Add or remove hydrogens so charges cancel. Quaternary nitrogens stay charged.",
    ),
    "remove_isotopes": ("Remove isotope labels", "Deuterium becomes hydrogen, ¹³C carbon, …"),
}
OPTIONAL = {
    "tautomer": (
        "Canonical tautomer",
        "One consistent tautomer per molecule (RDKit), so duplicates collapse. It is not "
        "necessarily the most stable or the bioactive form. Drawn stereocentres are kept.",
    ),
    "protonate": (
        "Protonate at pH",
        "Each ionisable site in its majority form at the pH, from MolGpKa's predicted pKa "
        "values (as the Major Species of Predict pKa).",
    ),
}


def output_name(source: str) -> str:
    return f"{source} (Standardized)"


class StandardizeDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Standardize Structures", guide="standardize-structures")
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

        optional = self.section("Optional")
        self.tautomer = QCheckBox(OPTIONAL["tautomer"][0])
        optional.row("", self.tautomer, help=OPTIONAL["tautomer"][1])
        protonation = QWidget()
        line = QHBoxLayout(protonation)
        line.setContentsMargins(0, 0, 0, 0)
        line.setSpacing(8)
        self.protonate = QCheckBox(OPTIONAL["protonate"][0])
        self.ph = QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setDecimals(1)
        self.ph.setSingleStep(0.1)
        self.ph.setValue(defaults.ph)
        self.ph.setEnabled(False)
        line.addWidget(self.protonate)
        line.addWidget(self.ph)
        line.addStretch(1)
        optional.row("", protonation, help=OPTIONAL["protonate"][1])
        self.protonate.toggled.connect(self.ph.setEnabled)

        output = self.section("Output")
        self.target = output.row(
            "Write to",
            Segmented((("new", "New column"), ("source", "Replace the structures")), "new"),
            help="Replacing keeps a row's structure when it cannot be standardized. "
            "Either way, Undo restores the table.",
        )
        self.notes = QCheckBox(f"Add a “{NOTES_COLUMN}” column")
        self.notes.setChecked(True)
        output.row(
            "",
            self.notes,
            help="What changed in each structure, and any undefined stereo to resolve "
            "before docking.",
        )
        self.if_exists = output.row("If they exist", OutputChoice())

        self.structure.currentTextChanged.connect(lambda *_: self.refresh())
        for signal in (self.ph.valueChanged, self.rows.changed, self.target.changed):
            signal.connect(lambda *_: self.refresh())
        for box in (*self.steps.values(), self.tautomer, self.protonate, self.notes):
            box.toggled.connect(lambda *_: self.refresh())
        self.refresh()

    # -- what will happen ---------------------------------------------------------------------
    def options(self) -> Options:
        return Options(
            **{key: box.isChecked() for key, box in self.steps.items()},
            tautomer=self.tautomer.isChecked(),
            protonate=self.protonate.isChecked(),
            ph=self.ph.value(),
        )

    def _replacing(self) -> bool:
        return self.target.value() == "source"

    def _columns(self) -> list[str]:
        source = self.structure.currentText()
        names = [source if self._replacing() else output_name(source)]
        if self.notes.isChecked():
            names.append(NOTES_COLUMN)
        return names

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

    def summary(self) -> str:
        steps = ", ".join(self.options().steps())
        return f"{steps[:1].upper()}{steps[1:]} → " + ", ".join(self._columns())

    def note(self) -> str:
        if self.protonate.isChecked() and self._rows() > 50_000:
            return f"{self._rows():,} rows may take a few minutes (you can cancel)."
        return ""

    def verb(self) -> str:
        return f"Standardize {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        ctx = self.ctx
        table = ctx.table
        source = self.structure.currentText()
        oids = self.rows.oids()
        column = table.column(source)
        slots = table.slots_of(oids)
        cells = [  # unparsed cells go as text; the workers parse them
            blob if blob is not None else (text or None)
            for blob, text in zip(column.blobs(slots), column.displays(slots), strict=True)
        ]
        options = self.options()
        replace = self._replacing()
        notes = self.notes.isChecked()
        overwrite = self.if_exists.overwrite()
        use_pool = len(cells) >= POOL_MIN_ROWS

        def work(job):
            return standardize_cells(
                cells,
                options,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(results: StandardizeResults) -> None:
            finish(ctx, oids[: len(results)], results, source, replace, notes, overwrite)

        return ctx.run("Standardizing structures…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    results: StandardizeResults,
    source: str,
    replace: bool,
    notes: bool,
    overwrite: bool,
) -> dict[str, str]:
    """Write the results for *oids* (the rows *results* covers)."""
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
        )
    }
    if notes:
        columns[NOTES_COLUMN] = (oids, results.notes)
    written = ctx.apply_results(
        columns, label="Standardize Structures", overwrite=overwrite or replace
    )
    message = f"Standardized {len(oids):,} row(s); {results.changed:,} changed"
    if results.failed:
        message += f"; {results.failed:,} could not be standardized"
    ctx.show_message(message, 8000)
    return written
