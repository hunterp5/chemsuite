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
"""The Generate Conformers and Export Conformers dialogs."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QLabel,
    QSpinBox,
    QVBoxLayout,
)

from chemsuite.app.tools import ToolContext
from chemsuite.conformers.etkdg import ConformerSettings, generate_conformers
from chemsuite.core.jobs.parallel import default_workers
from chemsuite.core.table.columns import NumberColumn, StructureColumn

#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 50
#: Past this many conformers in total, warn about memory and time.
WARN_CONFORMERS = 250_000
BYTES_PER_CONFORMER = 600  # measured on drug-like molecules, hydrogens included
MOLECULES_PER_SECOND = 2  # per worker process, ten conformers each
SUFFIX = " 3D"
COUNT_COLUMN = "Conformers"
SPAN_COLUMN = "Conformer Energy Range"
FORCE_FIELDS = {"mmff": "MMFF94 (UFF where MMFF lacks parameters)", "uff": "UFF", None: "None"}

_last = ConformerSettings()  # remembered for the session


def output_name(structure: str) -> str:
    return structure if structure.endswith(SUFFIX) else structure + SUFFIX


class ConformerDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Generate Conformers")
        self.ctx = ctx
        self.job = None

        self.structure = QComboBox()
        # Generating from a 3D column again would only re-embed; prefer the 2D source.
        names = ctx.structure_columns()
        self.structure.addItems(names)
        default = ctx.default_structure_column()
        if default and default.endswith(SUFFIX) and default[: -len(SUFFIX)] in names:
            default = default[: -len(SUFFIX)]
        if default:
            self.structure.setCurrentText(default)
        self.count = QSpinBox()
        self.count.setRange(1, 300)
        self.count.setToolTip("At most this many per molecule; rigid molecules get fewer.")
        self.prune = QDoubleSpinBox()
        self.prune.setRange(0.0, 5.0)
        self.prune.setSingleStep(0.1)
        self.prune.setSuffix(" Å")
        self.prune.setToolTip(
            "Conformers closer than this heavy-atom RMSD count as one (0: keep all)."
        )
        self.force_field = QComboBox()
        for key, title in FORCE_FIELDS.items():
            self.force_field.addItem(title, key)
        self.use_window = QCheckBox("Drop conformers more than")
        self.window = QDoubleSpinBox()
        self.window.setRange(0.5, 200.0)
        self.window.setSuffix(" kcal/mol above the lowest")
        self.seed = QSpinBox()
        self.seed.setRange(-1, 2**31 - 1)
        self.seed.setSpecialValueText("Random")
        self.seed.setToolTip("The same seed gives the same conformers every time.")
        self.selected_only = QCheckBox(f"Selected rows only ({len(ctx.selection):,})")
        self.selected_only.setEnabled(len(ctx.selection) > 0)
        self.overwrite = QCheckBox("Replace an existing 3D column of the same name")
        self.overwrite.setChecked(True)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Generate")

        form = QFormLayout()
        form.addRow("Structures:", self.structure)
        form.addRow("Conformers per molecule:", self.count)
        form.addRow("Distinct if RMSD above:", self.prune)
        form.addRow("Optimise with:", self.force_field)
        form.addRow(self.use_window, self.window)
        form.addRow("Random seed:", self.seed)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.selected_only)
        layout.addWidget(self.overwrite)
        layout.addWidget(self.summary)
        layout.addStretch(1)
        layout.addWidget(self.buttons)

        self.set_settings(_last)
        for signal in (
            self.structure.currentIndexChanged,
            self.count.valueChanged,
            self.selected_only.toggled,
            self.force_field.currentIndexChanged,
            self.use_window.toggled,
        ):
            signal.connect(lambda *_: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    # -- settings -----------------------------------------------------------------------------
    def set_settings(self, settings: ConformerSettings) -> None:
        self.count.setValue(settings.count)
        self.prune.setValue(settings.prune_rms)
        self.force_field.setCurrentIndex(max(0, self.force_field.findData(settings.optimize)))
        self.use_window.setChecked(settings.energy_window is not None)
        self.window.setValue(settings.energy_window or 10.0)
        self.seed.setValue(settings.seed)

    def settings(self) -> ConformerSettings:
        optimize = self.force_field.currentData()
        return ConformerSettings(
            count=self.count.value(),
            prune_rms=self.prune.value(),
            optimize=optimize,
            energy_window=self.window.value() if self.use_window.isChecked() and optimize else None,
            seed=self.seed.value(),
        )

    def _rows(self) -> int:
        return len(self.ctx.selection) if self.selected_only.isChecked() else len(self.ctx.table)

    def _update(self) -> None:
        optimize = self.force_field.currentData() is not None
        self.use_window.setEnabled(optimize)
        self.window.setEnabled(optimize and self.use_window.isChecked())
        rows = self._rows()
        structure = self.structure.currentText()
        text = f"Up to {self.count.value() * rows:,} conformers for {rows:,} row(s), "
        text += f"stored in “{output_name(structure)}”." if structure else ""
        total = self.count.value() * rows
        if total > WARN_CONFORMERS:
            megabytes = total * BYTES_PER_CONFORMER / 1e6
            hours = rows / (MOLECULES_PER_SECOND * default_workers()) / 3600
            text += (
                f"<br><b>That is a lot of 3D data</b>: up to about {megabytes:,.0f} MB of memory "
                f"and roughly {hours:.1f} hour(s). Consider fewer conformers or selected rows."
            )
        self.summary.setText(text)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            rows > 0 and bool(structure)
        )

    # -- running ------------------------------------------------------------------------------
    def accept(self) -> None:
        if not self._rows() or not self.structure.currentText():
            return
        global _last
        _last = self.settings()
        self.job = self.generate(_last)
        super().accept()

    def generate(self, settings: ConformerSettings):
        ctx = self.ctx
        table = ctx.table
        name = self.structure.currentText()
        oids = ctx.target_oids(self.selected_only.isChecked())
        column = table.column(name)
        slots = table.slots_of(oids)
        cells = [
            blob if blob is not None else (text or None)
            for blob, text in zip(column.blobs(slots), column.displays(slots), strict=True)
        ]
        overwrite = self.overwrite.isChecked()
        use_pool = len(cells) >= POOL_MIN_ROWS
        target = output_name(name)

        def work(job):
            return generate_conformers(
                cells,
                settings,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(batch) -> None:
            results = {
                target: (oids, StructureColumn.from_parsed(batch.blobs, batch.displays)),
                COUNT_COLUMN: (oids, NumberColumn.from_floats(batch.counts, 0)),
            }
            if settings.optimize:
                span = SPAN_COLUMN + " (kcal/mol)"
                results[span] = (oids, NumberColumn.from_floats(batch.spans, 2))
            names = ctx.apply_results(results, label="Generate conformers", overwrite=overwrite)
            made = int(batch.counts[batch.counts > 0].sum())
            message = f"{made:,} conformers for {len(oids):,} row(s) in “{names[target]}”"
            if batch.failed:
                message += f"; {batch.failed:,} structure(s) could not be embedded"
            ctx.show_message(message, 8000)

        return ctx.run(f"Generating conformers for {len(cells):,} row(s)…", work, done)


class ExportConformersDialog(QDialog):
    """Pick the 3D column and rows, then a file: every conformer becomes an SD record."""

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Export Conformers to SDF")
        self.ctx = ctx
        self.job = None
        self.structure = QComboBox()
        names = ctx.structure_columns()
        self.structure.addItems(names)
        three_d = [n for n in names if n.endswith(SUFFIX)]
        default = three_d[-1] if three_d else ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        self.selected_only = QCheckBox(f"Selected rows only ({len(ctx.selection):,})")
        self.selected_only.setEnabled(len(ctx.selection) > 0)
        self.selected_only.setChecked(len(ctx.selection) > 0)
        note = QLabel(
            "Each conformer is written as its own record, lowest energy first, with a "
            "“Conformer” field and the row's other columns as data fields."
        )
        note.setWordWrap(True)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Export…")
        form = QFormLayout()
        form.addRow("Structures:", self.structure)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.selected_only)
        layout.addWidget(note)
        layout.addWidget(buttons)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

    def accept(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Export Conformers", "", "SD files (*.sdf)")
        if not path:
            return
        if not path.lower().endswith(".sdf"):
            path += ".sdf"
        self.job = self.export(Path(path))
        super().accept()

    def export(self, path: Path):
        from chemsuite.core.io.sdf import write_sdf

        ctx = self.ctx
        oids = ctx.target_oids(self.selected_only.isChecked())
        structure = self.structure.currentText()
        snapshot = ctx.table.snapshot()

        def work(job) -> None:
            write_sdf(
                snapshot,
                path,
                oids=oids,
                structure=structure,
                all_conformers=True,
                progress=job.report,
                cancel=job.cancel_token,
            )

        return ctx.run(
            f"Exporting {path.name}…",
            work,
            lambda _r: ctx.show_message(f"Exported conformers to {path.name}"),
        )
