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
    QDoubleSpinBox,
    QFileDialog,
    QSpinBox,
)

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog, inline
from chemsuite.conformers.etkdg import ConformerSettings, generate_conformers
from chemsuite.core.table.columns import NumberColumn, StructureColumn
from chemsuite.core.wording import plural

#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 50
#: Past this many conformers in total, warn about memory and time.
SUFFIX = " 3D"
COUNT_COLUMN = "Conformers"
SPAN_COLUMN = "Conformer Energy Range"
FORCE_FIELDS = {"mmff": "MMFF94 (UFF where MMFF lacks parameters)", "uff": "UFF", None: "None"}

_last = ConformerSettings()  # remembered for the session


def output_name(structure: str) -> str:
    return structure if structure.endswith(SUFFIX) else structure + SUFFIX


class ConformerDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Generate Conformers", guide="generate-conformers")
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        # Generating from a 3D column again would only re-embed; prefer the 2D source.
        names = ctx.structure_columns()
        self.structure.addItems(names)
        default = ctx.default_structure_column()
        if default and default.endswith(SUFFIX) and default[: -len(SUFFIX)] in names:
            default = default[: -len(SUFFIX)]
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.count = QSpinBox()
        self.count.setRange(1, 300)
        method.row(
            "Conformers",
            self.count,
            help="At most this many per molecule; rigid molecules get fewer.",
        )
        self.prune = QDoubleSpinBox()
        self.prune.setRange(0.0, 5.0)
        self.prune.setSingleStep(0.1)
        self.prune.setSuffix(" Å")
        method.row(
            "Distinct above",
            self.prune,
            help="Conformers closer than this heavy-atom RMSD count as one (0: keep all).",
        )
        self.force_field = QComboBox()
        for key, title in FORCE_FIELDS.items():
            self.force_field.addItem(title, key)
        method.row("Optimise with", self.force_field)
        self.use_window = QCheckBox()
        self.window = QDoubleSpinBox()
        self.window.setRange(0.5, 200.0)
        self.window.setSuffix(" kcal/mol")
        method.row(
            "Energy window",
            inline(self.use_window, self.window),
            help="Drop conformers more than this above the molecule's lowest energy "
            "(needs a force field).",
        )

        more = self.section("More Options", collapsed=True)
        self.seed = QSpinBox()
        self.seed.setRange(-1, 2**31 - 1)
        self.seed.setSpecialValueText("Random")
        more.row("Seed", self.seed, help="The same seed gives the same conformers every time.")

        self.set_settings(_last)
        for signal in (
            self.structure.currentIndexChanged,
            self.count.valueChanged,
            self.rows.changed,
            self.force_field.currentIndexChanged,
            self.use_window.toggled,
        ):
            signal.connect(lambda *_: self._update())
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

    def _update(self) -> None:
        optimize = self.force_field.currentData() is not None
        self.use_window.setEnabled(optimize)
        self.window.setEnabled(optimize and self.use_window.isChecked())
        self.refresh()

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column."
        if not len(self.rows.oids()):
            return "There are no rows to work on."
        return ""

    def verb(self) -> str:
        return f"Generate for {len(self.rows.oids()):,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        global _last
        _last = self.settings()
        return self.generate(_last)

    def generate(self, settings: ConformerSettings):
        ctx = self.ctx
        name = self.structure.currentText()
        oids = self.rows.oids()
        cells = ctx.structure_inputs(name, oids)  # workers parse unparsed cells
        target = output_name(name)

        def work(job):
            return generate_conformers(
                cells,
                settings,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(batch) -> None:
            done_oids = oids[: len(batch.blobs)]  # after a cancel: only the first rows
            results = {
                target: (done_oids, StructureColumn.from_parsed(batch.blobs, batch.displays)),
                COUNT_COLUMN: (done_oids, NumberColumn.from_floats(batch.counts, 0)),
            }
            if settings.optimize:
                span = SPAN_COLUMN + " (kcal/mol)"
                results[span] = (done_oids, NumberColumn.from_floats(batch.spans, 2))
            names = ctx.apply_results(results, label="Generate conformers")
            made = int(batch.counts[batch.counts > 0].sum())
            message = (
                f"Generated {plural(made, 'conformer')} for "
                f"{plural(len(done_oids), 'row')} in “{names[target]}”"
            )
            if batch.failed:
                message += f"; {plural(batch.failed, 'structure')} could not be embedded"
            ctx.show_message(message)

        return ctx.run("Generating conformers…", work, done, keep_partial=True)


class ExportConformersDialog(ToolDialog):
    """Pick the 3D column and rows, then a file: every conformer becomes an SD record."""

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Export Conformers to SDF", guide="export-conformers-sdf")
        self.ctx = ctx
        inputs = self.section("Input")
        self.structure = QComboBox()
        names = ctx.structure_columns()
        self.structure.addItems(names)
        three_d = [n for n in names if n.endswith(SUFFIX)]
        default = three_d[-1] if three_d else ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row(
            "Structures",
            self.structure,
            help="Each conformer is written as its own record, lowest energy first, with a "
            "“Conformer” field and the row's other columns as data fields.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx, "selected" if len(ctx.selection) else "all"))
        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column."
        return ""

    def verb(self) -> str:
        return f"Export {len(self.rows.oids()):,} Rows…"

    def accept(self) -> None:
        if self._problem():
            self.refresh()
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Conformers", "", "SD files (*.sdf)")
        if not path:
            return
        if not path.lower().endswith(".sdf"):
            path += ".sdf"
        self.job = self.export(Path(path))
        QDialog.accept(self)

    def export(self, path: Path):
        from chemsuite.io.sdf import write_sdf

        ctx = self.ctx
        oids = self.rows.oids()
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
