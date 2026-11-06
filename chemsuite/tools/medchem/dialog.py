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
"""One dialog for the MedChem plots: compute what is missing, classify, then plot."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QVBoxLayout,
)

from chemsuite.analysis import medchem_space as ms
from chemsuite.app.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.descriptors.catalog import CATALOG
from chemsuite.descriptors.compute import compute_descriptors
from chemsuite.plotting.data import numbers
from chemsuite.tools.descriptors.dialog import POOL_MIN_ROWS
from chemsuite.tools.sar.common import activity_combo, structure_cells

LIPE_COLUMN = "LipE"


@dataclass(frozen=True)
class Kind:
    title: str
    keys: tuple[str, ...]  # descriptors the plot needs
    x: str  # descriptor key on X ("" for the activity)
    y: str
    column: str  # the classification (or LipE) column written
    about: str


KINDS = {
    "boiled_egg": Kind(
        "BOILED-Egg",
        ("TPSA_SP", "MolLogP"),
        "TPSA_SP",
        "MolLogP",
        "BOILED-Egg",
        "Compounds in the white are predicted to be absorbed from the gut; in the yolk, also to "
        "cross the blood-brain barrier (Daina & Zoete 2016). Writes a BOILED-Egg column.",
    ),
    "golden_triangle": Kind(
        "Golden Triangle",
        ("MolLogP", "MolWt"),
        "MolLogP",
        "MolWt",
        "Golden Triangle",
        "Compounds inside the triangle tend to combine good permeability with low clearance "
        "(Johnson et al. 2009; drawn with cLogP for logD). Writes a Golden Triangle column.",
    ),
    "lipe": Kind(
        "Lipophilic Efficiency",
        ("MolLogP",),
        "MolLogP",
        "",
        LIPE_COLUMN,
        "LipE = pActivity − cLogP: potency gained without adding lipophilicity. Dashed lines "
        "mark LipE 3, 5 and 7. Writes a LipE column.",
    ),
}


def title_of(key: str) -> str:
    return CATALOG[key].title


class MedChemPlotDialog(QDialog):
    def __init__(self, ctx: ToolContext, kind: str) -> None:
        super().__init__(ctx.parent)
        self.ctx = ctx
        self.kind_key = kind
        self.kind = KINDS[kind]
        self.setWindowTitle(self.kind.title)
        self.job = None
        self.plot = None
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        self.activity = activity_combo(ctx)
        self.recalculate = QCheckBox("Recalculate descriptors that are already columns")
        about = QLabel(self.kind.about)
        about.setWordWrap(True)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Plot")
        form = QFormLayout()
        form.addRow("Structures:", self.structure)
        if kind == "lipe":
            form.addRow("pActivity:", self.activity)
        form.addRow("", self.recalculate)
        layout = QVBoxLayout(self)
        layout.addWidget(about)
        layout.addLayout(form)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(480, 0)
        self.recalculate.toggled.connect(lambda *_: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    def missing(self) -> list[str]:
        table = self.ctx.table
        if self.recalculate.isChecked():
            return list(self.kind.keys)
        return [k for k in self.kind.keys if title_of(k) not in table]

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column."
        if self.kind_key == "lipe" and not self.activity.currentText():
            return "There is no number column for activity (pIC50, pKi…)."
        return ""

    def _update(self) -> None:
        problem = self.problem()
        missing = self.missing()
        note = problem or (
            f"Calculates {', '.join(title_of(k) for k in missing)} first." if missing else ""
        )
        self.status.setText(note)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if self.problem():
            return
        self.job = self.run()
        super().accept()

    def run(self):
        ctx = self.ctx
        keys = self.missing()
        if not keys:
            self.finish()
            return None
        oids = ctx.table.oids.copy()
        cells = structure_cells(ctx, self.structure.currentText(), oids)
        use_pool = len(cells) >= POOL_MIN_ROWS

        def work(job):
            return compute_descriptors(
                cells,
                keys,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(columns) -> None:
            finished = len(next(iter(columns.values())))  # fewer after a cancel
            ctx.apply_results(
                {title_of(k): (oids[:finished], columns[k]) for k in keys},
                label=f"Calculate {len(keys)} descriptor(s)",
                overwrite=True,
            )
            if finished == len(oids):  # classify and plot only a complete set
                self.finish()

        return ctx.run(f"Calculating {len(keys)} descriptor(s)…", work, done, keep_partial=True)

    def finish(self):
        """Write the classification column and open the plot."""
        ctx = self.ctx
        table = ctx.table
        kind = self.kind
        oids = table.oids.copy()
        x_name = title_of(kind.x)
        if self.kind_key == "lipe":
            y_name = self.activity.currentText()
            values = ms.lipe(numbers(table, y_name, oids), numbers(table, x_name, oids))
            result = NumberColumn.from_floats(values, 2)
        else:
            y_name = title_of(kind.y)
            x, y = numbers(table, x_name, oids), numbers(table, y_name, oids)
            result = (
                ms.boiled_egg(x, y) if self.kind_key == "boiled_egg" else ms.golden_triangle(x, y)
            )
        names = ctx.apply_results({kind.column: (oids, result)}, label=kind.title, overwrite=True)
        column = names[kind.column]
        self.plot = ctx.open_tool("scatter")
        self.plot.set_spec(x=x_name, y=y_name, color=column, guide=self.kind_key)
        ctx.show_message(self.summary(column, oids), 10_000)
        return self.plot

    def summary(self, column: str, oids: np.ndarray) -> str:
        table = self.ctx.table
        if self.kind_key == "lipe":
            values = numbers(table, column, oids)
            good = int(np.sum(values >= 5))
            return f"{good:,} compound(s) with LipE ≥ 5 (in “{column}”)"
        labels = table.column(column).displays(table.slots_of(oids))
        counts = {label: labels.count(label) for label in sorted(set(labels)) if label}
        return f"{self.kind.title}: " + ", ".join(f"{k} {v:,}" for k, v in counts.items())
