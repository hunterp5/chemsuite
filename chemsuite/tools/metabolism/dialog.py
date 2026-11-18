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

"""The Sites of Metabolism dialog."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QGridLayout, QVBoxLayout, QWidget

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural
from chemsuite.metabolism import COMMON_ISOFORMS, ISOFORMS
from chemsuite.metabolism.predict import SiteResults, predict_metabolism

#: output key -> (label in the dialog, tooltip)
OUTPUTS = {
    "top": (
        "Top Site, Score",
        "The atom most likely to be metabolized, as its element and index, and its score.",
    ),
    "sites": ("Ranked Sites", "The three highest-scoring atoms, best first, with their scores."),
    "domain": (
        "Applicability",
        "The share of the top atom's surroundings the model saw while training. Below 1 means "
        "the molecule is unlike the training set there, so trust the ranking less.",
    ),
}
DEFAULT_OUTPUTS = ("top", "sites")


def column_names(isoforms: list[str]) -> dict[str, list[str]]:
    """The columns each output writes, for *isoforms*."""
    return {
        "top": [name for i in isoforms for name in (f"Top Site ({i})", f"Score ({i})")],
        "sites": [f"Sites ({i})" for i in isoforms],
        "domain": [f"Applicability ({i})" for i in isoforms],
    }


class MetabolismDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Sites of Metabolism", guide="sites-of-metabolism")
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        self.isoforms: dict[str, QCheckBox] = {}
        method = self.section("Method")
        method.row(
            "Isoforms",
            self._isoform_boxes(COMMON_ISOFORMS),
            help="The isoforms to score. Each one is a column of its own; the model scores "
            "them all at once, so more of them cost little.",
        )
        more = self.section("More Options", collapsed=True)
        more.row(
            "Isoforms",
            self._isoform_boxes([i for i in ISOFORMS if i not in COMMON_ISOFORMS]),
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
            box.toggled.connect(lambda *_: self.refresh())
            layout.addWidget(box)
            self.outputs[key] = box
        output.row(
            "Columns",
            boxes,
            help="Each run writes new columns, one set per isoform. When a name is already "
            "taken, every column of the run gets the same suffix: (2), (3)…",
        )

        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def _isoform_boxes(self, isoforms: list[str] | tuple[str, ...]) -> QWidget:
        """A row of isoform checkboxes, three to a line."""
        holder = QWidget()
        grid = QGridLayout(holder)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(4)
        for position, isoform in enumerate(isoforms):
            box = QCheckBox(f"CYP{isoform}")
            box.setChecked(isoform in COMMON_ISOFORMS)
            box.toggled.connect(lambda *_: self.refresh())
            grid.addWidget(box, position // 3, position % 3)
            self.isoforms[isoform] = box
        return holder

    # -- what will happen ---------------------------------------------------------------------
    def chosen_isoforms(self) -> list[str]:
        """In the model's own order, not the order of the checkboxes."""
        return [name for name in ISOFORMS if self.isoforms[name].isChecked()]

    def chosen(self) -> list[str]:
        return [key for key, box in self.outputs.items() if box.isChecked()]

    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to predict from."
        if not self._rows():
            return "There are no rows to predict."
        if not self.chosen_isoforms():
            return "Choose at least one isoform."
        if not self.chosen():
            return "Choose at least one kind of column to write."
        return ""

    def verb(self) -> str:
        return f"Find Sites for {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        ctx = self.ctx
        oids = self.rows.oids()
        cells = ctx.structure_inputs(self.structure.currentText(), oids)
        isoforms = self.chosen_isoforms()
        chosen = self.chosen()

        def work(job):
            return predict_metabolism(cells, isoforms, progress=job.report, cancel=job.cancel_token)

        def done(results: SiteResults) -> None:
            finish(ctx, oids[: len(results)], results, isoforms, chosen)

        return ctx.run("Finding sites of metabolism…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    results: SiteResults,
    isoforms: list[str],
    chosen: list[str],
) -> dict[str, str]:
    """Write the chosen columns for *oids* (the rows *results* covers), always as new columns."""
    names = column_names(isoforms)
    values: dict[str, list[object]] = {
        "top": [
            value
            for isoform in isoforms
            for value in (
                results.top_site[isoform],
                NumberColumn.from_floats(np.array(results.top_score[isoform]), 2),
            )
        ],
        "sites": [results.sites[isoform] for isoform in isoforms],
        "domain": [
            NumberColumn.from_floats(np.array(results.domain[isoform]), 2) for isoform in isoforms
        ],
    }
    requested = [name for key in chosen for name in names[key]]
    flat = [value for key in chosen for value in values[key]]
    fresh = ctx.table.unique_names(requested)
    columns = {name: (oids, value) for name, value in zip(fresh, flat, strict=True)}
    written = ctx.apply_results(columns, label="Sites of Metabolism")
    message = (
        f"Found sites of metabolism for {plural(len(oids), 'row')}, "
        f"{plural(len(isoforms), 'isoform')}"
    )
    if results.failed:
        message += f"; {plural(results.failed, 'structure')} could not be processed"
    ctx.show_message(message)
    return written
