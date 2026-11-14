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

"""Chemistry ▸ Decompose ▸ Fragments: cut each compound at its BRICS or RECAP bonds."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QSpinBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.wording import plural
from chemsuite.sar.fragments import (
    MIN_SIZE,
    POOL_MIN_ROWS,
    FragmentRows,
    fragment_rows,
    fragment_stats,
)
from chemsuite.tools.sar.common import activity_combo, structure_cells

NO_ACTIVITY = "(none)"
#: the method's name is the prefix of its columns ("BRICS 1", "RECAP 1", …)
TITLES = {"brics": "BRICS", "recap": "RECAP"}
COUNT_SUFFIX = "Fragments"
MAX_COLUMNS = 20  # a compound cut into more pieces than this writes only the largest


def column_name(method: str, index: int) -> str:
    return f"{TITLES[method]} {index + 1}"


def count_name(method: str) -> str:
    return f"{TITLES[method]} {COUNT_SUFFIX}"


class FragmentsDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Fragments", guide="fragments")
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
        self.method = method.row(
            "Cut at",
            Segmented((("brics", "BRICS"), ("recap", "RECAP")), "brics"),
            help="BRICS cuts 16 kinds of bond and labels each open end with the environment "
            "it came from, which Recombine Fragments needs. RECAP cuts 11 retrosynthetic "
            "bonds and writes every open end as a plain *.",
        )
        self.min_size = QSpinBox()
        self.min_size.setRange(1, 30)
        self.min_size.setValue(MIN_SIZE)
        self.min_size.setSuffix(" heavy atoms")
        method.row(
            "Smallest piece",
            self.min_size,
            help="Fragments with fewer heavy atoms than this are left out. Attachment points "
            "are not counted.",
        )
        self.plain = QCheckBox("Plain attachment points")
        method.row(
            "",
            self.plain,
            help="Write every open end as * rather than a numbered BRICS label, so the same "
            "piece counts once however it was cut. Recombine Fragments needs the labels.",
        )

        output = self.section("Output")
        self.columns = QCheckBox("Fragment columns")
        self.columns.setChecked(True)
        output.row(
            "",
            self.columns,
            help="One structure column per fragment of a compound (BRICS 1, BRICS 2 …), "
            "largest first, and a count column.",
        )
        self.panel = QCheckBox("Fragment panel")
        self.panel.setChecked(True)
        output.row(
            "",
            self.panel,
            help="A panel listing each distinct fragment with how many compounds hold it, "
            "and the activity of those compounds.",
        )
        self.activity = activity_combo(ctx)
        self.activity.insertItem(0, NO_ACTIVITY)
        if self.activity.count() == 1:
            self.activity.setCurrentIndex(0)
        output.row("Activity", self.activity, help="Averaged over the compounds holding each "
                   "fragment, in the panel.")  # fmt: skip

        self.structure.currentTextChanged.connect(lambda *_: self.refresh())
        for signal in (self.rows.changed, self.method.changed, self.min_size.valueChanged):
            signal.connect(lambda *_: self.refresh())
        for box in (self.plain, self.columns, self.panel):
            box.toggled.connect(lambda *_: self.refresh())
        self.refresh()

    # -- what will happen ---------------------------------------------------------------------
    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to fragment."
        if not self._rows():
            return "There are no rows to fragment."
        if not self.columns.isChecked() and not self.panel.isChecked():
            return "Choose fragment columns, the panel, or both."
        return ""

    def verb(self) -> str:
        return f"Fragment {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        ctx = self.ctx
        source = self.structure.currentText()
        oids = self.rows.oids()
        cells = structure_cells(ctx, source, oids)
        method = self.method.value()
        min_size = self.min_size.value()
        plain = self.plain.isChecked()
        want_columns = self.columns.isChecked()
        want_panel = self.panel.isChecked()
        activity = self.activity.currentText()

        def work(job):
            return fragment_rows(
                cells,
                method,
                min_size=min_size,
                plain=plain,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(found: FragmentRows) -> None:
            finish(
                ctx,
                oids[: len(found)],
                found,
                method,
                columns=want_columns,
                panel=want_panel,
                activity=activity,
            )

        return ctx.run("Finding fragments…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    found: FragmentRows,
    method: str,
    *,
    columns: bool,
    panel: bool,
    activity: str = NO_ACTIVITY,
) -> dict[str, str]:
    """Write the fragment columns for *oids*, and open the panel."""
    written: dict[str, str] = {}
    if columns:
        wide = min(found.widest(), MAX_COLUMNS)
        results: dict[str, tuple[np.ndarray, object]] = {
            column_name(method, i): (
                oids,
                StructureColumn.from_texts([row[i] if i < len(row) else "" for row in found.rows]),
            )
            for i in range(wide)
        }
        results[count_name(method)] = (oids, [float(len(row)) for row in found.rows])
        written = ctx.apply_results(results, label="Fragments")
    if panel:
        open_panel(ctx, oids, found, method, activity)
    message = (
        f"Fragmented {plural(found.fragmented(), 'compound')} into "
        f"{plural(len({f for row in found.rows for f in row}), 'distinct fragment')}"
    )
    if found.failed:
        message += f"; {found.failed:,} could not be read"
    ctx.show_message(message)
    return written


def open_panel(
    ctx: ToolContext, oids: np.ndarray, found: FragmentRows, method: str, activity: str
) -> object:
    """Show the distinct fragments of *found* in a panel."""
    from chemsuite.plotting.data import numbers

    from .panel import FragmentPanel

    values = (
        numbers(ctx.table, activity, oids)
        if activity and activity != NO_ACTIVITY and activity in ctx.table
        else None
    )
    stats = fragment_stats(found, values)
    name = activity if values is not None else NO_ACTIVITY
    panel = FragmentPanel(ctx, oids, stats, method, name)
    ctx.add_panel(panel, panel.title())
    return panel
