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
from chemsuite.core.jobs import Cancelled
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.wording import plural
from chemsuite.sar.fragments import (
    MIN_SIZE,
    POOL_MIN_ROWS,
    FragmentRows,
    FragmentStat,
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
        recipe = {"structure": source, "min_size": min_size, "plain": plain}
        values = activity_values(ctx, activity, oids) if want_panel else None
        work = fragments_work(ctx, cells, method, min_size, plain, values, stats=want_panel)

        def done(result: tuple[FragmentRows, list[FragmentStat] | None]) -> None:
            found, stats = result
            finish(
                ctx,
                oids[: len(found)],
                found,
                method,
                columns=want_columns,
                panel=want_panel,
                activity=activity if values is not None else NO_ACTIVITY,
                recipe=recipe,
                stats=stats,
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
    recipe: dict | None = None,
    stats: list[FragmentStat] | None = None,
) -> dict[str, str]:
    """Write the fragment columns for *oids*, and open the panel (of *stats*, made in the job;
    worked out here when not given)."""
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
        if stats is None:
            values = activity_values(ctx, activity, oids)
            activity = activity if values is not None else NO_ACTIVITY
            stats = fragment_stats(found, values)
        shown = make_panel(ctx, oids, stats, method, activity, recipe)
        ctx.add_panel(shown, shown.title())
    message = (
        f"Fragmented {plural(found.fragmented(), 'compound')} into "
        f"{plural(len({f for row in found.rows for f in row}), 'distinct fragment')}"
    )
    if found.failed:
        message += f"; {found.failed:,} could not be read"
    ctx.show_message(message)
    return written


def activity_values(ctx: ToolContext, activity: str, oids: np.ndarray) -> np.ndarray | None:
    """*activity*'s numbers for *oids*, or None when there is no activity column."""
    from chemsuite.plotting.data import numbers

    if activity and activity != NO_ACTIVITY and activity in ctx.table:
        return numbers(ctx.table, activity, oids)
    return None


def fragments_work(
    ctx: ToolContext,
    cells: list,
    method: str,
    min_size: int,
    plain: bool,
    values: np.ndarray | None,
    *,
    stats: bool = True,
):
    """``work(job) -> (found, stats)``: the fragments and, with *stats*, their table (counted in
    the job: 17k rows took a fifth of a second on the GUI thread). A cancel keeps both for
    the rows finished."""
    executor = ctx.workers_for(len(cells), POOL_MIN_ROWS)

    def table_of(found: FragmentRows) -> list[FragmentStat] | None:
        if not stats:
            return None
        return fragment_stats(found, values[: len(found)] if values is not None else None)

    def work(job):
        try:
            found = fragment_rows(
                cells,
                method,
                min_size=min_size,
                plain=plain,
                executor=executor,
                progress=job.report,
                cancel=job.cancel_token,
            )
        except Cancelled as stopped:
            if stopped.partial is None:
                raise
            raise Cancelled((stopped.partial, table_of(stopped.partial))) from None
        return found, table_of(found)

    return work


def make_panel(
    ctx: ToolContext,
    oids: np.ndarray,
    stats: list[FragmentStat],
    method: str,
    activity: str,
    recipe: dict | None = None,
):
    """A panel of the distinct fragments *stats* lists."""
    from .panel import FragmentPanel

    return FragmentPanel(ctx, oids, stats, method, activity, recipe)


def restore(ctx: ToolContext, state: dict):
    """Reopen a saved panel: fragment the rows it covered again."""
    from chemsuite.tools.sar.common import RerunPanel, saved_rows

    from .panel import FragmentPanel

    structure = state["structure"]
    method = state.get("method", "brics")
    activity = state.get("activity", NO_ACTIVITY)
    oids = np.asarray(saved_rows(ctx, state, structure), dtype=np.int64)
    min_size, plain = int(state.get("min_size", MIN_SIZE)), bool(state.get("plain", False))
    recipe = {"structure": structure, "min_size": min_size, "plain": plain}
    holder = RerunPanel(FragmentPanel.tool_id, f"{method.upper()} Fragments", state)
    ctx.add_panel(holder, holder.title())
    cells = structure_cells(ctx, structure, oids)
    values = activity_values(ctx, activity, oids)
    work = fragments_work(ctx, cells, method, min_size, plain, values)
    name = activity if values is not None else NO_ACTIVITY

    def done(result: tuple[FragmentRows, list[FragmentStat] | None]) -> None:
        _found, stats = result
        holder.fill(make_panel(ctx, oids, stats or [], method, name, recipe))

    holder.job = ctx.run("Finding fragments…", work, done)
    return holder
