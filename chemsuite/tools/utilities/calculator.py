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

"""Tools ▸ Table ▸ Calculator…: a new column from a formula over other columns.

The formula language is :mod:`chemsuite.analysis.formula`. The dialog checks the formula as
it is typed and previews the result on the first rows; Insert menus add column names and
functions at the cursor. The whole column is computed in the background in chunks, so a
cancelled run keeps the rows it finished (formulas with statistics such as ``mean`` need
every row at once and are computed in one go).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QFontDatabase, QTextCursor
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QPlainTextEdit,
    QPushButton,
    QWidget,
)

from chemsuite.analysis import formula as fx
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import FIELD_WIDTH, LABEL_WIDTH, RowsChoice, ToolDialog
from chemsuite.core.jobs import Cancelled
from chemsuite.core.table import NUMBER
from chemsuite.core.table.columns import STRUCTURE
from chemsuite.core.table.filters import number_values

PREVIEW_ROWS = 5
CHUNK = 2000


def column_values(table, name: str, oids: np.ndarray, *, blobs: bool = False) -> fx.Value:
    """A column's cells for *oids* as a formula value (read on the GUI thread).

    With *blobs*, structures are their stored bytes, turned into molecules by
    :func:`materialise` in the background job.
    """
    column = table.column(name)
    slots = table.slots_of(oids)
    if column.kind == NUMBER:
        return fx.Value(fx.NUM, number_values(column)[slots])
    if column.kind == STRUCTURE:
        if blobs:
            for slot in slots[np.equal(column.blobs(slots), None)]:
                column.mol(int(slot))  # parse cells not parsed yet
            return fx.Value(BLOBS, np.array(column.blobs(slots), dtype=object))
        return fx.Value(fx.MOL, np.array([column.mol(int(s)) for s in slots], dtype=object))
    texts = column.displays(slots)
    return fx.Value(fx.TEXT, np.array([t or None for t in texts], dtype=object))


BLOBS = "blobs"


def materialise(data: dict[str, fx.Value]) -> dict[str, fx.Value]:
    """Structure bytes to molecules (in the worker thread)."""
    from chemsuite.core.mols import mol_from_blob

    out = {}
    for name, value in data.items():
        if value.kind == BLOBS:
            mols = [mol_from_blob(b) if b else None for b in value.data]  # type: ignore[union-attr]
            value = fx.Value(fx.MOL, np.array(mols, dtype=object))
        out[name] = value
    return out


def _slice(value: fx.Value, start: int, stop: int) -> fx.Value:
    data = value.data
    return fx.Value(value.kind, data[start:stop] if isinstance(data, np.ndarray) else data)


def evaluate(
    formula: fx.Formula, data: dict[str, fx.Value], n: int, job=None
) -> tuple[str, list[object]]:
    """Evaluate for *n* rows, in chunks when possible so a cancelled job keeps what it did."""
    if formula.aggregate or n <= CHUNK:
        return fx.result_values(formula.evaluate(data, n), n)
    kinds: set[str] = set()
    cells: list[object] = []
    for start in range(0, n, CHUNK):
        if job is not None and job.cancel_token.cancelled:
            raise Cancelled(partial=(_kind(kinds), cells))
        stop = min(n, start + CHUNK)
        chunk = {name: _slice(v, start, stop) for name, v in data.items()}
        kind, part = fx.result_values(formula.evaluate(chunk, stop - start), stop - start)
        kinds.add(kind)
        cells.extend(part)
        if job is not None:
            job.report(stop, n)
    return _kind(kinds), cells


def _kind(kinds: set[str]) -> str:
    return "text" if "text" in kinds else "number"


class CalculatorDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Calculator", guide="calculator")
        self.ctx = ctx
        self.formula: fx.Formula | None = None
        self.error = ""

        expression = self.section("Formula")
        self.edit = QPlainTextEdit()
        mono = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        self.edit.setFont(mono)
        self.edit.setPlaceholderText("e.g. 9 - log10([IC50 (nM)])   or   MW([Structure]) / 100")
        self.edit.setTabChangesFocus(True)
        self.edit.setFixedWidth(FIELD_WIDTH)
        self.edit.setFixedHeight(self.edit.fontMetrics().lineSpacing() * 4 + 12)
        expression.row(
            "Formula",
            self.edit,
            "Columns in brackets ([IC50 (nM)]) or by name (MW); + − * / ^, comparisons, "
            "and/or/not, if(condition, then, else). Blank cells give blank results.",
        )
        inserts = QWidget()
        row = QHBoxLayout(inserts)
        row.setContentsMargins(0, 0, 0, 0)
        self.column_button = QPushButton("Column")
        self.function_button = QPushButton("Function")
        for button in (self.column_button, self.function_button):
            button.setAutoDefault(False)
            row.addWidget(button)
        row.addStretch(1)
        self.column_menu = QMenu(self.column_button)
        self.column_menu.aboutToShow.connect(self._fill_columns)
        self.column_button.setMenu(self.column_menu)
        self.function_button.setMenu(self._function_menu())
        expression.row("Insert", inserts, "Add a column name or function at the cursor.")

        self.preview = QLabel()
        self.preview.setTextFormat(Qt.TextFormat.PlainText)
        self.preview.setFont(mono)
        self.preview.setFixedWidth(FIELD_WIDTH)
        self.preview.setWordWrap(True)
        self.preview.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        expression.row("Preview", self.preview, "The result on the first shown rows.")

        rows = self.section("Rows")
        self.rows = RowsChoice(ctx, "all")
        rows.row("Calculate for", self.rows, "Rows outside the choice are left blank.")

        output = self.section("Output")
        self.name = QLineEdit("Result")
        output.row("Column", self.name, "The new column's name.")

        self.edit.textChanged.connect(self._changed)
        self.rows.changed.connect(lambda *_: self.refresh())
        self.name.textChanged.connect(lambda *_: self.refresh())
        self._changed()
        self.setMinimumWidth(LABEL_WIDTH + FIELD_WIDTH + 60)

    # -- inserting ---------------------------------------------------------------------------
    def insert(self, text: str) -> None:
        cursor = self.edit.textCursor()
        cursor.insertText(text)
        if text.endswith("()"):
            cursor.movePosition(QTextCursor.MoveOperation.Left)
        self.edit.setTextCursor(cursor)
        self.edit.setFocus()

    def _fill_columns(self) -> None:
        self.column_menu.clear()
        for name in self.ctx.table.column_names:
            kind = self.ctx.table.kind(name)
            action = self.column_menu.addAction(f"{name}  ({kind})")
            action.triggered.connect(lambda _c=False, n=name: self.insert(f"[{n}]"))

    def _function_menu(self) -> QMenu:
        menu = QMenu(self)
        for group in fx.GROUPS:
            sub = menu.addMenu(group)
            for function in fx.FUNCTIONS.values():
                if function.group != group:
                    continue
                action = sub.addAction(function.signature)
                action.setToolTip(function.help)
                action.setStatusTip(function.help)
                action.triggered.connect(lambda _c=False, f=function: self.insert(f"{f.name}()"))
            sub.setToolTipsVisible(True)
        return menu

    # -- checking ----------------------------------------------------------------------------
    def _changed(self) -> None:
        text = self.edit.toPlainText()
        self.formula, self.error = None, ""
        try:
            self.formula = fx.parse(text, self.ctx.table.column_names)
            self.preview.setText(self._preview_text(self.formula))
        except fx.FormulaError as err:
            self.error = err.message
            self.preview.setText("")
        self.refresh()

    def _preview_text(self, formula: fx.Formula) -> str:
        from chemsuite.plotting.data import shown_oids

        oids = self.rows.oids()
        if formula.aggregate:
            sample = oids
        else:
            shown = shown_oids(self.ctx.table)
            sample = (shown if len(shown) else oids)[:PREVIEW_ROWS]
        sample_n = len(sample)
        if not sample_n:
            return "(no rows)"
        data = {name: column_values(self.ctx.table, name, sample) for name in formula.columns}
        kind, cells = fx.result_values(formula.evaluate(data, sample_n), sample_n)
        if formula.aggregate:
            positions = {int(o): i for i, o in enumerate(sample.tolist())}
            shown = [o for o in shown_oids(self.ctx.table).tolist() if o in positions]
            cells = [cells[positions[o]] for o in shown[:PREVIEW_ROWS]]
        lines = [
            fx.number_text(c) if isinstance(c, float) else ("" if c is None else str(c))
            for c in cells[:PREVIEW_ROWS]
        ]
        lines = [line if line else "(blank)" for line in lines]
        return f"{kind.capitalize()}: " + "  ·  ".join(lines)

    # -- ToolDialog --------------------------------------------------------------------------
    def problem(self) -> str:
        if self.error:
            return self.error
        if self.formula is None:
            return "Type a formula."
        if not self.name.text().strip():
            return "Name the new column."
        return ""

    def verb(self) -> str:
        return f"Calculate {len(self.rows.oids()):,} Rows"

    def run(self):
        formula = self.formula
        assert formula is not None
        ctx = self.ctx
        oids = self.rows.oids()
        table = ctx.table
        data = {name: column_values(table, name, oids, blobs=True) for name in formula.columns}
        name = self.name.text().strip()
        n = len(oids)

        def work(job):
            return evaluate(formula, materialise(data), n, job)

        def done(result) -> None:
            _kind, cells = result
            written = oids[: len(cells)]
            names = ctx.apply_results({name: (written, cells)}, label=f"Calculate {name}")
            partial = "" if len(cells) == n else f" (stopped after {len(cells):,} rows)"
            ctx.show_message(f"Calculated {names.get(name, name)}{partial}: {formula.text}")

        return ctx.run(f"Calculating {name}…", work, done, keep_partial=True)


def open_calculator(ctx: ToolContext) -> CalculatorDialog:
    dialog = CalculatorDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog
