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

"""Chemistry ▸ Split Reaction."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox

from chemsuite.analysis.reactions import POOL_MIN_ROWS, SPLIT_ROLES, split_cells
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.core.mols import parse_structure
from chemsuite.core.table.columns import StructureColumn

from .common import cells

#: role -> (title of one column per role, title of one column per component)
TITLES = {
    "reactant": ("Reactants", "Reactant"),
    "agent": ("Agents", "Agent"),
    "product": ("Products", "Product"),
}


def column_names(
    source: str, split: list[dict[str, list[str]] | None], each: bool
) -> dict[tuple[str, int], str]:
    """``(role, position) -> column name``; position -1 is the whole role. Roles nothing has
    (often agents) get no column."""
    out: dict[tuple[str, int], str] = {}
    for role in SPLIT_ROLES:
        most = max((len(row[role]) for row in split if row is not None), default=0)
        if not most:
            continue
        if each:
            for i in range(most):
                out[(role, i)] = f"{source} {TITLES[role][1]} {i + 1}"
        else:
            out[(role, -1)] = f"{source} {TITLES[role][0]}"
    return out


class SplitReactionDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Split Reaction", guide="split-reactions")
        self.ctx = ctx

        inputs = self.section("Input")
        self.column = QComboBox()
        self.column.addItems(ctx.reaction_columns())
        default = ctx.default_reaction_column()
        if default:
            self.column.setCurrentText(default)
        inputs.row("Reactions", self.column)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        output = self.section("Output")
        self.layout_choice = output.row(
            "Columns",
            Segmented((("role", "One per role"), ("each", "One per component")), "role"),
            help="One per role: Reactants, Agents and Products, each holding every component "
            "of that role (as reaction SMILES writes them, joined by “.”). One per component: "
            "Reactant 1, Reactant 2 … in the order the reaction lists them.",
        )
        self.keep_maps = QCheckBox("Keep atom map numbers")
        output.row(
            "",
            self.keep_maps,
            help="Off: map numbers are removed, so the same molecule compares equal across "
            "reactions (structure search, duplicates, descriptors).",
        )

        self.column.currentTextChanged.connect(lambda *_: self.refresh())
        for signal in (self.rows.changed, self.layout_choice.changed):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.column.currentText():
            return "There is no reaction column to split."
        if not len(self.rows.oids()):
            return "There are no rows to split."
        return ""

    def verb(self) -> str:
        return f"Split {len(self.rows.oids()):,} Reactions"

    def run(self):
        ctx = self.ctx
        source = self.column.currentText()
        oids = self.rows.oids()
        data = cells(ctx.table, source, oids)
        keep_maps = self.keep_maps.isChecked()
        each = self.layout_choice.value() == "each"
        use_pool = len(data) >= POOL_MIN_ROWS

        def work(job):
            return split_cells(
                data,
                keep_maps=keep_maps,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(split: list) -> None:
            finish(ctx, oids[: len(split)], split, source, each)

        return ctx.run("Splitting reactions…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    split: list[dict[str, list[str]] | None],
    source: str,
    each: bool,
) -> dict[str, str]:
    columns: dict[str, tuple[np.ndarray, StructureColumn]] = {}
    for (role, position), name in column_names(source, split, each).items():
        texts = []
        for row in split:
            parts = [] if row is None else row[role]
            if position < 0:
                texts.append(".".join(p for p in parts if p))
            else:
                texts.append(parts[position] if position < len(parts) else "")
        blobs = [parse_structure(t)[0] if t else b"" for t in texts]
        columns[name] = (oids, StructureColumn.from_parsed(blobs, texts))
    if not columns:
        ctx.show_message("No reaction in these rows could be read", 8000)
        return {}
    written = ctx.apply_results(columns, label="Split Reaction")
    failed = sum(1 for row in split if row is None)
    message = f"Split {len(oids) - failed:,} reaction(s) into {len(written)} column(s)"
    if failed:
        message += f"; {failed:,} blank or unreadable"
    ctx.show_message(message, 8000)
    return written
