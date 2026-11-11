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

"""Data ▸ Search ▸ Reaction Search."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QComboBox

from chemsuite.analysis.reactions import (
    POOL_MIN_ROWS,
    ROLES,
    QueryError,
    ReactionQuery,
    parse_reaction_query,
    search_reactions,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog

from .common import ReactionText, cells, current_reaction

COLUMN = "Reaction Search"

ROLE_HELP = (
    "Reaction contains: each reactant, agent and product you draw is a substructure of one "
    "of the reaction's own (draw agents above the arrow). The other choices find a structure "
    "in one part of the reaction."
)


class ReactionSearchDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Reaction Search", guide="reaction-search")
        self.ctx = ctx

        inputs = self.section("Input")
        self.column = QComboBox()
        self.column.addItems(ctx.reaction_columns())
        default = ctx.default_reaction_column()
        if default:
            self.column.setCurrentText(default)
        inputs.row("Reactions", self.column)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        query = self.section("Query")
        self.role = QComboBox()
        for key, title in ROLES.items():
            self.role.addItem(title, key)
        query.row("Find", self.role, help=ROLE_HELP)
        self.query = query.row(
            "Query",
            ReactionText(current=current_reaction(ctx, self.column.currentText)),
            help="Reaction SMILES or SMARTS (C(=O)O.N>>C(=O)N), or Draw…",
        )

        output = self.section("Output")
        self.output = output.row(
            "Matches",
            Segmented((("select", "Select them"), ("column", "Add a column")), "select"),
            help=f"Select the matching rows, or mark them Yes / No in a “{COLUMN}” column.",
        )

        self.role.currentIndexChanged.connect(lambda _i: self._role_changed())
        self.query.changed.connect(self.refresh)
        for signal in (self.rows.changed, self.output.changed):
            signal.connect(lambda *_: self.refresh())
        self.column.currentTextChanged.connect(lambda *_: self.refresh())
        self._role_changed()

    def _role_changed(self) -> None:
        self.query.set_reactions(self.role.currentData() == "reaction")
        self.refresh()

    def current_query(self) -> ReactionQuery:
        return ReactionQuery(self.query.text(), self.role.currentData())

    def problem(self) -> str:
        if not self.column.currentText():
            return "There is no reaction column to search."
        if not len(self.rows.oids()):
            return "There are no rows to search."
        try:
            parse_reaction_query(self.current_query())
        except QueryError as err:
            return str(err)
        return ""

    def verb(self) -> str:
        return f"Search {len(self.rows.oids()):,} Rows"

    def run(self):
        ctx = self.ctx
        oids = self.rows.oids()
        query = self.current_query()
        data = cells(ctx.table, self.column.currentText(), oids)
        select = self.output.value() == "select"
        use_pool = len(data) >= POOL_MIN_ROWS

        def work(job):
            return search_reactions(
                query,
                data,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(hits: np.ndarray) -> None:
            finish(ctx, oids[: len(hits)], hits, query, select, len(oids))

        return ctx.run("Searching reactions…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    hits: np.ndarray,
    query: ReactionQuery,
    select: bool,
    total: int,
) -> None:
    """Select or mark the matches among *oids* (the rows searched before any cancel)."""
    found = oids[hits]
    partial = f" (stopped after {len(oids):,} of {total:,} rows)" if len(oids) < total else ""
    if select:
        ctx.selection.set(found.tolist(), "reaction-search")
    else:
        values = ["Yes" if hit else "No" for hit in hits.tolist()]
        ctx.apply_results({COLUMN: (oids, values)}, label="Reaction Search")
    ctx.show_message(f"{len(found):,} reaction(s) match {query.describe()}{partial}", 8000)
