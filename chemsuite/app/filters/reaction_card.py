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

"""A filter card for a reaction column: keep the rows whose reaction matches a query.

The query is a reaction (each drawn reactant, agent and product must be a substructure of one
of the reaction's own) or a structure in one role (a reactant, product, agent or any
component); see :mod:`chemsuite.analysis.reactions`. As with structure cards, the search runs
in the background after a pause in typing, and the card keeps filtering by the previous matches
until a new search finishes.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QComboBox, QLabel, QWidget

from chemsuite.analysis.reactions import ROLES, QueryError, ReactionQuery, parse_reaction_query
from chemsuite.core.table import Table
from chemsuite.core.table.filters import Filter, OidFilter

from ..jobs import Job
from ..structure_search import StructureSearcher
from .cards import FilterCard

SEARCH_DELAY_MS = 400


class ReactionCard(FilterCard):
    kind = "reaction"

    def __init__(
        self,
        table: Table,
        column: str,
        searcher: StructureSearcher,
        parent: QWidget | None = None,
    ) -> None:
        from chemsuite.tools.reactions.common import ReactionText

        super().__init__(table, column, parent)
        self.searcher = searcher
        self.job: Job | None = None
        self._hits: frozenset[int] | None = None
        self._searched: ReactionQuery | None = None
        self.role = QComboBox()
        for key, title in ROLES.items():
            self.role.addItem(title, key)
        self.query = ReactionText(dialog_width=False)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: palette(placeholder-text)")
        for widget in (self.role, self.query, self.status):
            self.body.addWidget(widget)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(SEARCH_DELAY_MS)
        self._timer.timeout.connect(self.search_now)
        self.query.changed.connect(self._timer.start)
        self.role.currentIndexChanged.connect(lambda _i: self._role_changed())
        self._role_changed(search=False)

    def _role_changed(self, *, search: bool = True) -> None:
        self.query.set_reactions(self.role.currentData() == "reaction")
        if search:
            self._timer.start()

    def current_query(self) -> ReactionQuery:
        return ReactionQuery(self.query.text(), self.role.currentData())

    def set_query(self, query: ReactionQuery) -> None:
        """Fill in *query* and search at once."""
        for widget in (self.role, self.query):
            widget.blockSignals(True)
        try:
            index = self.role.findData(query.role)
            self.role.setCurrentIndex(max(index, 0))
            self.query.set_text(query.text)
        finally:
            for widget in (self.role, self.query):
                widget.blockSignals(False)
        self._role_changed(search=False)
        self.search_now()

    # -- searching ----------------------------------------------------------------------------
    def search_now(self) -> None:
        self._timer.stop()
        self._cancel()
        query = self.current_query()
        if not query.text.strip():
            self._set_hits(None, None, "")
            return
        try:
            parse_reaction_query(query)
        except QueryError as err:
            self._set_hits(None, None, str(err))
            return
        self.status.setText("Searching…")
        holder: list[Job] = []

        def done(oids: np.ndarray, hits: np.ndarray) -> None:
            if not holder or self.job is not holder[0]:
                return  # superseded by a newer search
            self.job = None
            found = frozenset(oids[hits].tolist())
            self._set_hits(found, query, f"{len(found):,} match(es)")

        self.job = self.searcher.start_reactions(self.column, query, done)
        holder.append(self.job)

    def _set_hits(
        self, hits: frozenset[int] | None, query: ReactionQuery | None, status: str
    ) -> None:
        self._hits = hits
        self._searched = query
        self.status.setText(status)
        self.changed.emit()

    def _cancel(self) -> None:
        if self.job is not None and self.job.running:
            self.job.cancel()
        self.job = None

    # -- FilterCard ---------------------------------------------------------------------------
    def body_state(self) -> dict:
        query = self.current_query()
        return {"query": query.text, "role": query.role}

    def restore_body(self, state: dict) -> None:
        role = state.get("role")
        self.set_query(
            ReactionQuery(str(state.get("query", "")), role if role in ROLES else "reaction")
        )

    @property
    def hits(self) -> frozenset[int] | None:
        return self._hits

    def build(self, negate: bool) -> Filter | None:
        if self._hits is None or self._searched is None:
            return None
        return OidFilter(self._hits, label=self._searched.describe(), negate=negate)

    def refresh(self) -> None:
        """The column's cells or the rows changed: search again."""
        if self.query.text():
            self.search_now()

    def dispose(self) -> None:
        self._timer.stop()
        self._cancel()
