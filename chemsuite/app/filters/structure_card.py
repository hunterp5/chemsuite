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
"""The structure filter card: keep rows matching a substructure, an exact structure or a
similarity threshold.

The search runs as a background job after a short pause in typing, and again when the column's
cells change. Until a new search finishes, the card keeps filtering by the previous matches, so
the table does not flash back to every row while you type.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QWidget,
)

from chemsuite.analysis.structure_search import (
    MODES,
    QueryError,
    SearchResult,
    StructureQuery,
    parse_query,
)
from chemsuite.core.table import Table
from chemsuite.core.table.filters import Filter, OidFilter

from ..jobs import Job
from ..structure_search import QueryInput, StructureSearcher, highlight_for
from .cards import FilterCard

SEARCH_DELAY_MS = 400


class StructureCard(FilterCard):
    kind = "structure"

    def __init__(
        self,
        table: Table,
        column: str,
        searcher: StructureSearcher,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(table, column, parent)
        self.searcher = searcher
        self.job: Job | None = None
        self._hits: frozenset[int] | None = None  # None: nothing to filter by
        self._searched: StructureQuery | None = None  # the query _hits belong to
        self._highlighting = False
        self.highlight_matches = True
        self.query = QueryInput()
        self.mode = QComboBox()
        for key, title in MODES.items():
            self.mode.addItem(title, key)
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.05, 1.0)
        self.threshold.setSingleStep(0.05)
        self.threshold.setValue(0.7)
        self.threshold.setPrefix("≥ ")
        self.stereo = QCheckBox("Stereo")
        self.stereo.setToolTip("Match stereochemistry")
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: palette(placeholder-text)")
        options = QHBoxLayout()
        options.setContentsMargins(0, 0, 0, 0)
        options.addWidget(self.mode, 1)
        options.addWidget(self.threshold)
        options.addWidget(self.stereo)
        holder = QWidget()
        holder.setLayout(options)
        for widget in (self.query, holder, self.status):
            self.body.addWidget(widget)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(SEARCH_DELAY_MS)
        self._timer.timeout.connect(self.search_now)
        self.query.changed.connect(self._timer.start)
        self.mode.currentIndexChanged.connect(lambda _i: self._mode_changed())
        self.threshold.valueChanged.connect(lambda _v: self._timer.start())
        self.stereo.toggled.connect(lambda _on: self._timer.start())
        self._mode_changed(search=False)

    # -- the query ----------------------------------------------------------------------------
    def current_query(self) -> StructureQuery:
        return StructureQuery(
            self.query.text(),
            self.mode.currentData(),
            stereo=self.stereo.isChecked(),
            threshold=self.threshold.value(),
        )

    def set_query(self, query: StructureQuery, *, highlight: bool = True) -> None:
        """Fill in *query* and search at once."""
        self.highlight_matches = highlight
        for widget in (self.query, self.mode, self.threshold, self.stereo):
            widget.blockSignals(True)
        try:
            self.query.set_text(query.text)
            self.mode.setCurrentIndex(self.mode.findData(query.mode))
            self.threshold.setValue(query.threshold)
            self.stereo.setChecked(query.stereo)
        finally:
            for widget in (self.query, self.mode, self.threshold, self.stereo):
                widget.blockSignals(False)
        self._mode_changed(search=False)
        self.search_now()

    def _mode_changed(self, *, search: bool = True) -> None:
        similar = self.mode.currentData() == "similarity"
        self.threshold.setVisible(similar)
        self.stereo.setVisible(not similar)
        if search:
            self._timer.start()

    # -- searching ----------------------------------------------------------------------------
    def search_now(self) -> None:
        self._timer.stop()
        self._cancel()
        query = self.current_query()
        if not query.text.strip():
            self._set_hits(None, None, "")
            return
        try:
            parse_query(query)
        except QueryError as err:
            self._set_hits(None, None, str(err))
            return
        self.status.setText("Searching…")
        holder: list[Job] = []

        def done(oids: np.ndarray, result: SearchResult) -> None:
            if not holder or self.job is not holder[0]:
                return  # superseded by a newer search
            self.job = None
            hits = frozenset(oids[result.hits].tolist())
            self._set_hits(hits, query, f"{len(hits):,} match(es)")

        self.job = self.searcher.start(self.column, query, done)
        holder.append(self.job)

    def _set_hits(
        self, hits: frozenset[int] | None, query: StructureQuery | None, status: str
    ) -> None:
        self._hits = hits
        self._searched = query
        self.status.setText(status)
        self._update_highlight()
        self.changed.emit()

    def _cancel(self) -> None:
        if self.job is not None and self.job.running:
            self.job.cancel()
        self.job = None

    def _update_highlight(self) -> None:
        query = self._searched
        wanted = (
            highlight_for(query)
            if query is not None and self.highlight_matches and self.is_enabled
            else None
        )
        if wanted is not None or self._highlighting:
            self.searcher.highlight(self.column, wanted)
        self._highlighting = wanted is not None

    # -- FilterCard ---------------------------------------------------------------------------
    def body_state(self) -> dict:
        query = self.current_query()
        return {
            "query": query.text,
            "mode": query.mode,
            "stereo": query.stereo,
            "threshold": query.threshold,
            "highlight": self.highlight_matches,
        }

    def restore_body(self, state: dict) -> None:
        mode = state.get("mode")
        threshold = state.get("threshold")
        query = StructureQuery(
            str(state.get("query", "")),
            mode if mode in ("substructure", "exact", "similarity") else "substructure",
            stereo=bool(state.get("stereo", False)),
            threshold=float(threshold) if isinstance(threshold, int | float) else 0.7,
        )
        self.set_query(query, highlight=bool(state.get("highlight", True)))

    @property
    def hits(self) -> frozenset[int] | None:
        return self._hits

    def build(self, negate: bool) -> Filter | None:
        if self._hits is None or self._searched is None:
            return None
        return OidFilter(self._hits, label=self._searched.describe(), negate=negate)

    def refresh(self) -> None:
        """The column's cells or the rows changed: search again."""
        if self.query.text().strip():
            self.search_now()

    def rename(self, new: str) -> None:
        if self._highlighting:
            self.searcher.highlight(self.column, None)
        super().rename(new)
        self._update_highlight()

    def _on_enabled(self, on: bool) -> None:
        super()._on_enabled(on)
        self._update_highlight()

    def dispose(self) -> None:
        self._timer.stop()
        self._cancel()
        if self._highlighting:
            self.searcher.highlight(self.column, None)
            self._highlighting = False
