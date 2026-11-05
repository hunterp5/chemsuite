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

"""The filter panel: a column of filter cards that together decide which rows are shown.

Each card owns one key in the table's filters (``panel:<n>``); the table ANDs them. Changes are
applied after a short pause so typing or dragging a slider does not refilter on every event.
Filters are view state, not edits: they are not undoable and do not mark the document modified.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QMenu,
    QScrollArea,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev

from .cards import FilterCard, card_for

APPLY_DELAY_MS = 120
_keys = itertools.count(1)


class FilterPanel(QWidget):
    def __init__(self, table: Table, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self._cards: dict[FilterCard, str] = {}  # card -> filter key
        self._dirty: set[FilterCard] = set()
        self._stale: set[FilterCard] = set()
        self._apply_timer = QTimer(self)
        self._apply_timer.setSingleShot(True)
        self._apply_timer.setInterval(APPLY_DELAY_MS)
        self._apply_timer.timeout.connect(self.apply_now)
        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.setInterval(250)
        self._refresh_timer.timeout.connect(self._refresh_stale)

        self.add_button = QToolButton()
        self.add_button.setText("Add Filter")
        self.add_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.add_menu = QMenu(self.add_button)
        self.add_menu.aboutToShow.connect(self._fill_add_menu)
        self.add_button.setMenu(self.add_menu)
        self.clear_button = QToolButton()
        self.clear_button.setText("Remove All")
        self.clear_button.clicked.connect(self.remove_all)
        self.summary = QLabel()
        self.summary.setStyleSheet("color: palette(placeholder-text)")

        bar = QHBoxLayout()
        bar.addWidget(self.add_button)
        bar.addWidget(self.clear_button)
        bar.addStretch(1)
        self._cards_box = QVBoxLayout()
        self._cards_box.setSpacing(8)
        self._cards_box.addStretch(1)
        host = QWidget()
        host.setLayout(self._cards_box)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(host)
        self.empty_hint = QLabel("No filters. Use Add Filter, or right-click a column header.")
        self.empty_hint.setWordWrap(True)
        self.empty_hint.setStyleSheet("color: palette(placeholder-text)")

        layout = QVBoxLayout(self)
        layout.addLayout(bar)
        layout.addWidget(self.summary)
        layout.addWidget(self.empty_hint)
        layout.addWidget(scroll, 1)
        self.setMinimumWidth(240)

        self._unsubscribe = table.subscribe(self._on_table_event)
        self._update_summary()

    # -- cards --------------------------------------------------------------------------------
    @property
    def cards(self) -> list[FilterCard]:
        return list(self._cards)

    def add_filter(self, column: str, kind: str | None = None) -> FilterCard:
        card = card_for(self.table, column, kind)
        self._cards[card] = f"panel:{next(_keys)}"
        self._cards_box.insertWidget(self._cards_box.count() - 1, card)
        card.changed.connect(lambda: self._mark_dirty(card))
        card.removeRequested.connect(lambda: self.remove_filter(card))
        self._update_summary()
        return card

    def remove_filter(self, card: FilterCard) -> None:
        key = self._cards.pop(card, None)
        if key is None:
            return
        self._dirty.discard(card)
        self._stale.discard(card)
        if key in self.table.filters:
            self.table.set_filter(key, None)
        card.setParent(None)
        card.deleteLater()
        self._update_summary()

    def remove_all(self) -> None:
        with self.table.batch():
            for card in list(self._cards):
                self.remove_filter(card)

    def set_all_enabled(self, on: bool) -> None:
        for card in self._cards:
            card.set_enabled(on)
        self.apply_now()

    def close(self) -> bool:
        self._unsubscribe()
        return super().close()

    # -- applying -----------------------------------------------------------------------------
    def _mark_dirty(self, card: FilterCard) -> None:
        self._dirty.add(card)
        self._apply_timer.start()

    def apply_now(self) -> None:
        """Push pending card changes to the table (one view update for all of them)."""
        self._apply_timer.stop()
        dirty = [card for card in self._dirty if card in self._cards]
        self._dirty.clear()
        with self.table.batch():
            for card in dirty:
                key = self._cards[card]
                flt = card.current_filter()
                if flt is not None or key in self.table.filters:
                    self.table.set_filter(key, flt)
        self._update_summary()

    def _update_summary(self) -> None:
        total = len(self.table)
        shown = len(self.table.view())
        active = sum(1 for key in self._cards.values() if key in self.table.filters)
        self.summary.setText(
            f"Showing {shown:,} of {total:,} rows" + (f" · {active} active" if active else "")
        )
        self.empty_hint.setVisible(not self._cards)

    # -- following the table ------------------------------------------------------------------
    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        for change in event.changes:
            if isinstance(change, ev.ColumnRemoved):
                self._drop_cards(lambda card, name=change.name: card.column == name)
            elif isinstance(change, ev.ColumnKindChanged):
                # The table dropped the old filter; a card built for the old kind cannot follow.
                self._drop_cards(lambda card, name=change.name: card.column == name)
            elif isinstance(change, ev.ColumnRenamed):
                for card in self._cards:
                    if card.column == change.old:
                        card.rename(change.new)
            elif isinstance(change, ev.CellsChanged | ev.RowsInserted | ev.RowsRemoved):
                touched = getattr(change, "column", None)
                for card in self._cards:
                    if touched is None or card.column == touched:
                        self._stale.add(card)
                self._refresh_timer.start()
        self._update_summary()

    def _drop_cards(self, predicate: Callable[[FilterCard], bool]) -> None:
        for card in [c for c in self._cards if predicate(c)]:
            self._cards.pop(card)  # the table already dropped this card's filter
            self._dirty.discard(card)
            self._stale.discard(card)
            card.setParent(None)
            card.deleteLater()

    def _refresh_stale(self) -> None:
        for card in [c for c in self._stale if c in self._cards]:
            card.refresh()
        self._stale.clear()

    def _fill_add_menu(self) -> None:
        menu = self.add_menu
        menu.clear()
        names = self.table.column_names
        if not names:
            menu.addAction("No columns").setEnabled(False)
            return
        for name in names:
            action = menu.addAction(name)
            action.triggered.connect(lambda _c=False, n=name: self.add_filter(n))
        menu.addSeparator()
        text = menu.addMenu("Text Match On")
        for name in names:
            action = text.addAction(name)
            action.triggered.connect(lambda _c=False, n=name: self.add_filter(n, "text"))
