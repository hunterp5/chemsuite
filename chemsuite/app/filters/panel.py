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
They are saved with the document (:meth:`FilterPanel.state`), as are named filter *sets*: saved
groups of cards you can switch between ("Actives, MW < 450").
"""

from __future__ import annotations

import itertools
import logging
from collections.abc import Callable
from typing import Any

from PySide6.QtCore import QEvent, QSize, Qt, QTimer
from PySide6.QtWidgets import (
    QFrame,
    QInputDialog,
    QMenu,
    QScrollArea,
    QToolBar,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev

from ..glyphs import glyph
from .cards import FILTER_TYPES, FilterCard, card_for, default_source, sources_for

APPLY_DELAY_MS = 120
EXTRAS_KEY = "filters"
_keys = itertools.count(1)
log = logging.getLogger(__name__)


class FilterPanel(QWidget):
    def __init__(self, table: Table, parent: QWidget | None = None, *, searcher=None) -> None:
        super().__init__(parent)
        self.table = table
        self.searcher = searcher  # runs structure cards' searches (None: no structure cards)
        self.sets: dict[str, list[dict]] = {}  # named filter sets: saved card states
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

        # Add, Sets and Remove All are icon buttons in the pane's title bar (``PanelHeader``
        # takes ``toolbar``), like a plot's tools, so the pane itself holds only the cards.
        self.toolbar = QToolBar()
        self.toolbar.setIconSize(QSize(13, 13))
        self.toolbar.setMovable(False)
        self.toolbar.setStyleSheet(
            "QToolBar { background: transparent; border: none; padding: 0px; spacing: 2px; }"
            " QToolButton { border: none; border-radius: 4px; padding: 2px; }"
            " QToolButton:hover { background: palette(midlight); }"
            " QToolButton::menu-indicator { image: none; }"
        )
        self.add_button, self.add_menu = self._tool("plus", "Add Filter: a card for a column")
        self.add_menu.aboutToShow.connect(self._fill_add_menu)
        self.sets_button, self.sets_menu = self._tool(
            "sets", "Filter Sets: save these filters under a name, or switch to a saved set"
        )
        self.sets_menu.aboutToShow.connect(self._fill_sets_menu)
        self.clear_button, _none = self._tool("trash", "Remove All Filters", menu=False)
        self.clear_button.clicked.connect(self.remove_all)
        self._cards_box = QVBoxLayout()
        self._cards_box.setContentsMargins(0, 0, 0, 0)
        self._cards_box.setSpacing(0)  # sections divide themselves with hairlines
        self._cards_box.addStretch(1)
        host = QWidget()
        host.setLayout(self._cards_box)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(host)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)  # striped cards reach the pane's edges
        layout.addWidget(scroll, 1)
        self.setMinimumWidth(240)

        self._unsubscribe = table.subscribe(self._on_table_event)
        self._style_tools()
        self._update_summary()

    # -- cards --------------------------------------------------------------------------------
    @property
    def cards(self) -> list[FilterCard]:
        return list(self._cards)

    def add_filter(self, column: str, kind: str | None = None) -> FilterCard:
        card = card_for(self.table, column, kind, self.searcher)
        self._cards[card] = f"panel:{next(_keys)}"
        self._cards_box.insertWidget(self._cards_box.count() - 1, card)
        self._stripe()
        self._connect(card)
        self._update_summary()
        return card

    def _connect(self, card: FilterCard) -> None:
        card.changed.connect(lambda: self._mark_dirty(card))
        card.removeRequested.connect(lambda: self.remove_filter(card))
        card.sourceRequested.connect(lambda name: self.change_source(card, name))

    def add_filter_of_type(self, kind: str) -> FilterCard | None:
        """A *kind* card on its first suitable column (choose another in the card)."""
        column = default_source(self.table, kind)
        if column is None:
            return None
        return self.add_filter(column, kind)

    def change_source(self, card: FilterCard, column: str) -> FilterCard:
        """Rebuild *card* (same kind, place, on/off and Not) for *column*."""
        if card not in self._cards:
            return card
        new = card_for(self.table, column, card.kind, self.searcher)
        new.set_enabled(card.is_enabled)
        new.not_button.setChecked(card.negated)
        position = self._cards_box.indexOf(card)
        key = self._cards[card]
        # Keep the card's filter key and place: rebuild the order with the new card.
        self._cards = {(new if c is card else c): k for c, k in self._cards.items()}
        self._dirty.discard(card)
        self._stale.discard(card)
        card.dispose()
        self._cards_box.insertWidget(position, new)
        card.setParent(None)
        card.deleteLater()
        self._connect(new)
        if key in self.table.filters:
            self.table.set_filter(key, None)
        self._stripe()
        self._mark_dirty(new)
        self._update_summary()
        return new

    def remove_filter(self, card: FilterCard) -> None:
        key = self._cards.pop(card, None)
        if key is None:
            return
        self._dirty.discard(card)
        self._stale.discard(card)
        card.dispose()
        if key in self.table.filters:
            self.table.set_filter(key, None)
        card.setParent(None)
        card.deleteLater()
        self._stripe()
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
        for card in self._cards:
            card.dispose()
        return super().close()

    # -- saving and sets ----------------------------------------------------------------------
    def card_states(self) -> list[dict]:
        return [card.state() for card in self._cards]

    def load_cards(self, states: list[Any]) -> int:
        """Replace the cards with saved ones; returns how many could not be restored."""
        skipped = 0
        with self.table.batch():
            self.remove_all()
            for state in states:
                column = state.get("column") if isinstance(state, dict) else None
                kind = state.get("kind") if isinstance(state, dict) else None
                if column not in self.table or kind not in (
                    "range",
                    "text",
                    "category",
                    "structure",
                    "reaction",
                ):
                    skipped += 1
                    continue
                if kind in ("structure", "reaction") and (
                    self.searcher is None or self.table.kind(column) != kind
                ):
                    skipped += 1
                    continue
                if kind == "range" and self.table.kind(column) != "number":
                    skipped += 1
                    continue
                try:
                    self.add_filter(column, kind).restore(state)
                except (TypeError, ValueError):
                    log.exception("could not restore filter %r", state)
                    skipped += 1
            self.apply_now()
        return skipped

    def state(self) -> dict[str, Any]:
        """The cards and filter sets, to save with the document ({} when there are none)."""
        out: dict[str, Any] = {}
        if self._cards:
            out["cards"] = self.card_states()
        if self.sets:
            out["sets"] = {name: list(cards) for name, cards in self.sets.items()}
        return out

    def load_state(self, saved: dict[str, Any]) -> None:
        sets = saved.get("sets")
        if isinstance(sets, dict):
            self.sets = {
                str(name): cards for name, cards in sets.items() if isinstance(cards, list)
            }
        cards = saved.get("cards")
        if isinstance(cards, list) and cards:
            skipped = self.load_cards(cards)
            if skipped:
                log.warning("%d saved filter(s) no longer fit the table", skipped)

    def save_set(self, name: str) -> None:
        name = name.strip()
        if not name:
            raise ValueError("A filter set needs a name.")
        self.sets[name] = self.card_states()

    def apply_set(self, name: str) -> int:
        """Switch to the saved set *name*; returns how many of its filters no longer fit."""
        return self.load_cards(self.sets[name])

    def delete_set(self, name: str) -> None:
        self.sets.pop(name, None)

    def ask_save_set(self) -> None:
        name, ok = QInputDialog.getText(
            self, "Save Filter Set", "Name for the current filters (an existing name is replaced):"
        )
        if ok and name.strip():
            self.save_set(name)

    def _fill_sets_menu(self) -> None:
        menu = self.sets_menu
        menu.clear()
        save = menu.addAction("Save Current Filters as Set…")
        save.setEnabled(bool(self._cards))
        save.triggered.connect(lambda _c=False: self.ask_save_set())
        if not self.sets:
            return
        menu.addSeparator()
        for name, cards in self.sets.items():
            action = menu.addAction(f"{name}  ({len(cards)})")
            action.setToolTip(", ".join(str(c.get("column")) for c in cards))
            action.triggered.connect(lambda _c=False, n=name: self.apply_set(n))
        delete = menu.addMenu("Delete Set")
        for name in self.sets:
            delete.addAction(name).triggered.connect(lambda _c=False, n=name: self.delete_set(n))

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

    def summary(self) -> str:
        """ "Showing 3 of 8 rows · 1 active": the Remove All button's tooltip (the status bar
        shows the row counts too)."""
        total = len(self.table)
        shown = len(self.table.view())
        active = sum(1 for key in self._cards.values() if key in self.table.filters)
        return f"Showing {shown:,} of {total:,} rows" + (f" · {active} active" if active else "")

    def _update_summary(self) -> None:
        self.clear_button.setEnabled(bool(self._cards))
        self.clear_button.setToolTip(f"Remove All Filters ({self.summary()})")

    def _tool(self, icon: str, tip: str, *, menu: bool = True):
        """An icon button on the toolbar; with *menu*, it opens a menu: ``(button, menu)``."""
        button = QToolButton()
        button.setToolTip(tip)
        button.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        button.setProperty("glyph", icon)
        popup = None
        if menu:
            popup = QMenu(button)
            popup.setToolTipsVisible(True)
            button.setMenu(popup)
            button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.toolbar.addWidget(button)
        return button, popup

    def _style_tools(self) -> None:
        ink = self.palette().text().color()
        for button in self.toolbar.findChildren(QToolButton):
            name = button.property("glyph")
            if name:
                button.setIcon(glyph(name, ink))

    def changeEvent(self, event) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "toolbar"):
            self._style_tools()

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
        if event.structural:
            for card in self._cards:
                card.update_sources()  # columns added, removed, renamed or converted
        self._update_summary()

    def _drop_cards(self, predicate: Callable[[FilterCard], bool]) -> None:
        for card in [c for c in self._cards if predicate(c)]:
            self._cards.pop(card)  # the table already dropped this card's filter
            self._dirty.discard(card)
            self._stale.discard(card)
            card.dispose()
            card.setParent(None)
            card.deleteLater()
        self._stripe()

    def _stripe(self) -> None:
        """Alternate the cards' backgrounds, like table rows, so neighbours stand apart."""
        for i, card in enumerate(self._cards):
            card.set_alternate(i % 2 == 1)

    def _refresh_stale(self) -> None:
        for card in [c for c in self._stale if c in self._cards]:
            card.refresh()
        self._stale.clear()

    def _fill_add_menu(self) -> None:
        """One entry per kind of filter; the card then picks its column."""
        menu = self.add_menu
        menu.clear()
        for kind, title in FILTER_TYPES.items():
            if kind in ("structure", "reaction") and self.searcher is None:
                continue
            action = menu.addAction(title)
            action.triggered.connect(lambda _c=False, k=kind: self.add_filter_of_type(k))
            if not sources_for(self.table, kind):
                action.setEnabled(False)
                action.setToolTip(f"No column suits a {title} filter")
