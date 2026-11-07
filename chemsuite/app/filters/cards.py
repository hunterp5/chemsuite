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

"""Filter cards: one widget per filter in the filter panel.

A card reads its column through the table, builds an immutable core filter from its controls
(:meth:`FilterCard.current_filter`) and says ``changed`` when that might differ. It never
touches the table's filters itself; the panel applies them.
"""

from __future__ import annotations

import math
import re

import numpy as np
from PySide6.QtCore import QEvent, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.table import NUMBER, Table
from chemsuite.core.table.filters import (
    CategoryFilter,
    Filter,
    RangeFilter,
    TextFilter,
    number_values,
    value_counts,
)
from chemsuite.core.table.numbers import format_general

from ..glyphs import flat_buttons, glyph
from .range_slider import RangeSlider

#: Text columns with at most this many distinct values get a category card by default.
CATEGORY_LIMIT = 200
#: A category list shows at most this many values (the most frequent).
CATEGORY_LIST_LIMIT = 1000
HISTOGRAM_BINS = 40
#: how far a striped card's background moves from the base colour toward the text colour
STRIPE = 0.055
#: a category list grows to this many rows before it scrolls
LIST_ROWS = 8
BLANK_LABEL = "(blank)"

TEXT_MODES = {
    "contains": "Contains",
    "equals": "Equals",
    "starts": "Starts with",
    "ends": "Ends with",
    "regex": "Regular expression",
}


class FilterCard(QFrame):
    """Header (on/off, title, Not, remove) plus a body supplied by subclasses."""

    kind = ""
    changed = Signal()
    removeRequested = Signal()
    #: another column was chosen as the card's source (the panel swaps in a new card)
    sourceRequested = Signal(str)

    def __init__(self, table: Table, column: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self.column = column
        self.setObjectName("FilterCard")
        self.setFrameShape(QFrame.Shape.NoFrame)
        self._alternate = False  # striped like table rows (the panel sets it)
        self._glyphs: dict[QToolButton, str] = {}
        self._restyling = False

        self.enabled_box = QCheckBox()
        self.enabled_box.setChecked(True)
        self.enabled_box.setToolTip("Apply this filter")
        #: the card's source column, chosen here (the card is rebuilt for another column)
        self.title = QComboBox()
        self.title.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.title.setMinimumContentsLength(6)
        font = self.title.font()
        font.setWeight(QFont.Weight.DemiBold)
        self.title.setFont(font)
        self.not_button = self.glyph_button("not", "Not: keep the rows this filter would hide")
        self.not_button.setCheckable(True)
        self.remove_button = self.glyph_button("close", "Remove filter")

        header = QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.setSpacing(4)
        header.addWidget(self.enabled_box)
        header.addSpacing(4)  # the check box reports no gap after its indicator
        header.addWidget(self.title)
        header.addStretch(1)
        header.addWidget(self.not_button)
        header.addWidget(self.remove_button)
        self.body = QVBoxLayout()
        self.body.setContentsMargins(0, 0, 0, 0)
        self.body.setSpacing(6)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 6, 8, 10)
        outer.setSpacing(6)
        outer.addLayout(header)
        outer.addLayout(self.body)

        self.enabled_box.toggled.connect(self._on_enabled)
        self.not_button.toggled.connect(lambda _on: self.changed.emit())
        self.remove_button.clicked.connect(self.removeRequested.emit)
        self.update_sources()
        self.title.activated.connect(self._source_chosen)
        self._restyle()

    # -- look ---------------------------------------------------------------------------------
    def glyph_button(self, name: str, tip: str) -> QToolButton:
        """A flat glyph button, drawn in the palette's colours."""
        button = QToolButton()
        button.setAutoRaise(True)
        button.setIconSize(QSize(14, 14))
        button.setFixedSize(22, 22)
        button.setToolTip(tip)
        button.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._glyphs[button] = name
        return button

    def set_alternate(self, on: bool) -> None:
        """Use the alternate background, like every other table row."""
        if on != self._alternate:
            self._alternate = on
            self._restyle()

    def _restyle(self) -> None:
        palette = self.palette()
        ink = palette.text().color()
        for button, name in self._glyphs.items():
            button.setIcon(glyph(name, ink))
        background = palette.base().color()
        if self._alternate:  # base a shade toward the text, as the native table stripes rows
            ink_rgb, base_rgb = ink.getRgbF()[:3], background.getRgbF()[:3]
            background = QColor.fromRgbF(
                *(b + (i - b) * STRIPE for b, i in zip(base_rgb, ink_rgb, strict=True))
            )
        line = palette.mid().color().name()
        sheet = (
            f"QFrame#FilterCard {{ background: {background.name()}; border: none;"
            f" border-bottom: 1px solid {line}; }} "
            + flat_buttons(palette, "QFrame#FilterCard QToolButton")
        )
        if sheet != self.styleSheet():  # a new sheet changes the palette: no loop
            self.setStyleSheet(sheet)

    def changeEvent(self, event) -> None:
        super().changeEvent(event)
        palette = event.type() == QEvent.Type.PaletteChange
        if palette and hasattr(self, "_glyphs") and not self._restyling:
            self._restyling = True
            try:
                self._restyle()
            finally:
                self._restyling = False

    # -- state --------------------------------------------------------------------------------
    @property
    def is_enabled(self) -> bool:
        return self.enabled_box.isChecked()

    def set_enabled(self, on: bool) -> None:
        self.enabled_box.setChecked(on)

    @property
    def negated(self) -> bool:
        return self.not_button.isChecked()

    def current_filter(self) -> Filter | None:
        """The filter to apply, or None when the card is off or does not constrain anything."""
        return self.build(self.negated) if self.is_enabled else None

    def build(self, negate: bool) -> Filter | None:
        raise NotImplementedError

    def refresh(self) -> None:
        """The column's data changed: update summaries, keeping the user's choices."""

    def dispose(self) -> None:
        """The card is being removed: stop background work and undo side effects."""

    # -- saving -------------------------------------------------------------------------------
    def state(self) -> dict:
        """The card's settings as JSON-ready data (see :meth:`restore`)."""
        return {
            "column": self.column,
            "kind": self.kind,
            "enabled": self.is_enabled,
            "negate": self.negated,
            **self.body_state(),
        }

    def restore(self, state: dict) -> None:
        """Apply settings saved by :meth:`state` (unknown or bad values are ignored)."""
        self.not_button.setChecked(bool(state.get("negate", False)))
        self.restore_body(state)
        self.set_enabled(bool(state.get("enabled", True)))
        self.changed.emit()

    def body_state(self) -> dict:
        return {}

    def restore_body(self, state: dict) -> None:
        pass

    def rename(self, new: str) -> None:
        self.column = new
        self.update_sources()

    def update_sources(self) -> None:
        """List the columns this kind of card can filter (the table's columns changed)."""
        names = sources_for(self.table, self.kind)
        if self.column not in names:
            names = [self.column, *names]
        if names != [self.title.itemText(i) for i in range(self.title.count())]:
            self.title.blockSignals(True)
            self.title.clear()
            self.title.addItems(names)
            self.title.blockSignals(False)
        self.title.setCurrentText(self.column)
        self._update_title()

    def _source_chosen(self, index: int) -> None:
        name = self.title.itemText(index)
        if name and name != self.column:
            self.sourceRequested.emit(name)

    #: what the card covers ("7 numbers from 5 to 9"), shown in the title's tooltip: cards
    #: carry no text of their own under the controls
    details = ""

    def set_details(self, text: str) -> None:
        self.details = text
        self._update_title()

    def _update_title(self) -> None:
        tip = f"{FILTER_TYPES.get(self.kind, self.kind)} filter on {self.column}"
        self.title.setToolTip(f"{tip}\n{self.details}" if self.details else tip)

    def _on_enabled(self, on: bool) -> None:
        for i in range(self.body.count()):
            item = self.body.itemAt(i)
            widget = item.widget() if item is not None else None
            if widget is not None:
                widget.setEnabled(on)
        self.changed.emit()


# -- numbers ----------------------------------------------------------------------------------
class RangeCard(FilterCard):
    """Keep rows whose number lies between two bounds; a histogram shows the distribution."""

    kind = "range"

    def __init__(self, table: Table, column: str, parent: QWidget | None = None) -> None:
        super().__init__(table, column, parent)
        self._lo = 0.0  # data bounds
        self._hi = 0.0
        self._low: float | None = None  # chosen bounds (None = open)
        self._high: float | None = None
        self.slider = RangeSlider()
        self.low_edit = QLineEdit()
        self.high_edit = QLineEdit()
        for edit, text in ((self.low_edit, "minimum"), (self.high_edit, "maximum")):
            edit.setPlaceholderText(text)
            edit.setAlignment(Qt.AlignmentFlag.AlignRight)
            edit.setMinimumWidth(48)  # shrink with the pane rather than clip
        self.blanks = QCheckBox("Include blanks")
        self.blanks.setChecked(True)
        bounds = QHBoxLayout()
        bounds.addWidget(self.low_edit)
        bounds.addWidget(QLabel("to"))
        bounds.addWidget(self.high_edit)
        holder = QWidget()
        holder.setLayout(bounds)
        bounds.setContentsMargins(0, 0, 0, 0)
        for widget in (self.slider, holder, self.blanks):
            self.body.addWidget(widget)

        self.slider.rangeChanged.connect(self._from_slider)
        self.low_edit.editingFinished.connect(self._from_edits)
        self.high_edit.editingFinished.connect(self._from_edits)
        self.blanks.toggled.connect(lambda _on: self.changed.emit())
        self.refresh()

    def refresh(self) -> None:
        values = number_values(self.table.column(self.column))
        finite = values[np.isfinite(values)]
        blanks = len(values) - len(finite)
        if not len(finite):
            self._lo = self._hi = 0.0
            self.slider.set_histogram(None)
            self.slider.setEnabled(False)
            self.set_details("No numbers in this column")
        else:
            self._lo, self._hi = float(finite.min()), float(finite.max())
            span = (self._lo, self._hi) if self._hi > self._lo else (self._lo - 0.5, self._lo + 0.5)
            counts, _edges = np.histogram(finite, bins=HISTOGRAM_BINS, range=span)
            self.slider.set_histogram(counts)
            self.slider.setEnabled(self.is_enabled)
            self.set_details(
                f"{len(finite):,} numbers from {format_general(self._lo)} to "
                f"{format_general(self._hi)}" + (f" · {blanks:,} blank" if blanks else "")
            )
        # keep chosen bounds, but drop ones the data no longer reaches past
        if self._low is not None and self._low <= self._lo:
            self._low = None
        if self._high is not None and self._high >= self._hi:
            self._high = None
        self._show_bounds()

    def set_bounds(self, low: float | None, high: float | None) -> None:
        if low is not None and high is not None and low > high:
            low, high = high, low
        self._low = None if low is None or low <= self._lo else float(low)
        self._high = None if high is None or high >= self._hi else float(high)
        self._show_bounds()
        self.changed.emit()

    def bounds(self) -> tuple[float | None, float | None]:
        return self._low, self._high

    def body_state(self) -> dict:
        return {"low": self._low, "high": self._high, "blanks": self.blanks.isChecked()}

    def restore_body(self, state: dict) -> None:
        def number(value) -> float | None:
            return float(value) if isinstance(value, int | float) else None

        self.blanks.setChecked(bool(state.get("blanks", True)))
        self.set_bounds(number(state.get("low")), number(state.get("high")))

    def build(self, negate: bool) -> Filter | None:
        low, high = self._low, self._high
        include_blank = self.blanks.isChecked()
        if low is None and high is None and include_blank:
            return None
        return RangeFilter(self.column, low, high, include_blank=include_blank, negate=negate)

    # -- syncing controls -----------------------------------------------------------------------
    def _value_of(self, position: int) -> float:
        return self._lo + (self._hi - self._lo) * position / self.slider.steps

    def _position_of(self, value: float) -> int:
        if self._hi <= self._lo:
            return 0
        return round((value - self._lo) / (self._hi - self._lo) * self.slider.steps)

    def _from_slider(self, low: int, high: int) -> None:
        steps = self.slider.steps
        self._low = None if low <= 0 else self._value_of(low)
        self._high = None if high >= steps else self._value_of(high)
        self._show_edits()
        self.changed.emit()

    def _from_edits(self) -> None:
        def parse(edit: QLineEdit) -> float | None:
            text = edit.text().strip()
            try:
                value = float(text)
            except ValueError:
                return None
            return value if math.isfinite(value) else None

        self.set_bounds(parse(self.low_edit), parse(self.high_edit))

    def _show_bounds(self) -> None:
        low = 0 if self._low is None else self._position_of(self._low)
        high = self.slider.steps if self._high is None else self._position_of(self._high)
        self.slider.set_values(low, high, emit=False)
        self._show_edits()

    def _show_edits(self) -> None:
        self.low_edit.setText("" if self._low is None else format_general(self._low))
        self.high_edit.setText("" if self._high is None else format_general(self._high))


# -- text -------------------------------------------------------------------------------------
class TextCard(FilterCard):
    kind = "text"

    def __init__(self, table: Table, column: str, parent: QWidget | None = None) -> None:
        super().__init__(table, column, parent)
        self.query = QLineEdit()
        self.query.setPlaceholderText("Text to match")
        self.query.setClearButtonEnabled(True)
        self.mode = QComboBox()
        for key, title in TEXT_MODES.items():
            self.mode.addItem(title, key)
        self.case = QCheckBox("Match case")
        options = QHBoxLayout()
        options.setContentsMargins(0, 0, 0, 0)
        options.addWidget(self.mode, 1)
        options.addWidget(self.case)
        holder = QWidget()
        holder.setLayout(options)
        self.body.addWidget(self.query)
        self.body.addWidget(holder)
        self.query.textChanged.connect(lambda _t: self.changed.emit())
        self.mode.currentIndexChanged.connect(lambda _i: self.changed.emit())
        self.case.toggled.connect(lambda _on: self.changed.emit())

    def body_state(self) -> dict:
        return {
            "query": self.query.text(),
            "mode": self.mode.currentData(),
            "case": self.case.isChecked(),
        }

    def restore_body(self, state: dict) -> None:
        self.query.setText(str(state.get("query", "")))
        index = self.mode.findData(state.get("mode"))
        if index >= 0:
            self.mode.setCurrentIndex(index)
        self.case.setChecked(bool(state.get("case", False)))

    def build(self, negate: bool) -> Filter | None:
        text = self.query.text()
        mode = self.mode.currentData()
        valid = True
        if mode == "regex":
            try:
                re.compile(text)
            except re.error as err:
                valid = False
                self.query.setToolTip(f"Invalid regular expression: {err}")
        self.query.setStyleSheet("" if valid else "border: 1px solid #d33")
        if valid:
            self.query.setToolTip("")
        if not text or not valid:
            return None
        return TextFilter(self.column, text, mode, self.case.isChecked(), negate)


# -- categories -------------------------------------------------------------------------------
class CategoryCard(FilterCard):
    """Tick the values to keep; counts show how many rows have each value."""

    kind = "category"

    def __init__(self, table: Table, column: str, parent: QWidget | None = None) -> None:
        super().__init__(table, column, parent)
        self._unchecked: set[str] = set()  # remembered across refreshes; new values start ticked
        self._hidden_values = 0
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find a value")
        self.search.setClearButtonEnabled(True)
        self.values = QListWidget()
        self.values.setUniformItemSizes(True)
        self.all_button = self.glyph_button("tick_all", "Tick every value")
        self.none_button = self.glyph_button("tick_none", "Untick every value")
        all_button, none_button = self.all_button, self.none_button
        find = QHBoxLayout()
        find.setContentsMargins(0, 0, 0, 0)
        find.setSpacing(2)
        find.addWidget(self.search, 1)
        find.addSpacing(4)
        find.addWidget(all_button)
        find.addWidget(none_button)
        holder = QWidget()
        holder.setLayout(find)
        for widget in (holder, self.values):
            self.body.addWidget(widget)

        self.search.textChanged.connect(self._apply_search)
        self.values.itemChanged.connect(self._on_item)
        all_button.clicked.connect(lambda: self.set_checked(None))
        none_button.clicked.connect(lambda: self.set_checked(set()))
        self.refresh()

    def refresh(self) -> None:
        counts = value_counts(self.table.column(self.column))
        shown = counts[:CATEGORY_LIST_LIMIT]
        self._hidden_values = len(counts) - len(shown)
        self.values.blockSignals(True)
        try:
            self.values.clear()
            for value, count in shown:
                item = QListWidgetItem(f"{value or BLANK_LABEL}    {count:,}")
                item.setData(Qt.ItemDataRole.UserRole, value)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                checked = value not in self._unchecked
                item.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)
                self.values.addItem(item)
        finally:
            self.values.blockSignals(False)
        self._fit_list()
        self.set_details(
            f"Listing the {len(shown):,} most common of {len(counts):,} values"
            if self._hidden_values
            else f"{len(counts):,} values"
        )
        self._apply_search(self.search.text())

    def _fit_list(self) -> None:
        """As tall as its values, up to LIST_ROWS rows (then it scrolls)."""
        rows = min(max(self.values.count(), 1), LIST_ROWS)
        row = self.values.sizeHintForRow(0) if self.values.count() else 0
        row = row if row > 0 else self.fontMetrics().height() + 4
        self.values.setFixedHeight(rows * row + 2 * self.values.frameWidth() + 2)

    def listed_values(self) -> list[str]:
        return [
            self.values.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.values.count())
        ]

    def checked_values(self) -> set[str]:
        return {v for v in self.listed_values() if v not in self._unchecked}

    def set_checked(self, values: set[str] | None) -> None:
        """Tick exactly *values* (None ticks everything)."""
        listed = self.listed_values()
        self._unchecked = set() if values is None else {v for v in listed if v not in values}
        self.values.blockSignals(True)
        try:
            for i, value in enumerate(listed):
                state = (
                    Qt.CheckState.Unchecked if value in self._unchecked else Qt.CheckState.Checked
                )
                self.values.item(i).setCheckState(state)
        finally:
            self.values.blockSignals(False)
        self.changed.emit()

    def body_state(self) -> dict:
        return {"unticked": sorted(self._unchecked)}

    def restore_body(self, state: dict) -> None:
        unticked = state.get("unticked")
        if isinstance(unticked, list):
            self._unchecked = {str(v) for v in unticked}
            self.refresh()

    def build(self, negate: bool) -> Filter | None:
        # Expressed as "hide the unticked values", so values too rare to be listed stay visible.
        unticked = self._unchecked & set(self.listed_values())
        if not unticked:
            return None
        return CategoryFilter(self.column, frozenset(unticked), negate=not negate)

    def _on_item(self, item: QListWidgetItem) -> None:
        value = item.data(Qt.ItemDataRole.UserRole)
        if item.checkState() == Qt.CheckState.Checked:
            self._unchecked.discard(value)
        else:
            self._unchecked.add(value)
        self.changed.emit()

    def _apply_search(self, text: str) -> None:
        needle = text.casefold()
        for i in range(self.values.count()):
            item = self.values.item(i)
            value = item.data(Qt.ItemDataRole.UserRole) or BLANK_LABEL
            item.setHidden(bool(needle) and needle not in value.casefold())


#: the kinds of filter card, as the panel's Add menu offers them
FILTER_TYPES = {
    "range": "Number Range",
    "category": "Category",
    "text": "Text Match",
    "structure": "Structure",
    "reaction": "Reaction",
}


def sources_for(table: Table, kind: str) -> list[str]:
    """The columns a *kind* card can filter: numbers for a range, structures for a structure
    card, reactions for a reaction card, text and numbers for categories, any column for a
    text match."""

    def fits(name: str) -> bool:
        column_kind = table.kind(name)
        if kind == "range":
            return column_kind == NUMBER
        if kind in ("structure", "reaction"):
            return column_kind == kind
        if kind == "category":
            return column_kind not in ("structure", "reaction")
        return True

    return [name for name in table.column_names if fits(name)]


def default_source(table: Table, kind: str) -> str | None:
    """The column a new *kind* card starts on: for categories, the first column with a short
    list of values (not an ID column); for text, the first text column; else the first that
    fits."""
    names = sources_for(table, kind)
    if kind == "category":
        for name in names:
            distinct = len(value_counts(table.column(name)))
            if 2 <= distinct <= CATEGORY_LIMIT and distinct < len(table):  # not an ID column
                return name
    if kind == "text":
        text = [n for n in names if table.kind(n) == "text"]
        if text:
            return text[0]
    return names[0] if names else None


def card_for(table: Table, column: str, kind: str | None = None, searcher=None) -> FilterCard:
    """A card for *column*: *kind* ("range", "category", "text", "structure", "reaction") or
    the best fit.

    Structure and reaction cards need a
    :class:`~chemsuite.app.structure_search.StructureSearcher`; without one, those columns get
    a text card.
    """
    if kind is None:
        if table.kind(column) == NUMBER:
            kind = "range"
        elif table.kind(column) in ("structure", "reaction"):
            kind = table.kind(column) if searcher is not None else "text"
        else:
            distinct = len(value_counts(table.column(column)))
            kind = "category" if distinct <= CATEGORY_LIMIT else "text"
    if kind == "structure":
        if searcher is None:
            raise ValueError("a structure card needs a searcher")
        from .structure_card import StructureCard

        return StructureCard(table, column, searcher)
    if kind == "reaction":
        if searcher is None:
            raise ValueError("a reaction card needs a searcher")
        from .reaction_card import ReactionCard

        return ReactionCard(table, column, searcher)
    cls = {"range": RangeCard, "category": CategoryCard, "text": TextCard}[kind]
    return cls(table, column)
