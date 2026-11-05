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
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
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

from .range_slider import RangeSlider

#: Text columns with at most this many distinct values get a category card by default.
CATEGORY_LIMIT = 200
#: A category list shows at most this many values (the most frequent).
CATEGORY_LIST_LIMIT = 1000
HISTOGRAM_BINS = 40
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

    def __init__(self, table: Table, column: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self.column = column
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setObjectName("FilterCard")

        self.enabled_box = QCheckBox()
        self.enabled_box.setChecked(True)
        self.enabled_box.setToolTip("Apply this filter")
        self.title = QLabel()
        self.title.setStyleSheet("font-weight: 600")
        self.not_button = QToolButton()
        self.not_button.setText("Not")
        self.not_button.setCheckable(True)
        self.not_button.setToolTip("Keep the rows this filter would hide instead")
        self.remove_button = QToolButton()
        self.remove_button.setText("✕")
        self.remove_button.setAutoRaise(True)
        self.remove_button.setToolTip("Remove filter")

        header = QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.addWidget(self.enabled_box)
        header.addWidget(self.title, 1)
        header.addWidget(self.not_button)
        header.addWidget(self.remove_button)
        self.body = QVBoxLayout()
        self.body.setContentsMargins(0, 0, 0, 0)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 6, 8, 8)
        outer.setSpacing(6)
        outer.addLayout(header)
        outer.addLayout(self.body)

        self.enabled_box.toggled.connect(self._on_enabled)
        self.not_button.toggled.connect(lambda _on: self.changed.emit())
        self.remove_button.clicked.connect(self.removeRequested.emit)
        self._update_title()

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

    def rename(self, new: str) -> None:
        self.column = new
        self._update_title()

    def _update_title(self) -> None:
        self.title.setText(self.column)
        self.title.setToolTip(f"{self.column} ({self.kind})")

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
        self.blanks = QCheckBox("Include blanks")
        self.blanks.setChecked(True)
        self.summary = QLabel()
        self.summary.setStyleSheet("color: palette(placeholder-text)")
        bounds = QHBoxLayout()
        bounds.addWidget(self.low_edit)
        bounds.addWidget(QLabel("to"))
        bounds.addWidget(self.high_edit)
        holder = QWidget()
        holder.setLayout(bounds)
        bounds.setContentsMargins(0, 0, 0, 0)
        for widget in (self.slider, holder, self.blanks, self.summary):
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
            self.summary.setText("No numbers in this column")
        else:
            self._lo, self._hi = float(finite.min()), float(finite.max())
            span = (self._lo, self._hi) if self._hi > self._lo else (self._lo - 0.5, self._lo + 0.5)
            counts, _edges = np.histogram(finite, bins=HISTOGRAM_BINS, range=span)
            self.slider.set_histogram(counts)
            self.slider.setEnabled(self.is_enabled)
            self.summary.setText(
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
        self.values.setMaximumHeight(180)
        self.note = QLabel()
        self.note.setStyleSheet("color: palette(placeholder-text)")
        all_button = QPushButton("All")
        none_button = QPushButton("None")
        buttons = QHBoxLayout()
        buttons.setContentsMargins(0, 0, 0, 0)
        buttons.addWidget(all_button)
        buttons.addWidget(none_button)
        buttons.addStretch(1)
        holder = QWidget()
        holder.setLayout(buttons)
        for widget in (self.search, self.values, holder, self.note):
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
        self.note.setText(
            f"Showing the {len(shown):,} most common of {len(counts):,} values"
            if self._hidden_values
            else ""
        )
        self.note.setVisible(bool(self._hidden_values))
        self._apply_search(self.search.text())

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


def card_for(table: Table, column: str, kind: str | None = None) -> FilterCard:
    """A card for *column*: *kind* ("range", "category", "text") or the best fit for the data."""
    if kind is None:
        if table.kind(column) == NUMBER:
            kind = "range"
        elif table.kind(column) == "structure":
            kind = "text"
        else:
            distinct = len(value_counts(table.column(column)))
            kind = "category" if distinct <= CATEGORY_LIMIT else "text"
    cls = {"range": RangeCard, "category": CategoryCard, "text": TextCard}[kind]
    return cls(table, column)
