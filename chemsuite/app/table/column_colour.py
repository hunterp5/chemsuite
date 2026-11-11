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

"""Colour Column…: cell backgrounds from a column's values (header menu).

Numbers can run along a two- or three-colour gradient; any column can give each distinct value
its own colour. Range boxes left empty follow the data, so the colours stay right as values
change. OK applies the colours as one undoable step.
"""

from __future__ import annotations

from PySide6.QtGui import QColor, QDoubleValidator
from PySide6.QtWidgets import (
    QColorDialog,
    QComboBox,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QWidget,
)

from chemsuite.core.table import NUMBER, ColourRule, Table
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.table.style import CATEGORIES, GRADIENT, GRADIENT3, HIGH, LOW, MID

from ..widgets.dialogs import ToolDialog

OFF = "off"
MODES = (
    (OFF, "None"),
    (GRADIENT, "Low → high"),
    (GRADIENT3, "Low → middle → high"),
    (CATEGORIES, "One colour per value"),
)


class ColourButton(QPushButton):
    """A button showing a colour; clicking it opens the colour picker."""

    def __init__(self, title: str, rgb: int) -> None:
        super().__init__(title)
        self.title = title
        self.setAutoDefault(False)
        self.set_rgb(rgb)
        self.clicked.connect(self._pick)

    def set_rgb(self, rgb: int) -> None:
        self.rgb = rgb
        colour = QColor.fromRgb(rgb)
        text = "#000000" if colour.lightness() > 140 else "#ffffff"
        self.setStyleSheet(
            f"QPushButton {{ background: {colour.name()}; color: {text};"
            " border: 1px solid palette(mid); border-radius: 4px; padding: 3px 10px; }"
        )

    def _pick(self) -> None:
        chosen = QColorDialog.getColor(QColor.fromRgb(self.rgb), self, f"{self.title} Colour")
        if chosen.isValid():
            self.set_rgb(chosen.rgb() & 0xFFFFFF)


class ColumnColourDialog(ToolDialog):
    def __init__(self, parent: QWidget | None, table: Table, name: str) -> None:
        super().__init__(parent, f"Colour Column — {name}", guide="columns")
        self.table = table
        self.name = name
        column = table.column(name)
        self.is_number = column.kind == NUMBER
        self.data_range = column.finite_range() if isinstance(column, NumberColumn) else None
        rule = column.style.colour or ColourRule(mode=GRADIENT if self.is_number else CATEGORIES)
        current = rule.mode  # a column without colours starts on the suggested mode

        self.colours = colours = self.section("Colours")
        self.mode = QComboBox()
        for key, title in MODES:
            if key in (GRADIENT, GRADIENT3) and not self.is_number:
                continue
            self.mode.addItem(title, key)
        self.mode.setCurrentIndex(max(0, self.mode.findData(current)))
        colours.row(
            "Colour by",
            self.mode,
            "Gradients follow numbers (Number columns only); one colour per value works for any "
            "column, and the same value always gets the same colour.",
        )
        self.low = ColourButton("Low", rule.low if rule.numeric else LOW)
        self.mid = ColourButton("Middle", rule.mid if rule.numeric else MID)
        self.high = ColourButton("High", rule.high if rule.numeric else HIGH)
        swatches = QWidget()
        row = QHBoxLayout(swatches)
        row.setContentsMargins(0, 0, 0, 0)
        for button in (self.low, self.mid, self.high):
            row.addWidget(button)
        row.addStretch(1)
        self.swatches = colours.row("Colours", swatches, "Click a colour to change it.")
        self.opacity = QSpinBox()
        self.opacity.setRange(10, 100)
        self.opacity.setSuffix(" %")
        self.opacity.setValue(round(rule.opacity * 100 / 255))
        colours.row("Strength", self.opacity, "How strongly the colour covers the cell.")

        self.range = self.section("Range")
        self.minimum = self._bound(rule.minimum, 0)
        self.middle = self._bound(rule.middle, None)
        self.maximum = self._bound(rule.maximum, 1)
        tip = "Leave empty to follow the data. Values beyond the range take the end colour."
        self.range.row("Low at", self.minimum, tip)
        self.range.row("Middle at", self.middle, "Leave empty for halfway between low and high.")
        self.range.row("High at", self.maximum, tip)

        self.mode.currentIndexChanged.connect(self.refresh)
        for field in (self.minimum, self.middle, self.maximum):
            field.textChanged.connect(self.refresh)
        self.refresh()

    def _bound(self, value: float | None, end: int | None) -> QLineEdit:
        field = QLineEdit("" if value is None else f"{value:g}")
        field.setValidator(QDoubleValidator())
        if self.data_range is not None and end is not None:
            field.setPlaceholderText(f"Automatic ({self.data_range[end]:g})")
        else:
            field.setPlaceholderText("Automatic")
        return field

    # -- state ----------------------------------------------------------------------------------
    def chosen_mode(self) -> str:
        return self.mode.currentData()

    @staticmethod
    def _number(field: QLineEdit) -> float | None:
        text = field.text().strip()
        if not text:
            return None
        try:
            return float(text.replace(",", "."))
        except ValueError:
            return float("nan")

    def rule(self) -> ColourRule | None:
        mode = self.chosen_mode()
        if mode == OFF:
            return None
        opacity = round(self.opacity.value() * 255 / 100)
        if mode == CATEGORIES:
            return ColourRule(mode=CATEGORIES, opacity=opacity)
        return ColourRule(
            mode=mode,
            low=self.low.rgb,
            mid=self.mid.rgb,
            high=self.high.rgb,
            minimum=self._number(self.minimum),
            middle=self._number(self.middle) if mode == GRADIENT3 else None,
            maximum=self._number(self.maximum),
            opacity=opacity,
        )

    # -- ToolDialog -------------------------------------------------------------------------------
    def refresh(self) -> None:
        mode = self.chosen_mode()
        numeric = mode in (GRADIENT, GRADIENT3)
        self.range.setVisible(numeric)
        self.range.set_row_visible(self.middle, mode == GRADIENT3)
        self.colours.set_row_visible(self.swatches, numeric)
        self.colours.set_row_visible(self.opacity, mode != OFF)
        self.mid.setVisible(mode == GRADIENT3)
        super().refresh()

    def problem(self) -> str:
        if self.chosen_mode() not in (GRADIENT, GRADIENT3):
            return ""
        values = [self._number(f) for f in (self.minimum, self.middle, self.maximum)]
        if any(v is not None and v != v for v in values):  # NaN: not a number
            return "Range values must be numbers."
        lo, mid, hi = values
        data = self.data_range or (None, None)
        lo = lo if lo is not None else data[0]
        hi = hi if hi is not None else data[1]
        if lo is None or hi is None:
            return "The column has no numbers; give both ends of the range."
        if lo >= hi:
            return f"Low ({lo:g}) must be below high ({hi:g})."
        if self.chosen_mode() == GRADIENT3 and mid is not None and not lo <= mid <= hi:
            return "Middle must lie between low and high."
        return ""

    def verb(self) -> str:
        return "Remove Colours" if self.chosen_mode() == OFF else "Apply Colours"

    def run(self):
        column = self.table.column(self.name)
        style = column.style.with_(colour=self.rule())
        self.table.set_style(self.name, style, label=f"Colour {self.name}")
        return None
