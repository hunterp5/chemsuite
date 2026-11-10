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

"""Number Format…: how a Number column shows its values (display only, undoable)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QComboBox, QSpinBox, QWidget

from chemsuite.core.table import Table
from chemsuite.core.table.style import (
    FIXED,
    MAX_PRECISION,
    PERCENT,
    SCIENTIFIC,
    SIGNIFICANT,
    THOUSANDS,
    format_with,
)
from chemsuite.core.table.summary import FORMAT_TITLES

from .dialogs import ToolDialog

AS_TYPED = "as_typed"
#: (key, title, example value for the tooltip)
CHOICES = (
    (AS_TYPED, "As typed", None),
    (FIXED, FORMAT_TITLES[FIXED], 1234.5678),
    (THOUSANDS, FORMAT_TITLES[THOUSANDS], 1234567.891),
    (SCIENTIFIC, FORMAT_TITLES[SCIENTIFIC], 0.000123456),
    (SIGNIFICANT, FORMAT_TITLES[SIGNIFICANT], 1234.5678),
    (PERCENT, FORMAT_TITLES[PERCENT], 0.4567),
)
DEFAULT_DIGITS = {FIXED: 2, THOUSANDS: 0, SCIENTIFIC: 2, SIGNIFICANT: 3, PERCENT: 1}


class NumberFormatDialog(ToolDialog):
    """Choose a format and its digits for one Number column."""

    def __init__(self, parent: QWidget | None, table: Table, name: str) -> None:
        super().__init__(parent, f"Number Format — {name}", guide="columns")
        self.table = table
        self.name = name
        style = table.column(name).style

        section = self.section("Display")
        self.format = QComboBox()
        for key, title, example in CHOICES:
            self.format.addItem(title, key)
            if example is not None:
                sample = format_with(example, key, DEFAULT_DIGITS[key])
                self.format.setItemData(
                    self.format.count() - 1,
                    f"{example:g} is shown as {sample}",
                    Qt.ItemDataRole.ToolTipRole,
                )
        section.row(
            "Format",
            self.format,
            help="How the column shows its numbers; the values themselves are unchanged.",
        )
        self.digits = QSpinBox()
        self.digits.setRange(0, MAX_PRECISION)
        section.row("Decimals", self.digits)
        self.digits_label = section.caption(self.digits)
        current = style.number_format if style.formats_numbers else AS_TYPED
        self.format.setCurrentIndex(self.format.findData(current))
        self.digits.setValue(style.digits if style.formats_numbers else 2)
        self.format.currentIndexChanged.connect(self._format_changed)
        self._show_digits()
        self.refresh()

    def verb(self) -> str:
        return "Apply"

    def run(self):
        self.apply()

    def chosen(self) -> str:
        return self.format.currentData()

    def _format_changed(self, _index: int) -> None:
        key = self.chosen()
        if key in DEFAULT_DIGITS:
            self.digits.setValue(DEFAULT_DIGITS[key])
        self._show_digits()

    def _show_digits(self) -> None:
        key = self.chosen()
        significant = key == SIGNIFICANT
        self.digits_label.setText("Significant figures:" if significant else "Decimals:")
        self.digits.setMinimum(1 if significant else 0)
        self.digits.setEnabled(key != AS_TYPED)

    def apply(self) -> bool:
        """Set the chosen format on the column as one undoable step; False if unchanged."""
        column = self.table.column(self.name)
        key = self.chosen()
        if key == AS_TYPED:
            style = column.style.with_(precision=None, number_format=FIXED)
        else:
            style = column.style.with_(precision=self.digits.value(), number_format=key)
        if style == column.style:
            return False
        self.table.set_style(self.name, style, label=f"Number format of {self.name}")
        return True
