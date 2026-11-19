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
"""A dialog that opens with the values it was last run with.

:class:`Memory` reads and writes named fields: spin boxes, combo boxes (by item data, or
text when an item has none), check boxes, segmented choices and :class:`FitChoice`. Values are
kept per user (``ToolContext.tool_options``). A value that no longer fits (a list that was
deleted, an option a newer version dropped) is skipped, so the field keeps its default.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox, QWidget

from .kit import FIELD_WIDTH, Segmented

Value = str | int | float | bool | None


class FitChoice(QComboBox):
    """Which rows a map is fitted on: All, Shown, Selected, or a named list. The map still
    covers the dialog's Rows; the rest are placed onto what these rows made."""

    def __init__(self, ctx) -> None:
        super().__init__()
        self.ctx = ctx
        self.setFixedWidth(FIELD_WIDTH)
        for key, title in (("all", "All rows"), ("shown", "Shown rows"),
                           ("selected", "Selected rows")):  # fmt: skip
            self.addItem(title, key)
        for name in ctx.lists.names:
            self.addItem(f"List: {name}", f"list:{name}")

    def oids(self) -> np.ndarray | None:
        """The rows to fit on, or None for every row the map covers."""
        from chemsuite.plotting.data import scope_oids

        key = self.currentData()
        if key == "all":
            return None
        if key.startswith("list:"):
            name = key[len("list:") :]
            return self.ctx.lists.oids(name) if name in self.ctx.lists else np.empty(0, int)
        return scope_oids(self.ctx.table, self.ctx.selection.as_array(), key)

    def rows(self, oids: np.ndarray) -> np.ndarray | None:
        """Positions in *oids* (the rows mapped) to fit on; None for all of them."""
        chosen = self.oids()
        if chosen is None:
            return None
        return np.flatnonzero(np.isin(oids, chosen))


def value_of(widget: QWidget) -> Value:
    """A field's value as plain JSON."""
    if isinstance(widget, QCheckBox):
        return widget.isChecked()
    if isinstance(widget, QSpinBox):
        return widget.value()
    if isinstance(widget, QDoubleSpinBox):
        return widget.value()
    if isinstance(widget, QComboBox):
        data = widget.currentData()
        return data if isinstance(data, (str, int, float, bool)) else widget.currentText()
    if isinstance(widget, Segmented):
        return widget.value()
    raise TypeError(f"cannot remember a {type(widget).__name__}")


def set_value(widget: QWidget, value: Value) -> bool:
    """Put *value* into a field; False (and nothing changes) when it does not fit."""
    if isinstance(widget, QCheckBox):
        if not isinstance(value, bool):
            return False
        widget.setChecked(value)
        return True
    if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
        number = isinstance(value, (int, float)) and not isinstance(value, bool)
        if not number or not widget.minimum() <= value <= widget.maximum():
            return False
        widget.setValue(int(value) if isinstance(widget, QSpinBox) else float(value))
        return True
    if isinstance(widget, QComboBox):
        index = widget.findData(value)
        if index < 0 and isinstance(value, str):
            index = widget.findText(value)
        if index < 0:
            return False
        flags = widget.model().flags(widget.model().index(index, 0))
        if not flags & Qt.ItemFlag.ItemIsEnabled:
            return False
        widget.setCurrentIndex(index)
        return True
    if isinstance(widget, Segmented):
        if value not in widget.buttons or not widget.buttons[value].isEnabled():
            return False
        widget.set_value(value)
        return True
    return False


class Memory:
    """The named *fields* of one dialog, remembered under *key*."""

    def __init__(self, ctx, key: str, fields: Mapping[str, QWidget]) -> None:
        self.ctx = ctx
        self.key = key
        self.fields = dict(fields)
        self.defaults = {name: value_of(w) for name, w in self.fields.items()}

    def restore(self) -> None:
        """Put back the values last run with (twice: a field can enable another's choices).
        A field that is off for this table (no structures for a fingerprint) keeps its value."""
        saved = self.ctx.tool_options(self.key)
        for _round in range(2):
            for name, widget in self.fields.items():
                if name in saved and widget.isEnabled():
                    set_value(widget, saved[name])

    def save(self) -> None:
        self.ctx.set_tool_options(self.key, self.values())

    def values(self) -> dict[str, Value]:
        return {name: value_of(w) for name, w in self.fields.items()}

    def reset(self) -> None:
        """Every field back to the value it had when the dialog was built."""
        for _round in range(2):
            for name, widget in self.fields.items():
                set_value(widget, self.defaults[name])
