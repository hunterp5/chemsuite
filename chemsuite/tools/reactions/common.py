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

"""What the reaction tools share: the current row's reaction, reaction cells, the named
template chooser and the one-off building block."""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.reaction_templates import GROUPS, Template, by_id, in_group
from chemsuite.app.widgets.reaction_text import ReactionText as ReactionText

CUSTOM = ""  # the chooser's data for "Custom", a template typed or drawn by hand


class TemplateChoice(QComboBox):
    """Pick a named reaction template, or Custom to type or draw one.

    The named templates are grouped (Coupling, Ring Formation, …) with their group as a
    disabled heading, the way a menu separates its sections.
    """

    picked = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.addItem("Custom", CUSTOM)
        for group in GROUPS:
            templates = in_group(group)
            if not templates:
                continue
            self.addItem(group, None)
            index = self.count() - 1
            self.model().item(index).setEnabled(False)
            for template in templates:
                self.addItem(f"    {template.name}", template.id)
                self.setItemData(self.count() - 1, template.describe(), 3)  # Qt::ToolTipRole
        self.activated.connect(lambda _index: self.picked.emit())

    def template(self) -> Template | None:
        """The chosen template, or None for Custom."""
        chosen = self.currentData()
        return by_id(chosen) if chosen else None

    def set_template(self, template_id: str) -> None:
        index = self.findData(template_id or CUSTOM)
        if index >= 0:
            self.setCurrentIndex(index)

    def show_custom(self) -> None:
        """Fall back to Custom, without emitting :attr:`picked` (the text already stands)."""
        self.setCurrentIndex(0)


class StructureDialog(QDialog):
    """One building block, typed as SMILES or drawn, for a component that is a single
    structure rather than a column."""

    def __init__(self, parent: QWidget | None, role: str, text: str = "") -> None:
        super().__init__(parent)
        self.setWindowTitle(role)
        self.field = ReactionText("SMILES", dialog_width=False)
        self.field.set_reactions(False)
        self.field.set_text(text)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(self.field)
        layout.addWidget(buttons)

    def text(self) -> str:
        return self.field.text()


def ask_structure(parent: QWidget | None, role: str, text: str = "") -> str:
    """A single building block from the user, or "" when cancelled."""
    dialog = StructureDialog(parent, role, text)
    return dialog.text() if dialog.exec() else ""


def current_reaction(ctx, column: Callable[[], str]) -> Callable[[], str]:
    """A function giving the reaction text of *column()* in the table's current row."""

    def text() -> str:
        name, oid = column(), ctx.current_oid()
        table = ctx.table
        if not name or name not in table or oid is None or oid not in table.oids:
            return ""
        return table.cell(name, oid)

    return text


def cells(table, column: str, oids):
    """Cells to send to workers: the parsed blob, else the text (workers parse it)."""
    data = table.column(column)
    slots = table.slots_of(oids)
    return [
        blob if blob is not None else (text or None)
        for blob, text in zip(data.blobs(slots), data.displays(slots), strict=True)
    ]
