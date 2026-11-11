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

"""A reaction (or a structure) typed, drawn in the sketcher or taken from the current row: the
field the reaction tools and the reaction filter card share."""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QLineEdit, QPushButton, QVBoxLayout, QWidget

from chemsuite.app.dialogs.kit import FIELD_WIDTH
from chemsuite.core.reactions import reaction_from_text, reaction_smiles


class ReactionText(QWidget):
    """A reaction typed as reaction SMILES or SMARTS, or drawn in the sketcher (which comes
    back as reaction SMILES with its atom maps). With ``reactions=False`` it takes a structure
    instead, for a structure query. *current* gives the current row's reaction, for a
    "Current Row" button."""

    changed = Signal()

    def __init__(
        self,
        placeholder: str = "Reaction SMILES or SMARTS",
        current: Callable[[], str] | None = None,
        parent: QWidget | None = None,
        *,
        dialog_width: bool = True,
    ) -> None:
        super().__init__(parent)
        self.reactions = True
        self._current = current
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.setClearButtonEnabled(True)
        self.draw_button = QPushButton("Draw…")
        self.draw_button.setAutoDefault(False)
        # The text takes the dialog's field width; the buttons sit under it, so the field's
        # right edge lines up with every other field.
        column = QVBoxLayout(self)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(4)
        column.addWidget(self.edit)
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(self.draw_button)
        self.current_button: QPushButton | None = None
        if current is not None:
            self.current_button = QPushButton("Current Row")
            self.current_button.setAutoDefault(False)
            self.current_button.setToolTip("Use the reaction in the table's current row")
            self.current_button.clicked.connect(self.use_current_row)
            row.addWidget(self.current_button)
        row.addStretch(1)
        column.addLayout(row)
        if dialog_width:  # in a tool dialog; a filter card sets its own width
            self.setFixedWidth(FIELD_WIDTH)
        self.edit.textChanged.connect(lambda _t: self.changed.emit())
        self.draw_button.clicked.connect(self.draw)

    def text(self) -> str:
        return self.edit.text().strip()

    def set_text(self, text: str) -> None:
        self.edit.setText(text)

    def set_reactions(self, on: bool) -> None:
        """Take a reaction (True) or a structure (False)."""
        self.reactions = on
        self.edit.setPlaceholderText("Reaction SMILES or SMARTS" if on else "SMILES or SMARTS")
        if self.current_button is not None:
            self.current_button.setVisible(on)

    def use_current_row(self) -> None:
        text = self._current() if self._current else ""
        if text:
            self.edit.setText(text)

    def draw(self) -> None:
        from chemsuite.app.sketcher import SketcherDialog

        if not self.reactions:
            return self._draw_structure()
        dialog = SketcherDialog(self.text(), "Draw Reaction", self, reaction_only=True)
        if not dialog.exec() or not dialog.rxn:
            return None
        rxn = reaction_from_text(dialog.rxn)
        if rxn is not None:
            self.edit.setText(reaction_smiles(rxn))
        return None

    def _draw_structure(self) -> None:
        from chemsuite.analysis.structure_search import drawn_query_text
        from chemsuite.app.sketcher import SketcherDialog
        from chemsuite.core.structure_input import read_structure

        dialog = SketcherDialog(self.text(), "Query Structure", self)
        if not dialog.exec() or not dialog.molblock:
            return
        parsed = read_structure(dialog.molblock)
        if parsed.mol is not None:
            self.edit.setText(drawn_query_text(parsed.mol))
