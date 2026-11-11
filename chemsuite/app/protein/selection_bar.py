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
"""The selection field and its quick picks.

Return selects what the text says (:mod:`chemsuite.protein.select`); a mistake marks the field
and puts the cursor where reading stopped, with the reason in its tooltip. The Select menu
holds the usual picks (ligands, binding site, waters…), saved selections and Save Selection.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMenu,
    QToolButton,
    QWidget,
)

from chemsuite.protein.select import SelectionError

from .session import SET, ProteinSession

#: the Select menu's quick picks: title, expression
PICKS = (
    ("Ligands", "ligand"),
    ("Binding Site", "byres (polymer within 5 of ligand)"),
    ("Ligands and Binding Site", "ligand or byres (polymer within 5 of ligand)"),
    (None, None),
    ("Protein", "protein"),
    ("Backbone", "backbone"),
    ("Side Chains", "sidechain"),
    ("Waters", "water"),
    ("Ions", "ion"),
    ("Hetero Atoms", "hetero"),
    ("Hydrogens", "hydrogen"),
)
SITE = "byres (polymer within 5 of ligand)"
#: Select ▸ Around: distances (Å), and how the surroundings are taken
AROUND = (3, 5, 8, 10)
AROUND_KINDS = (
    ("Residues", "byres (within {} of $selection)"),
    ("Atoms", "within {} of $selection"),
)
_SYNTAX = (
    "chain A, resi 10-50, resn HIS+ASP, name CA, elem Fe; protein, ligand, water, ion, "
    "backbone, sidechain, helix, sheet; byres S, within 5 of S; and, or, not, ( ); $saved"
)


class SelectionBar(QWidget):
    selected = Signal()  # a selection was made from here

    def __init__(self, session: ProteinSession, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session = session
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        self.menu_button = QToolButton()
        self.menu_button.setText("Select")
        self.menu_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.menu = QMenu(self.menu_button)
        self.menu.aboutToShow.connect(self._fill_menu)
        self.menu_button.setMenu(self.menu)
        self.field = QLineEdit()
        self.field.setClearButtonEnabled(True)
        self.field.setPlaceholderText("Selection")
        self.field.setToolTip(_SYNTAX)
        self.field.returnPressed.connect(self.apply)
        self.field.textEdited.connect(lambda _text: self._mark(None))
        self.count = QLabel()
        self.count.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        layout.addWidget(self.menu_button)
        layout.addWidget(self.field, 1)
        layout.addWidget(self.count)
        session.selectionChanged.connect(self._update_count)
        session.proteinsChanged.connect(self._update_count)
        self._update_count()

    def apply(self) -> bool:
        text = self.field.text().strip()
        if not text:
            self.session.clear_selection()
            return True
        return self.select(text)

    def select(self, text: str, mode: str = SET) -> bool:
        try:
            self.session.select_text(text, mode)
        except SelectionError as error:
            self._mark(error)
            return False
        self._mark(None)
        self.selected.emit()
        return True

    def _mark(self, error: SelectionError | None) -> None:
        if error is None:
            self.field.setStyleSheet("")
            self.field.setToolTip(_SYNTAX)
            return
        self.field.setStyleSheet("QLineEdit { border: 1px solid palette(highlight); }")
        self.field.setToolTip(f"⚠ {error}")
        self.field.setCursorPosition(error.position)

    def _update_count(self) -> None:
        atoms = 0
        residues = 0
        for shown in self.session.displays:
            if shown.visible and shown.selected.any():
                atoms += int(shown.selected.sum())
                residues += len(set(shown.structure.residue_of[shown.selected].tolist()))
        self.count.setText(f"{atoms:,} atoms, {residues:,} residues" if atoms else "")

    def _fill_menu(self) -> None:
        self.menu.clear()
        for title, expression in PICKS:
            if title is None:
                self.menu.addSeparator()
                continue
            action = QAction(title, self.menu)
            action.triggered.connect(lambda _=False, e=expression: self._pick(e))
            self.menu.addAction(action)
        chains = self.menu.addMenu("Chain")
        names = sorted({str(n) for d in self.session.displays if d.visible
                        for n in d.structure.chain_name})  # fmt: skip
        for name in names:
            chains.addAction(name, lambda n=name: self._pick(f"chain {n}"))
        chains.setEnabled(bool(names))
        around = self.menu.addMenu("Around")
        for kind, template in AROUND_KINDS:
            if around.actions():
                around.addSeparator()
            for distance in AROUND:
                expression = template.format(distance)
                around.addAction(f"{kind} Within {distance} Å",
                                 lambda e=expression: self.grow(e))  # fmt: skip
        around.setToolTip("Add what lies near the selection to it")
        around.setEnabled(self.session.selected_count() > 0)
        self.menu.addSeparator()
        saved = self.menu.addMenu("Saved Selections")
        for name in sorted(self.session.named):
            saved.addAction(name, lambda n=name: self._pick(f"${n}"))
        saved.setEnabled(bool(self.session.named))
        save = self.menu.addAction("Save Selection…", self.save)
        save.setEnabled(self.session.selected_count() > 0)
        self.menu.addSeparator()
        clear = self.menu.addAction("Clear Selection", self.session.clear_selection)
        clear.setEnabled(self.session.selected_count() > 0)

    def grow(self, expression: str) -> None:
        """Select by an expression of the selection itself (``$selection``): the field is
        cleared, as its text no longer says what is selected."""
        if self.select(expression):
            self.field.clear()

    def _pick(self, expression: str) -> None:
        self.field.setText(expression)
        self.select(expression)

    def save(self) -> None:
        name, ok = QInputDialog.getText(self, "Save Selection", "Name:")
        if not ok or not name.strip():
            return
        try:
            self.session.save_selection(name.strip())
        except SelectionError as error:
            self._mark(error)
