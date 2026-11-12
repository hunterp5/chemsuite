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
"""The proteins list: each protein (its box shows or hides it), its chains, its ligands, and
its ions, waters and additives.

Clicking an entry selects what it stands for (Shift adds, Ctrl or Cmd toggles); a double-click
also frames it in the view (``focusRequested``).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QKeySequence, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QMenu,
    QTreeWidget,
    QTreeWidgetItem,
    QWidget,
)

from chemsuite.protein.chem.ligands import GUESSED
from chemsuite.protein.structure.model import (
    ADDITIVE,
    AMINO,
    GLYCAN,
    ION,
    LIGAND,
    NUCLEIC,
    WATER,
    Structure,
)
from chemsuite.protein.structure.select import select

from .session import ADD, SET, TOGGLE, ProteinSession

_PROTEIN = Qt.ItemDataRole.UserRole
_SELECTS = Qt.ItemDataRole.UserRole + 1
_GROUPS = ((ION, "Ions", "ion"), (WATER, "Waters", "water"), (ADDITIVE, "Additives", "additive"),
           (GLYCAN, "Glycans", "glycan"))  # fmt: skip


def residue_expression(structure: Structure, residue: int) -> str:
    chain = structure.chain_name[structure.chain_of[residue]]
    number = f"{structure.res_seq[residue]}{structure.icode[residue]}"
    return f"chain {chain} and resi {number} and resn {structure.res_name[residue]}"


class ProteinTree(QTreeWidget):
    focusRequested = Signal()
    removeRequested = Signal(str)  # a protein's id: Remove Protein from its menu or Delete

    def __init__(self, session: ProteinSession, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session = session
        self.setHeaderHidden(True)
        self.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.setUniformRowHeights(True)
        self.setExpandsOnDoubleClick(False)
        self._filling = False
        self._ticking = False  # applying a protein's check box: the tree already shows it
        self.itemChanged.connect(self._item_changed)
        self.itemClicked.connect(self._item_clicked)
        self.itemDoubleClicked.connect(self._item_double_clicked)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._context_menu)
        self._masks: dict[int, np.ndarray] = {}  # id(item) -> the atoms it selects, when asked
        session.proteinsChanged.connect(self.reload)
        session.selectionChanged.connect(self.show_selection)
        self.reload()

    def reload(self) -> None:
        if self._ticking:
            # Qt is still inside the ticked item's setData: clearing would delete it under Qt
            return
        expanded = {
            self.topLevelItem(i).data(0, _PROTEIN)
            for i in range(self.topLevelItemCount())
            if self.topLevelItem(i).isExpanded()
        }
        self._filling = True
        self.clear()
        for shown in self.session.displays:
            top = QTreeWidgetItem([shown.name])
            top.setData(0, _PROTEIN, shown.id)
            top.setData(0, _SELECTS, "all")
            top.setFlags(top.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            top.setCheckState(
                0, Qt.CheckState.Checked if shown.visible else Qt.CheckState.Unchecked
            )
            meta = shown.structure.meta
            facts = [meta.get("title", ""), meta.get("method", "").title(),
                     f"{meta['resolution']} Å" if "resolution" in meta else ""]  # fmt: skip
            top.setToolTip(0, "\n".join(f for f in facts if f))
            self._fill(top, shown.id, shown.structure)
            self.addTopLevelItem(top)
            top.setExpanded(shown.id in expanded or len(self.session.displays) == 1)
        self._filling = False
        self._masks = {}
        self.show_selection()

    def show_selection(self) -> None:
        """Mark the items whose atoms are all selected (highlighted) or some (bold)."""
        palette = self.palette()
        fill = palette.color(QPalette.ColorRole.Highlight)
        ink = palette.color(QPalette.ColorRole.HighlightedText)
        self._filling = True
        blocked = self.blockSignals(True)
        try:
            for i in range(self.topLevelItemCount()):
                top = self.topLevelItem(i)
                for item in [top] + [top.child(k) for k in range(top.childCount())]:
                    self._mark(item, fill, ink)
        finally:
            self.blockSignals(blocked)
            self._filling = False

    def _mark(self, item: QTreeWidgetItem, fill: QColor, ink: QColor) -> None:
        shown = self.session.display(str(item.data(0, _PROTEIN)))
        some = every = False
        if shown is not None and shown.selected.any():
            mask = self._masks.get(id(item))
            if mask is None:
                mask = select(shown.structure, str(item.data(0, _SELECTS))) & shown.active
                self._masks[id(item)] = mask
            chosen = shown.selected[mask]
            some = bool(chosen.any())
            every = some and bool(chosen.all())
        font = item.font(0)
        font.setBold(some and not every)
        item.setFont(0, font)
        item.setBackground(0, QBrush(fill) if every else QBrush())
        item.setForeground(0, QBrush(ink) if every else QBrush())

    def _fill(self, top: QTreeWidgetItem, protein: str, structure: Structure) -> None:
        kinds = structure.res_kind
        for c in range(structure.n_chains):
            span = structure.residues_of(c)
            polymer = int(np.isin(kinds[span], (AMINO, NUCLEIC)).sum())
            if polymer:
                name = structure.chain_name[c]
                self._child(top, protein, f"Chain {name} ({polymer})",
                            f"chain {name} and polymer")  # fmt: skip
        asset = self.session.table.assets.get(protein)
        guessed = set(asset.meta.get(GUESSED, "").splitlines()) if asset is not None else set()
        for residue in np.flatnonzero(kinds == LIGAND):
            chain = structure.chain_name[structure.chain_of[residue]]
            label = f"{structure.res_name[residue]} {chain}{structure.res_seq[residue]}"
            ligand_id = structure.residue_label(int(residue))
            if ligand_id in guessed:
                label += " ⚠"
                ligand_id += "\nBond orders guessed from the geometry"
            item = self._child(top, protein, label, residue_expression(structure, int(residue)))
            item.setToolTip(0, ligand_id)
        for kind, title, expression in _GROUPS:
            count = int((kinds == kind).sum())
            if count:
                self._child(top, protein, f"{title} ({count})", expression)

    def _child(self, top: QTreeWidgetItem, protein: str, text: str, selects: str):
        item = QTreeWidgetItem([text])
        item.setData(0, _PROTEIN, protein)
        item.setData(0, _SELECTS, selects)
        top.addChild(item)
        return item

    def _item_changed(self, item: QTreeWidgetItem, _column: int) -> None:
        if self._filling or item.parent() is not None:
            return
        visible = item.checkState(0) == Qt.CheckState.Checked
        self._ticking = True
        try:
            self.session.set_visible(item.data(0, _PROTEIN), visible)
        finally:
            self._ticking = False

    def _item_clicked(self, item: QTreeWidgetItem, _column: int) -> None:
        modifiers = QApplication.keyboardModifiers()
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            mode = ADD
        elif modifiers & Qt.KeyboardModifier.ControlModifier:
            mode = TOGGLE
        else:
            mode = SET
        self.select_item(item, mode)

    def select_item(self, item: QTreeWidgetItem, mode: str = SET) -> None:
        protein = item.data(0, _PROTEIN)
        shown = self.session.display(protein)
        if shown is None or not shown.visible:
            return
        self.session.select_text(item.data(0, _SELECTS), mode, protein=protein)

    def _item_double_clicked(self, item: QTreeWidgetItem, _column: int) -> None:
        self.select_item(item, SET)
        self.focusRequested.emit()

    def _context_menu(self, where) -> None:
        item = self.itemAt(where)
        if item is None:
            return
        protein = str(item.data(0, _PROTEIN))
        menu = QMenu(self)
        menu.addAction("Remove Protein…", lambda: self.removeRequested.emit(protein))
        menu.exec(self.viewport().mapToGlobal(where))

    def keyPressEvent(self, event) -> None:
        if event.matches(QKeySequence.StandardKey.SelectAll) or event.matches(
            QKeySequence.StandardKey.Copy
        ):
            event.ignore()  # the viewer's: select every atom, copy the selection
            return
        item = self.currentItem()
        if item is not None and event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.removeRequested.emit(str(item.data(0, _PROTEIN)))
            return
        super().keyPressEvent(event)

    def changeEvent(self, event) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange:
            self.show_selection()  # the highlight's colours, light or dark
