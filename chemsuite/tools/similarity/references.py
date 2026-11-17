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
"""The references a similarity search compares with: a short list of structures (typed or
drawn, the current row, a protein's ligand), or the selected rows, followed live."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from PySide6.QtCore import QObject, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QLineEdit,
    QListWidget,
    QMenu,
    QPushButton,
    QToolButton,
)
from rdkit import Chem

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.shell.glyphs import glyph
from chemsuite.app.widgets.dialogs import Segmented, follow_rows, inline
from chemsuite.core.mols import smiles_of
from chemsuite.core.structure_input import read_structure
from chemsuite.similarity.model import FUSIONS, Reference
from chemsuite.superpose.run import mol_of_cell

from .common import labels

LISTED, SELECTION = "listed", "selection"


@dataclass(frozen=True)
class Entry:
    """A reference as the dialog holds it: the structure, its name, the row it is (None: not a
    row of the table) and what the reference asset of a superposed pose records about it."""

    mol: Chem.Mol
    name: str
    oid: int | None = None
    meta: dict[str, str] = field(default_factory=dict)


class References(QObject):
    """The References section's fields. *structure* gives the structure column rows are read
    from (the search's own)."""

    changed = Signal()

    def __init__(self, ctx: ToolContext, section, structure: Callable[[], str]) -> None:
        super().__init__(section)
        self.ctx = ctx
        self.section = section
        self.structure = structure
        self._entries: list[Entry] = []
        self._selection_cache: tuple[tuple, list[Entry]] | None = None

        self.source = Segmented(((LISTED, "Structures"), (SELECTION, "Selected Rows")))
        section.row(
            "Compare with",
            self.source,
            help="Structures: the list below. Selected Rows: every selected row, followed as "
            "the selection changes; a selected row is not compared with itself.",
        )
        follow_rows(
            self.source,
            ctx,
            lambda: self.source.set_text(SELECTION, f"Selected Rows ({len(ctx.selection):,})"),
        )
        self.smiles = QLineEdit()
        self.smiles.setPlaceholderText("SMILES")
        self.smiles.setClearButtonEnabled(True)
        draw = QToolButton()
        draw.setAutoRaise(True)
        draw.setFixedSize(22, 22)
        draw.setIcon(glyph("pencil", draw.palette().text().color()))
        draw.setToolTip("Draw a structure in the sketcher")
        self.add_button = QPushButton("Add")
        self.add_button.setAutoDefault(False)
        self.add_field = inline(self.smiles, draw, self.add_button, stretch=0)
        section.row("Structure", self.add_field, help="Type or paste SMILES, or draw it, then Add.")

        self.list = QListWidget()
        self.list.setFixedHeight(76)
        section.row("References", self.list)
        self.current = QPushButton("Current Row")
        self.current.setAutoDefault(False)
        self.current.setToolTip("Add the table's current row")
        self.ligand = QPushButton("Ligand")
        self.ligand.setAutoDefault(False)
        self.ligand.setToolTip("Add a ligand of a protein in the document, in its bound pose")
        self._ligands = self._ligand_refs()
        menu = QMenu(self.ligand)
        for ref in self._ligands:
            menu.addAction(ref.title, lambda r=ref: self.add_ligand(r))
        self.ligand.setMenu(menu)
        self.ligand.setVisible(bool(self._ligands))
        self.remove = QPushButton("Remove")
        self.remove.setAutoDefault(False)
        self.buttons = inline(self.current, self.ligand, self.remove)
        section.row("", self.buttons)

        self.fuse = QComboBox()
        for key, title in FUSIONS.items():
            self.fuse.addItem(title, key)
        section.row(
            "Fuse references",
            self.fuse,
            help="How a row's scores against several references become one: the best of them, "
            "or their mean. Nearest Reference names the one it is most like.",
        )

        draw.clicked.connect(self.draw)
        self.add_button.clicked.connect(self.add_typed)
        self.smiles.returnPressed.connect(self.add_typed)
        self.smiles.textChanged.connect(lambda _t: self._update())
        self.current.clicked.connect(self.add_current_row)
        self.remove.clicked.connect(self.remove_chosen)
        self.list.itemSelectionChanged.connect(self._update)
        self.source.changed.connect(lambda _k: self._update())
        self.fuse.currentIndexChanged.connect(lambda _i: self.changed.emit())
        follow_rows(self.list, ctx, self._selection_changed)
        if ctx.current_oid() is not None:
            self.add_current_row()
        self._update()

    def _ligand_refs(self) -> list:
        from chemsuite.tools.macromolecule.ligand_refs import ligand_refs

        return ligand_refs(self.ctx)

    # -- adding and removing ------------------------------------------------------------------
    def add(self, entry: Entry) -> None:
        self._entries.append(entry)
        self.list.addItem(entry.name)
        self._update()

    def typed_mol(self) -> tuple[Chem.Mol | None, str]:
        """The structure in the SMILES field and, if it does not read, why."""
        found = read_structure(self.smiles.text())
        return found.mol, ("" if found.mol is not None else found.error or "Enter a structure.")

    def add_typed(self) -> None:
        mol, _why = self.typed_mol()
        if mol is None:
            return
        smiles = smiles_of(mol)
        self.add(Entry(mol, smiles, None, {"source": "structure", "smiles": smiles}))
        self.smiles.clear()

    def add_current_row(self) -> None:
        oid, column = self.ctx.current_oid(), self.structure()
        if oid is None or not column:
            return
        mol = mol_of_cell(self.ctx.structure_inputs(column, [oid])[0])
        if mol is None:
            return
        name = labels(self.ctx, np.array([oid]))[0]
        self.add(Entry(mol, name, oid, {"source": "row", "oid": str(oid), "column": column}))

    def add_ligand(self, ref) -> None:
        from chemsuite.tools.macromolecule.ligand_refs import ligand_mol

        mol = ligand_mol(self.ctx, ref)
        if mol is None:
            return
        meta = {"source": "protein", "protein": ref.protein, "ligand": ref.label}
        self.add(Entry(mol, ref.title, None, meta))

    def remove_chosen(self) -> None:
        for row in sorted((self.list.row(i) for i in self.list.selectedItems()), reverse=True):
            self.list.takeItem(row)
            del self._entries[row]
        self._update()

    def draw(self) -> None:
        from chemsuite.app.sketcher import SketcherDialog

        dialog = SketcherDialog(self.smiles.text(), "Reference", self.list.window())
        if not dialog.exec():
            return
        found = read_structure(dialog.molblock) if dialog.molblock else None
        if found is not None and found.mol is not None:
            self.smiles.setText(smiles_of(found.mol))

    # -- state --------------------------------------------------------------------------------
    def from_selection(self) -> bool:
        return self.source.value() == SELECTION

    def _selection_changed(self) -> None:
        if self.from_selection():
            self.changed.emit()

    def _update(self) -> None:
        listed = not self.from_selection()
        for widget in (self.add_field, self.list, self.buttons):
            self.section.set_row_visible(widget, listed)
        self.add_button.setEnabled(self.typed_mol()[0] is not None)
        self.current.setEnabled(self.ctx.current_oid() is not None)
        self.remove.setEnabled(bool(self.list.selectedItems()))
        self.changed.emit()

    def set_fuse_visible(self, visible: bool) -> None:
        self.section.set_row_visible(self.fuse, visible)

    def fusion(self) -> str:
        return self.fuse.currentData()

    def entries(self) -> list[Entry]:
        """The references as they are now: the list, or every selected row that reads."""
        if not self.from_selection():
            return list(self._entries)
        column = self.structure()
        oids = self.ctx.selection.as_array()
        key = (column, oids.tobytes())
        if self._selection_cache is not None and self._selection_cache[0] == key:
            return self._selection_cache[1]
        found: list[Entry] = []
        if column and len(oids):
            names = labels(self.ctx, oids)
            cells = self.ctx.structure_inputs(column, oids)
            for oid, name, cell in zip(oids.tolist(), names, cells, strict=True):
                mol = mol_of_cell(cell)
                if mol is not None:
                    found.append(Entry(mol, name, oid, {"source": "row", "oid": str(oid),
                                                        "column": column}))  # fmt: skip
        self._selection_cache = (key, found)
        return found

    def references(
        self, targets: np.ndarray, entries: list[Entry] | None = None
    ) -> list[Reference]:
        """The references for a search of the rows *targets* (in that order). A selected row
        that is one of them knows its position there, so it is not compared with itself; a
        listed one (the current row, say) is, and scores as itself, as a typed structure does."""
        entries = self.entries() if entries is None else entries
        if not self.from_selection():
            return [Reference(e.mol, e.name) for e in entries]
        position = {oid: i for i, oid in enumerate(np.asarray(targets).tolist())}
        return [Reference(e.mol, e.name, position.get(e.oid, -1)) for e in entries]

    def problem(self) -> str:
        if self.from_selection():
            if not len(self.ctx.selection):
                return "No rows are selected. Select the compounds to compare with."
            if not self.entries():
                return "None of the selected rows has a structure."
            return ""
        if not self._entries:
            return "Add a reference: a structure (type or draw it, then Add) or the current row."
        return ""
