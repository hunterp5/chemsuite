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
"""The Browser's conformer list: the shown row's conformers with their stored energies.

One line per conformer: its number, ΔE above the lowest, Boltzmann population and energy
(:mod:`chemsuite.conformers.energies`). It is linked to the 3D view both ways: picking a line
shows that conformer, picking several overlays just those, and stepping through conformers
moves the highlight. Sorting a column changes only the order of the lines.
"""

from __future__ import annotations

import math
from collections.abc import Callable

from PySide6.QtCore import QItemSelectionModel, Qt, Signal
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QLabel,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QToolBar,
    QWidget,
)
from rdkit import Chem

from chemsuite.conformers import energies

HEADERS = ("#", "ΔE (kcal/mol)", "Population (%)", "Energy (kcal/mol)")
_POSITION = Qt.ItemDataRole.UserRole  # the conformer's position in the molecule (0 first)
_VALUE = Qt.ItemDataRole.UserRole + 1  # the number a column sorts by


class _Number(QTableWidgetItem):
    """A cell that sorts by its number, not its text."""

    def __lt__(self, other: QTableWidgetItem) -> bool:
        return float(self.data(_VALUE)) < float(other.data(_VALUE))


def counter_text(total: int, conformer: int, overlay: bool, chosen: tuple[int, ...],
                 delta: float) -> str:  # fmt: skip
    """The toolbar's counter: "3 / 50", with the shown conformer's ΔE when it has one."""
    if overlay:
        return f"{total} Conformers"
    if len(chosen) > 1:
        return f"{len(chosen)} of {total} Conformers"
    text = f"{conformer + 1} / {total}"
    return text if math.isnan(delta) else f"{text}   ΔE {delta:+.2f}"


class ConformerList(QTableWidget):
    """The conformers of the row shown, with their energies."""

    #: positions picked (0 first): one shows that conformer, several are overlaid
    picked = Signal(tuple)

    def __init__(self) -> None:
        super().__init__(0, len(HEADERS))
        self.setHorizontalHeaderLabels(list(HEADERS))
        self.verticalHeader().setVisible(False)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.setAlternatingRowColors(True)
        self.setWordWrap(False)
        header = self.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        header.setHighlightSections(False)
        header.setSortIndicator(0, Qt.SortOrder.AscendingOrder)  # as stored, lowest first
        self.verticalHeader().setDefaultSectionSize(self.fontMetrics().height() + 6)
        self.setMinimumHeight(4 * (self.fontMetrics().height() + 6))
        self._key: object = None
        self._filling = False
        self.itemSelectionChanged.connect(self._selection_changed)

    def has_energies(self) -> bool:
        return self.rowCount() > 0

    def show_conformers(self, key: object, mol: Chem.Mol | None) -> bool:
        """List *mol*'s conformers, unless *key* (which row and column) is listed already.
        Whether a new molecule replaced the one listed before (not the first one shown)."""
        if key == self._key:
            return False
        replaced = self._key is not None
        self._key = key
        rows = energies.rows(mol)
        self._filling = True
        self.setSortingEnabled(False)
        self.clearContents()
        self.setRowCount(len(rows))
        for line, row in enumerate(rows):
            values = (
                (str(row.number), row.number),
                (f"{row.delta:.2f}", row.delta),
                (f"{100 * row.population:.1f}", row.population),
                (f"{row.energy:.2f}", row.energy),
            )
            for column, (text, value) in enumerate(values):
                item = _Number(text)
                item.setData(_POSITION, row.number - 1)
                item.setData(_VALUE, value)
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.setItem(line, column, item)
        self.setSortingEnabled(True)
        self._filling = False
        self._tooltips(energies.force_field(mol))
        return replaced

    def _tooltips(self, force_field: str) -> None:
        source = f"{force_field}, " if force_field else ""
        tips = (
            "The conformer's number, lowest energy first as stored",
            f"Energy above the lowest conformer ({source}kcal/mol)",
            f"Boltzmann population at {energies.TEMPERATURE} K over the kept conformers "
            f"({source}gas phase)",
            f"The conformer's energy ({source}kcal/mol); compare within one molecule only",
        )
        for column, tip in enumerate(tips):
            item = self.horizontalHeaderItem(column)
            if item is not None:
                item.setToolTip(tip)

    def set_current(self, position: int, chosen: tuple[int, ...] = ()) -> None:
        """Highlight the conformer at *position* (or the *chosen* ones) without signalling."""
        self._filling = True
        model = self.selectionModel()
        model.clearSelection()
        wanted = set(chosen) if len(chosen) > 1 else {position}
        flags = QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows
        for line in range(self.rowCount()):
            item = self.item(line, 0)
            if item is not None and item.data(_POSITION) in wanted:
                model.select(self.model().index(line, 0), flags)
                if item.data(_POSITION) == position:
                    self.scrollToItem(item)
        self._filling = False

    def _selection_changed(self) -> None:
        if self._filling:
            return
        lines = sorted({index.row() for index in self.selectionModel().selectedRows()})
        positions = tuple(sorted(int(self.item(line, 0).data(_POSITION)) for line in lines))
        if positions:
            self.picked.emit(positions)


class ConformerBar:
    """The Browser's conformer controls: Previous Conformer, the counter, Next Conformer and
    the Conformer Energies toggle on its toolbar (shown for a row with several conformers, in
    3D), and the list under the view. *chosen* holds the conformers picked to overlay."""

    def __init__(
        self,
        bar: QToolBar,
        add: Callable[..., QAction],
        step: Callable[[int], None],
        picked: Callable[[tuple[int, ...]], None],
        changed: Callable[..., None],
    ) -> None:
        self.previous = add("previous", "Previous Conformer", "one back (↑)", lambda: step(-1))
        self.counter = QLabel()
        self.counter.setToolTip("Conformer shown, of how many, and its energy above the lowest")
        self.counter.setContentsMargins(4, 0, 4, 0)
        self.actions = [self.previous, bar.addWidget(self.counter)]
        self.next = add("next", "Next Conformer", "one on (↓)", lambda: step(1))
        self.energies = add("energies", "Conformer Energies", "list the row's conformers with "
                            "their energies and populations; pick one to show it, several to "
                            "overlay them", changed, checkable=True)  # fmt: skip
        blocked = self.energies.blockSignals(True)  # on from the start, before the panel is
        self.energies.setChecked(True)
        self.energies.blockSignals(blocked)
        self.actions += [self.next, self.energies, bar.addSeparator()]
        self.list = ConformerList()
        self.list.picked.connect(picked)
        self.chosen: tuple[int, ...] = ()
        self._mol: Chem.Mol | None = None

    def beside(self, view: QWidget) -> QSplitter:
        """*view* above the list, the boundary between them dragged to resize."""
        split = QSplitter(Qt.Orientation.Vertical)
        split.setChildrenCollapsible(False)
        split.addWidget(view)
        split.addWidget(self.list)
        split.setStretchFactor(0, 1)
        return split

    def show(self, key: object, mol: Chem.Mol | None) -> None:
        """*mol* is the row shown (*key*: which row and column); a new one drops the picks."""
        if self.list.show_conformers(key, mol):
            self.chosen = ()
        self._mol = mol

    def update(self, total: int, conformer: int, overlay: bool, rows: bool) -> None:
        """Show *total* conformers with the one at *conformer* current (*total* is 0 in 2D);
        *rows* is Overlay Rows, where one row's energies mean nothing."""
        delta = math.nan if rows else energies.delta_of(self._mol, conformer)
        self.counter.setText(counter_text(total, conformer, overlay, self.chosen, delta))
        for action in self.actions:
            action.setVisible(total > 1)
        self.list.setVisible(total > 1 and self.energies.isChecked() and not rows
                             and self.list.has_energies())  # fmt: skip
        self.list.set_current(conformer, () if overlay else self.chosen)
        stepping = total > 1 and not overlay
        self.previous.setEnabled(stepping and conformer > 0)
        self.next.setEnabled(stepping and conformer < total - 1)

    def state(self) -> dict:
        return {"chosen": list(self.chosen), "energies": self.energies.isChecked()}

    def set_state(self, state: dict) -> None:
        chosen = state.get("chosen")
        self.chosen = tuple(int(i) for i in chosen) if isinstance(chosen, list) else ()
        self.energies.setChecked(bool(state.get("energies", True)))
