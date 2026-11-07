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
"""The template library window (Shift+T, or the last button of the ring bar): pick a ring
system to draw with, by group or by name."""

from __future__ import annotations

from functools import lru_cache

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLineEdit,
    QListView,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.sketch import library

ALL = "All"
TILE = QSize(104, 84)
ZOOM = 1.35  # drawings a little larger than in table cells


@lru_cache(maxsize=256)
def _drawing(smiles: str, dark: bool, ratio: float) -> QPixmap:
    from rdkit import Chem

    from chemsuite.core.depict import fitted_size, render_png

    # a little larger than the table's bond length, so small rings stay small but readable;
    # the list centres each drawing in its tile
    mol = Chem.MolFromSmiles(smiles)
    width, height = fitted_size(mol, round(TILE.width() / ZOOM), round(TILE.height() / ZOOM))
    width, height = round(width * ZOOM), round(height * ZOOM)
    drawing = QPixmap()
    drawing.loadFromData(render_png(mol, width, height, scale=ratio, dark=dark), "PNG")
    drawing.setDevicePixelRatio(ratio)
    return drawing


class TemplateLibrary(QDialog):
    """Groups on the left, drawings on the right; double-click (or OK) picks one."""

    chosen = Signal(str)  # the entry's name

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Template Library")
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find a template")
        self.search.setClearButtonEnabled(True)
        self.groups = QListWidget()
        self.groups.addItems([ALL, *library.GROUPS])
        self.groups.setFixedWidth(190)
        self.groups.setCurrentRow(0)
        self.tiles = QListWidget()
        self.tiles.setViewMode(QListView.ViewMode.IconMode)
        self.tiles.setResizeMode(QListView.ResizeMode.Adjust)
        self.tiles.setMovement(QListView.Movement.Static)
        self.tiles.setIconSize(TILE)
        self.tiles.setGridSize(QSize(TILE.width() + 20, TILE.height() + 36))
        self.tiles.setWordWrap(True)
        self.tiles.setUniformItemSizes(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        body = QHBoxLayout()
        body.addWidget(self.groups)
        body.addWidget(self.tiles, 1)
        layout = QVBoxLayout(self)
        layout.addWidget(self.search)
        layout.addLayout(body, 1)
        layout.addWidget(self.buttons)
        self.resize(760, 520)

        self.search.textChanged.connect(lambda _t: self.fill())
        self.groups.currentRowChanged.connect(lambda _r: self.fill())
        self.tiles.itemDoubleClicked.connect(lambda _item: self.accept())
        self.tiles.currentItemChanged.connect(lambda *_: self._update_ok())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.fill()

    def fill(self) -> None:
        group = self.groups.currentItem().text() if self.groups.currentItem() else ALL
        entries = library.search(self.search.text())
        if group != ALL and not self.search.text().strip():  # a search looks in every group
            entries = [e for e in entries if e.group == group]
        dark = self.palette().base().color().lightness() < 128
        ratio = self.devicePixelRatioF()
        self.tiles.clear()
        for entry in entries:
            item = QListWidgetItem(QIcon(_drawing(entry.smiles, dark, ratio)), entry.name)
            item.setData(Qt.ItemDataRole.UserRole, entry.name)
            item.setToolTip(f"{entry.name} · {entry.group}")
            item.setTextAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
            self.tiles.addItem(item)
        if self.tiles.count():
            self.tiles.setCurrentRow(0)
        self._update_ok()

    def names(self) -> list[str]:
        return [
            self.tiles.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.tiles.count())
        ]

    def current(self) -> str | None:
        item = self.tiles.currentItem()
        return item.data(Qt.ItemDataRole.UserRole) if item is not None else None

    def _update_ok(self) -> None:
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            self.current() is not None
        )

    def accept(self) -> None:
        name = self.current()
        if name is None:
            return
        super().accept()
        self.chosen.emit(name)
