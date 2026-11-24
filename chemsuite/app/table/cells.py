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

"""How a table cell is drawn: structures as pictures, everything else as Qt's text.

:class:`StructureImages` renders and caches the pictures (one cache per window, within a byte
budget), :class:`CellDelegate` draws them and answers the hover preview, and :class:`Alignment`
is a column's chosen layout core. What a drawing *marks* is in :mod:`.highlights`.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import QEvent, QModelIndex, QRect, QSize, Qt
from PySide6.QtGui import QColor, QHelpEvent, QImage, QPainter, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QToolTip,
    QWidget,
)
from rdkit import Chem

from chemsuite.core.depict import aligned_mol, render_png, render_reaction_png
from chemsuite.core.table.columns import REACTION, STRUCTURE, ParsedColumn

from .highlights import AtomHighlight, AtomSource, Highlight, marked, resolve
from .structure_preview import preview_html
from .table_model import KIND_ROLE, STATUS_ROLE, STRUCTURE_ROLE

#: Column kinds whose cells are drawn as pictures.
DRAWN = (STRUCTURE, REACTION)


@dataclass(frozen=True)
class Alignment:
    """Drawings in a structure column laid out to match *core* (display only)."""

    text: str  # the core as the user gave it (SMILES or SMARTS); also its cache key
    core: Chem.Mol


class StructureImages:
    """Rendered structure pictures, keyed by molecule, size and highlight, within a byte budget."""

    def __init__(self, max_bytes: int = 192 * 1024**2) -> None:
        self._cache: OrderedDict[tuple, QImage] = OrderedDict()
        self._bytes = 0
        self.max_bytes = max_bytes

    def image(
        self,
        column: ParsedColumn,
        slot: int,
        size: QSize,
        dpr: float,
        dark: bool,
        highlight: Highlight | AtomHighlight | None = None,
        alignment: Alignment | None = None,
    ):
        mol = column.obj(slot)  # parses a pending cell on first use
        if mol is None:
            return None
        blob = column.blobs(np.array([slot]))[0]
        key = (
            blob,
            size.width(),
            size.height(),
            round(dpr, 2),
            dark,
            highlight.key if highlight else None,
            alignment.text if alignment and column.kind != REACTION else None,
        )
        image = self._cache.get(key)
        if image is not None:
            self._cache.move_to_end(key)
            return image
        width, height = int(size.width() * dpr), int(size.height() * dpr)
        if column.kind == REACTION:
            png = render_reaction_png(mol, width, height, scale=dpr, dark=dark)
            if png is None:
                return None
        else:
            atoms: list[int] = []
            bonds: list[int] = []
            colors: dict[int, tuple[float, float, float]] = {}
            if alignment is not None:
                mol = aligned_mol(mol, alignment.core)
            if highlight is not None:
                mol, atoms, bonds, colors = marked(mol, highlight, dark)
            png = render_png(
                mol,
                width,
                height,
                scale=dpr,
                dark=dark,
                highlight_atoms=atoms,
                highlight_bonds=bonds,
                highlight_colors=colors,
            )
        image = QImage.fromData(png, "PNG")
        image.setDevicePixelRatio(dpr)
        self._cache[key] = image
        self._bytes += image.sizeInBytes()
        while self._bytes > self.max_bytes and self._cache:
            _, old = self._cache.popitem(last=False)
            self._bytes -= old.sizeInBytes()
        return image

    def clear(self) -> None:
        self._cache.clear()
        self._bytes = 0


class CellDelegate(QStyledItemDelegate):
    """Draws structure cells as pictures; every other cell is drawn by Qt as text."""

    def __init__(self, images: StructureImages, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._images = images
        #: column name -> what to highlight in that column's drawings
        self.highlights: dict[str, Highlight | AtomSource] = {}
        #: column name -> the core its drawings are aligned to
        self.alignments: dict[str, Alignment] = {}

    def paint(self, painter: QPainter, option: QStyleOptionViewItem, index: QModelIndex) -> None:
        if index.data(KIND_ROLE) not in DRAWN:
            super().paint(painter, option, index)
            return
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        opt.text = ""
        widget = option.widget
        style = widget.style() if widget is not None else QApplication.style()
        style.drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, painter, widget)

        status = index.data(STATUS_ROLE)
        rect = option.rect.adjusted(3, 3, -3, -3)
        if status == "blank" or rect.width() < 8 or rect.height() < 8:
            return
        column, slot = index.data(STRUCTURE_ROLE)
        if status != "invalid":
            dark = option.palette.color(QPalette.ColorRole.Base).lightness() < 128
            dpr = painter.device().devicePixelRatioF() if painter.device() else 1.0
            highlight = self.highlight_for(index, slot)
            alignment = self.alignment_for(index)
            image = self._images.image(column, slot, rect.size(), dpr, dark, highlight, alignment)
            if image is not None:
                painter.drawImage(rect.topLeft(), image)
                return
        self._paint_invalid(painter, rect, column.display(slot))

    def highlight_for(self, index: QModelIndex, slot: int) -> Highlight | AtomHighlight | None:
        """What to mark in the drawing of *index*'s cell, resolved for its row."""
        if not self.highlights:
            return None
        model = index.model()
        name = model.column_name(index.column())
        highlight = self.highlights.get(name)
        if isinstance(highlight, AtomSource):
            return resolve(model.table, name, highlight, slot)
        return highlight

    def alignment_for(self, index: QModelIndex) -> Alignment | None:
        """The core *index*'s column lays its drawings out on, if it has one."""
        if not self.alignments:
            return None
        return self.alignments.get(index.model().column_name(index.column()))

    def helpEvent(
        self,
        event: QHelpEvent,
        view: QAbstractItemView,
        option: QStyleOptionViewItem,
        index: QModelIndex,
    ) -> bool:
        """The pointer resting on a structure cell enlarges the drawing the cell shows: its
        marks and its alignment with it.

        A cell with neither falls through to the model's ``Qt.ToolTipRole``, which is where
        the plain preview comes from.
        """
        if (
            event.type() != QEvent.Type.ToolTip
            or index.data(KIND_ROLE) != STRUCTURE
            or index.data(STATUS_ROLE) != "ok"
        ):
            return super().helpEvent(event, view, option, index)
        column, slot = index.data(STRUCTURE_ROLE)
        highlight = self.highlight_for(index, slot)
        alignment = self.alignment_for(index)
        if highlight is None and alignment is None:
            return super().helpEvent(event, view, option, index)
        mol = column.obj(slot)
        if mol is None or not mol.GetNumAtoms():
            return super().helpEvent(event, view, option, index)
        QToolTip.showText(
            event.globalPos(),
            preview_html(mol, column.display(slot), highlight, alignment),
            view,
        )
        return True

    @staticmethod
    def _paint_invalid(painter: QPainter, rect: QRect, text: str) -> None:
        painter.save()
        painter.setPen(QColor(200, 40, 40))
        first_line = next((line for line in text.splitlines() if line.strip()), text)
        painter.drawText(
            rect,
            Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap,
            f"Invalid structure\n{first_line[:60]}",
        )
        painter.restore()
