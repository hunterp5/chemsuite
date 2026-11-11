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
"""The sequence strip: one row of one-letter codes per polymer chain, under the 3D view.

Positions follow the full (SEQRES) sequence, so residues missing from the model show faint.
Selected residues are highlighted and the hovered one framed. Click selects a residue, drag or
Shift-click a range, Ctrl-click (Cmd on macOS) toggles; hovering a letter lights the residue up
in 3D. Only the visible part is painted, so long chains stay fast.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import QPoint, QRect, QSize, Qt, Signal
from PySide6.QtGui import (
    QFontDatabase,
    QKeyEvent,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QPalette,
)
from PySide6.QtWidgets import QScrollArea, QWidget

from chemsuite.protein.structure.sequence import ChainSequence, chain_sequences

from .session import ADD, SET, TOGGLE, ProteinSession

_PAD = 4
_EVERY = 10  # a residue number above every tenth letter


@dataclass(frozen=True)
class _Row:
    protein: str
    label: str
    sequence: ChainSequence
    letters: str
    numbers: np.ndarray  # author residue number at each modelled position


class SequenceStrip(QWidget):
    """The painted rows; lives in a :class:`SequenceView` (a scroll area)."""

    residueHovered = Signal(str, int)  # protein, residue (-1: none)
    cursorMoved = Signal(QRect)  # where the keyboard's cursor is now (to scroll to it)

    def __init__(self, session: ProteinSession, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session = session
        self.setMouseTracking(True)
        self.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.rows: list[_Row] = []
        self._cache: dict[int, list[ChainSequence]] = {}
        self._hover: tuple[int, int] | None = None  # row, position
        self._anchor: tuple[int, int] | None = None
        self._cursor: tuple[int, int] | None = None  # the keyboard's row and position
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._press_mode = SET
        session.proteinsChanged.connect(self.reload)
        session.selectionChanged.connect(self.update)
        self.reload()

    # layout
    @property
    def cell(self) -> int:
        return self.fontMetrics().horizontalAdvance("W") + 2

    @property
    def line(self) -> int:
        return self.fontMetrics().height() + 2

    @property
    def row_height(self) -> int:
        return 2 * self.line + _PAD

    @property
    def label_width(self) -> int:
        widths = [self.fontMetrics().horizontalAdvance(r.label) for r in self.rows]
        return max(widths, default=0) + 3 * _PAD

    def reload(self) -> None:
        rows = []
        for shown in self.session.displays:
            if not shown.visible:
                continue
            structure = shown.structure
            key = id(structure)
            if key not in self._cache:
                self._cache[key] = chain_sequences(structure)
            for sequence in self._cache[key]:
                chain = structure.chain_name[sequence.chain]
                numbers = np.where(
                    sequence.modelled, structure.res_seq[np.maximum(sequence.residue_at, 0)], 0
                )
                rows.append(_Row(shown.id, f"{shown.name} {chain}", sequence,
                                 sequence.letters, numbers))  # fmt: skip
        self.rows = rows
        self._hover = None
        self._cursor = self._anchor = None
        self.updateGeometry()
        self.adjustSize()
        self.update()

    def sizeHint(self) -> QSize:
        longest = max((len(r.letters) for r in self.rows), default=0)
        return QSize(self.label_width + longest * self.cell + _PAD,
                     max(len(self.rows), 1) * self.row_height + _PAD)  # fmt: skip

    def position_at(self, point: QPoint) -> tuple[int, int] | None:
        row = point.y() // self.row_height
        position = (point.x() - self.label_width) // self.cell
        if 0 <= row < len(self.rows) and 0 <= position < len(self.rows[row].letters):
            return int(row), int(position)
        return None

    def residue_at(self, row: int, position: int) -> int:
        return int(self.rows[row].sequence.residue_at[position])

    # painting
    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        palette = self.palette()
        text = palette.color(QPalette.ColorRole.Text)
        faint = palette.color(QPalette.ColorRole.PlaceholderText)
        selected_fill = palette.color(QPalette.ColorRole.Highlight)
        selected_text = palette.color(QPalette.ColorRole.HighlightedText)
        cell, line, label_width = self.cell, self.line, self.label_width
        area = event.rect()
        first = max((area.left() - label_width) // cell, 0)
        for r, row in enumerate(self.rows):
            top = r * self.row_height
            if top > area.bottom() or top + self.row_height < area.top():
                continue
            painter.setPen(faint)
            painter.drawText(QRect(_PAD, top + line, label_width, line),
                             Qt.AlignmentFlag.AlignVCenter, row.label)  # fmt: skip
            shown = self.session.display(row.protein)
            chosen = (
                np.zeros(0, np.int64)
                if shown is None
                else np.unique(shown.structure.residue_of[shown.selected])
            )
            last = min(len(row.letters), (area.right() - label_width) // cell + 2)
            residue_at = row.sequence.residue_at
            for k in range(first, last):
                x = label_width + k * cell
                box = QRect(x, top + line, cell, line)
                residue = int(residue_at[k])
                if residue >= 0 and np.isin(residue, chosen):
                    painter.fillRect(box, selected_fill)
                    painter.setPen(selected_text)
                else:
                    painter.setPen(text if residue >= 0 else faint)
                painter.drawText(box, Qt.AlignmentFlag.AlignCenter, row.letters[k])
                if k % _EVERY == 0 and residue >= 0:
                    painter.setPen(faint)
                    number = str(int(row.numbers[k]))
                    painter.drawText(QRect(x, top, cell * _EVERY, line),
                                     Qt.AlignmentFlag.AlignLeft, number)  # fmt: skip
            marked = self._hover if self._hover is not None else self._focused_cursor()
            if marked is not None and marked[0] == r:
                painter.setPen(text)
                painter.drawRect(QRect(label_width + marked[1] * cell, top + line,
                                       cell - 1, line - 1))  # fmt: skip
        painter.end()

    # mouse
    def mousePressEvent(self, event: QMouseEvent) -> None:
        hit = self.position_at(event.position().toPoint())
        if hit is None:
            return
        modifiers = event.modifiers()
        if modifiers & Qt.KeyboardModifier.ShiftModifier and self._anchor is not None:
            self._select_range(self._anchor, hit, ADD)
            return
        self._press_mode = TOGGLE if modifiers & Qt.KeyboardModifier.ControlModifier else SET
        self._anchor = self._cursor = hit
        self._select_range(hit, hit, self._press_mode)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        hit = self.position_at(event.position().toPoint())
        dragging = event.buttons() & Qt.MouseButton.LeftButton and hit and self._anchor
        if dragging and hit[0] == self._anchor[0] and self._press_mode == SET:
            self._select_range(self._anchor, hit, SET)
        if hit != self._hover:
            self._hover = hit
            self.update()
            if hit is None:
                self.residueHovered.emit("", -1)
            else:
                self.residueHovered.emit(self.rows[hit[0]].protein, self.residue_at(*hit))

    def leaveEvent(self, event) -> None:
        super().leaveEvent(event)
        if self._hover is not None:
            self._hover = None
            self.update()
            self.residueHovered.emit("", -1)

    # keyboard: arrows move a cursor, Shift+arrows select from where it started, Return
    # selects the residue under it
    def _focused_cursor(self) -> tuple[int, int] | None:
        return self._cursor if self.hasFocus() else None

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if not self.rows:
            super().keyPressEvent(event)
            return
        row, position = self._cursor or (0, 0)
        key = event.key()
        steps = {Qt.Key.Key_Left: (0, -1), Qt.Key.Key_Right: (0, 1),
                 Qt.Key.Key_Up: (-1, 0), Qt.Key.Key_Down: (1, 0)}  # fmt: skip
        if key in steps:
            row = min(max(row + steps[key][0], 0), len(self.rows) - 1)
            position += steps[key][1]
        elif key == Qt.Key.Key_Home:
            position = 0
        elif key == Qt.Key.Key_End:
            position = len(self.rows[row].letters) - 1
        elif key in (Qt.Key.Key_Return, Qt.Key.Key_Enter) and self._cursor is not None:
            self._anchor = self._cursor
            self._select_range(self._cursor, self._cursor, SET)
            return
        else:
            super().keyPressEvent(event)
            return
        position = min(max(position, 0), len(self.rows[row].letters) - 1)
        extend = event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        if not extend or self._anchor is None or self._anchor[0] != row:
            self._anchor = (row, position)
        self._cursor = (row, position)
        if extend:
            self._select_range(self._anchor, self._cursor, SET)
        self.update()
        self.residueHovered.emit(self.rows[row].protein, self.residue_at(row, position))
        self.cursorMoved.emit(QRect(self.label_width + position * self.cell,
                                    row * self.row_height, self.cell, self.row_height))  # fmt: skip

    def focusOutEvent(self, event) -> None:
        super().focusOutEvent(event)
        self.update()
        if self._hover is None:
            self.residueHovered.emit("", -1)

    def _select_range(self, start: tuple[int, int], end: tuple[int, int], mode: str) -> None:
        row = self.rows[start[0]]
        if end[0] != start[0]:
            end = (start[0], len(row.letters) - 1 if end[0] > start[0] else 0)
        low, high = sorted((start[1], end[1]))
        residues = row.sequence.residue_at[low : high + 1]
        residues = residues[residues >= 0]
        shown = self.session.display(row.protein)
        if shown is None:
            return
        mask = np.isin(shown.structure.residue_of, residues) & shown.active
        self.session.select({row.protein: mask}, mode)

    def show_residue(self, protein: str, residue: int) -> QRect | None:
        """Where a residue's letter is (to scroll to it)."""
        for r, row in enumerate(self.rows):
            if row.protein == protein:
                found = np.flatnonzero(row.sequence.residue_at == residue)
                if len(found):
                    x = self.label_width + int(found[0]) * self.cell
                    return QRect(x, r * self.row_height, self.cell, self.row_height)
        return None


class SequenceView(QScrollArea):
    def __init__(self, session: ProteinSession, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.strip = SequenceStrip(session)
        self.setWidget(self.strip)
        self.setWidgetResizable(False)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setFrameShape(QScrollArea.Shape.NoFrame)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)  # the strip inside takes the keys
        self.strip.cursorMoved.connect(
            lambda box: self.ensureVisible(box.center().x(), box.center().y(), 60, 0)
        )

    def sizeHint(self) -> QSize:
        hint = self.strip.sizeHint()
        return QSize(400, min(hint.height(), 4 * self.strip.row_height) + 16)

    def scroll_to(self, protein: str, residue: int) -> None:
        box = self.strip.show_residue(protein, residue)
        if box is not None:
            self.ensureVisible(box.center().x(), box.center().y(), 60, 0)
