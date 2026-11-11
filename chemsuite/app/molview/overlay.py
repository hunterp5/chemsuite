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
"""What the 3D Viewer draws on top of the molecules: the atom under the pointer, the atoms
picked for a measurement, and the measurement itself."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget

from chemsuite.molview.measure import describe

HIGHLIGHT = QColor(230, 70, 20)
RING_PAD = 3.0  # pixels outside the atom


class MolOverlay(QWidget):
    def __init__(self, renderer) -> None:
        super().__init__()
        self.renderer = renderer
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        self.hover = -1
        self.marked: list[int] = []

    def set_hover(self, atom: int) -> None:
        if atom != self.hover:
            self.hover = atom
            self.update()

    def set_marked(self, atoms: list[int]) -> None:
        self.marked = list(atoms)
        self.update()

    def measurement(self) -> str:
        scene = self.renderer.scene
        if not 2 <= len(self.marked) <= 4 or max(self.marked) >= len(scene):
            return ""
        return describe([scene.positions[i] for i in self.marked])

    def paintEvent(self, _event) -> None:
        scene = self.renderer.scene
        atoms = [a for a in (*self.marked, self.hover) if 0 <= a < len(scene)]
        if not atoms:
            return
        camera = self.renderer.camera
        width, height = self.width(), self.height()
        xy, depth = camera.project(scene.positions[atoms], width, height)
        rings = {
            atom: (QPointF(*xy[k]), scene.radii[atom] * camera.pixels_per_unit(depth[k], height))
            for k, atom in enumerate(atoms)
            if depth[k] > 0
        }
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        marked = [a for a in self.marked if a in rings]
        if len(marked) > 1:
            dashed = QPen(HIGHLIGHT, 1.5, Qt.PenStyle.DashLine)
            painter.setPen(dashed)
            path = QPainterPath(rings[marked[0]][0])
            for atom in marked[1:]:
                path.lineTo(rings[atom][0])
            painter.drawPath(path)
        for atom, (centre, radius) in rings.items():
            width_px = 2.0 if atom in self.marked else 1.5
            colour = QColor(HIGHLIGHT)
            if atom not in self.marked:
                colour.setAlphaF(0.7)
            painter.setPen(QPen(colour, width_px))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawEllipse(centre, radius + RING_PAD, radius + RING_PAD)
        text = self.measurement()
        if text and len(marked) == len(self.marked):
            x = float(np.mean([rings[a][0].x() for a in marked]))
            top = min(rings[a][0].y() - rings[a][1] - RING_PAD for a in marked)
            self._label(painter, text, x, top)
        painter.end()

    def _label(self, painter: QPainter, text: str, x: float, bottom: float) -> None:
        """*text* centred on *x*, just above *bottom*: clear of the atoms and lines it names."""
        font = QFont(self.font())
        font.setBold(True)
        painter.setFont(font)
        metrics = painter.fontMetrics()
        box = QRectF(metrics.boundingRect(text)).adjusted(-5, -2, 5, 2)
        box.moveCenter(QPointF(x, bottom - 4 - box.height() / 2))
        box.moveTop(max(box.top(), 2.0))
        fill = QColor(self.palette().base().color())
        fill.setAlphaF(0.85)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(fill)
        painter.drawRoundedRect(box, 4, 4)
        painter.setPen(HIGHLIGHT)
        painter.drawText(box, Qt.AlignmentFlag.AlignCenter, text)
