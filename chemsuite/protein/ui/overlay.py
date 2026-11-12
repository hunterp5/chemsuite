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
"""What the protein view draws on top: residue and atom labels (:mod:`.labels`),
measurements (a dashed line and its value) and the atoms picked for the measurement being
made."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget

from . import labels as labelling
from .session import SessionSignals

MEASURE = QColor(230, 70, 20)
PICK_RADIUS = 9.0  # pixels


class ProteinOverlay(QWidget):
    def __init__(self, view, signals: SessionSignals) -> None:
        super().__init__()
        self.view = view
        self.session = signals.session
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        view.cameraMoved.connect(self.update)
        signals.measurementsChanged.connect(self.update)
        signals.sceneChanged.connect(self._labels_moved)
        signals.labelsChanged.connect(self._labels_moved)
        self._labelled: tuple[np.ndarray, list[str]] | None = None

    def _labels_moved(self) -> None:
        self._labelled = None
        self.update()

    def placed_labels(self) -> tuple[np.ndarray, list[str]]:
        """The labels drawn now: their boxes (n, 4) and texts (for tests and painting)."""
        if self._labelled is None:
            self._labelled = self.session.labelled()
        positions, texts = self._labelled
        if not texts:
            return np.zeros((0, 4)), []
        xy, depth = self.view.camera.project(positions, self.width(), self.height())
        font = labelling.label_font(self.font())
        boxes = labelling.boxes(xy, texts, font)
        on_screen = (boxes[:, 0] + boxes[:, 2] > 0) & (boxes[:, 0] < self.width()) & (
            boxes[:, 1] + boxes[:, 3] > 0) & (boxes[:, 1] < self.height())  # fmt: skip
        depth = np.where(on_screen, depth, -1.0)
        near, far = self.view.slab
        chosen = labelling.layout(boxes, depth, near, far)
        return boxes[chosen], [texts[i] for i in chosen.tolist()]

    def labels(self) -> list[tuple[QPointF, str]]:
        """Where each measurement's value goes, and the value (for tests and painting)."""
        out = []
        camera = self.view.camera
        for points, text in self.session.measured():
            xy, depth = camera.project(np.array(points), self.width(), self.height())
            if (depth > 0).all():
                out.append((QPointF(*xy.mean(axis=0)), text))
        return out

    def paintEvent(self, _event) -> None:
        camera = self.view.camera
        width, height = self.width(), self.height()
        measured = self.session.measured()
        picked = [p for p in (self.session.where(q) for q in self.session.picked) if p]
        boxes, texts = self.placed_labels()
        if not measured and not picked and not texts:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        if texts:
            painter.setFont(labelling.label_font(self.font()))
            labelling.paint(painter, boxes, texts, self.palette())
        for points, text in measured:
            xy, depth = camera.project(np.array(points), width, height)
            if (depth <= 0).any():
                continue
            painter.setPen(QPen(MEASURE, 1.5, Qt.PenStyle.DashLine))
            path = QPainterPath(QPointF(*xy[0]))
            for x, y in xy[1:]:
                path.lineTo(QPointF(x, y))
            painter.drawPath(path)
            self._label(painter, text, QPointF(*xy.mean(axis=0)))
        if picked:
            xy, depth = camera.project(np.array([p[1] for p in picked]), width, height)
            painter.setPen(QPen(MEASURE, 2.0))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            for (x, y), d in zip(xy, depth, strict=True):
                if d > 0:
                    painter.drawEllipse(QPointF(x, y), PICK_RADIUS, PICK_RADIUS)
        painter.end()

    def _label(self, painter: QPainter, text: str, at: QPointF) -> None:
        font = QFont(self.font())
        font.setBold(True)
        painter.setFont(font)
        box = QRectF(painter.fontMetrics().boundingRect(text)).adjusted(-5, -2, 5, 2)
        box.moveCenter(at - QPointF(0, box.height()))
        fill = QColor(self.palette().base().color())
        fill.setAlphaF(0.85)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(fill)
        painter.drawRoundedRect(box, 4, 4)
        painter.setPen(MEASURE)
        painter.drawText(box, Qt.AlignmentFlag.AlignCenter, text)
