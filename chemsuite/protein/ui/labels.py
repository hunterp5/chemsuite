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
"""Residue and atom labels over the protein view, kept legible.

:func:`layout` chooses which labels to draw: front to back, skipping one that would overlap a
label already placed, and at most :data:`MAX_LABELS`, so a crowded view thins out instead of
turning into a smear. Labels behind the camera, off screen or outside the slab are dropped.
:func:`paint` draws each with a halo in the window's background colour, so it reads over any
atom, in light and dark.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QFont, QFontMetricsF, QPainter, QPainterPath, QPalette, QPen

#: the most labels drawn at once, and the most tried (front first) to find them
MAX_LABELS = 300
MAX_TRIED = 3000
#: pixels between a label's point and its text
LIFT = 6.0
HALO = 3.0


def layout(
    boxes: np.ndarray, depth: np.ndarray, near: float = 0.0, far: float = np.inf
) -> np.ndarray:
    """Which labels to draw, front first: ``boxes`` (n, 4) as x, y, width, height in pixels,
    ``depth`` their view depth; ``near`` and ``far`` the slab kept."""
    boxes = np.asarray(boxes, np.float64).reshape(-1, 4)
    depth = np.asarray(depth, np.float64)
    inside = np.flatnonzero((depth > max(near, 0.0)) & (depth < far))
    order = inside[np.argsort(depth[inside], kind="stable")][:MAX_TRIED]
    placed = np.zeros((MAX_LABELS, 4))
    out: list[int] = []
    for i in order.tolist():
        x, y, w, h = boxes[i]
        p = placed[: len(out)]
        overlaps = (
            (p[:, 0] < x + w)
            & (x < p[:, 0] + p[:, 2])
            & (p[:, 1] < y + h)
            & (y < p[:, 1] + p[:, 3])
        )
        if overlaps.any():
            continue
        placed[len(out)] = (x, y, w, h)
        out.append(i)
        if len(out) == MAX_LABELS:
            break
    return np.array(out, np.int64)


def boxes(xy: np.ndarray, texts: list[str], font: QFont) -> np.ndarray:
    """Each label's box: its text centred a little above its point."""
    metrics = QFontMetricsF(font)
    out = np.zeros((len(texts), 4))
    for i, text in enumerate(texts):
        width, height = metrics.horizontalAdvance(text), metrics.height()
        out[i] = (xy[i, 0] - width / 2, xy[i, 1] - LIFT - height, width, height)
    return out


def paint(painter: QPainter, places: np.ndarray, texts: list[str], palette: QPalette) -> None:
    """Draw the texts in the boxes ``places`` (n, 4), with a halo."""
    font = painter.font()
    ascent = QFontMetricsF(font).ascent()
    halo = QPen(palette.base().color(), HALO)
    halo.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    ink = palette.text().color()
    for (x, y, _w, _h), text in zip(places, texts, strict=True):
        path = QPainterPath()
        path.addText(QPointF(x, y + ascent), font, text)
        painter.strokePath(path, halo)
        painter.fillPath(path, ink)


def label_font(base: QFont) -> QFont:
    font = QFont(base)
    font.setWeight(QFont.Weight.DemiBold)
    return font
