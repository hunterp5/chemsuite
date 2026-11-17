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
"""The Structure Viewer's 2D page: one structure drawn large, redrawn at its own size."""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QImage, QPainter, QPixmap
from PySide6.QtWidgets import QLabel, QSizePolicy
from rdkit import Chem

from chemsuite.core.depict import render_png

SCALE = 1.6  # bonds drawn larger than in table cells (shrunk to fit if needed)
RESIZE_MS = 60


def structure_image(mol: Chem.Mol, width: int, height: int, dpr: float, dark: bool) -> QImage:
    """*mol* drawn in a ``width x height`` (device-independent) box at *dpr*."""
    png = render_png(mol, int(width * dpr), int(height * dpr), scale=dpr * SCALE, dark=dark)
    image = QImage.fromData(png, "PNG")
    image.setDevicePixelRatio(dpr)
    return image


def flat_picture(mol: Chem.Mol | None, width: int, height: int, dpr: float, base: QColor) -> QImage:
    """The 2D page as an image: *mol* (if any) drawn on the *base* color."""
    image = QImage(int(width * dpr), int(height * dpr), QImage.Format.Format_ARGB32)
    image.setDevicePixelRatio(dpr)
    image.fill(base)
    if mol is not None:
        drawing = structure_image(mol, width, height, dpr, base.lightness() < 128)
        painter = QPainter(image)
        painter.drawImage(0, 0, drawing)
        painter.end()
    return image


class Depiction2D(QLabel):
    """A structure, or a message when there is none; redrawn once resizing settles."""

    def __init__(self, on_resized: Callable[[], None]) -> None:
        super().__init__()
        self.setObjectName("ViewerDepiction")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumSize(160, 140)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.mol: Chem.Mol | None = None
        self._image = QImage()  # the last drawing, scaled to fit until it is redrawn
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(RESIZE_MS)
        self._timer.timeout.connect(on_resized)

    def show_message(self, text: str) -> None:
        self.mol = None
        self._image = QImage()
        self.setPixmap(QPixmap())
        self.setText(text)

    def show_mol(self, mol: Chem.Mol) -> None:
        self.mol = mol
        dark = self.palette().base().color().lightness() < 128
        self._image = structure_image(
            mol, self.width(), self.height(), self.devicePixelRatioF(), dark
        )
        self.setText("")
        self.setPixmap(QPixmap.fromImage(self._image))

    def resizeEvent(self, event) -> None:
        """Fit the last drawing at once (never cropped), and redraw it sharp once resizing
        settles."""
        super().resizeEvent(event)
        if not self._image.isNull():
            dpr = self.devicePixelRatioF()
            fitted = self._image.scaled(
                int(self.width() * dpr),
                int(self.height() * dpr),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            fitted.setDevicePixelRatio(dpr)
            self.setPixmap(QPixmap.fromImage(fitted))
        self._timer.start()
