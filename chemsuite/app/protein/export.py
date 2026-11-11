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
"""File ▸ Export Image… and Edit ▸ Copy Image: the view as a picture, for slides and papers.

The picture keeps the view's shape (so labels and measurements land where they are on screen)
and is drawn at the size asked for, not scaled up from the screen: in pixels, or in centimetres
or inches at a resolution, which is also written into the file. The background can be left
see-through.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QGuiApplication, QImage, QPainter, QRegion
from PySide6.QtWidgets import QCheckBox, QDoubleSpinBox, QFileDialog, QSpinBox, QWidget

from chemsuite.app.dialogs import Segmented, ToolDialog

if TYPE_CHECKING:
    from .panel import ProteinPanel

PIXELS, CENTIMETRES, INCHES = "px", "cm", "in"
UNITS = ((PIXELS, "Pixels"), (CENTIMETRES, "cm"), (INCHES, "Inches"))
CM_PER_INCH = 2.54
MAX_PIXELS = 12_000  # a side; GPUs draw at least 16k, and a bigger PNG is no use
COPY_SCALE = 2  # Copy Image draws at twice the view's size
FORMATS = "PNG (*.png);;TIFF (*.tif *.tiff);;JPEG (*.jpg *.jpeg)"


def picture(panel: ProteinPanel, width: int, height: int, transparent: bool = False) -> QImage:
    """The view drawn at *width* x *height*, with its labels and measurements on top (null
    when the view cannot draw)."""
    view, overlay = panel.view, panel.overlay
    if view is None or overlay is None:
        return QImage()
    image = view.picture(width, height, transparent)
    if image.isNull():
        return image
    image.setDevicePixelRatio(1.0)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
    painter.scale(width / max(view.width(), 1), height / max(view.height(), 1))
    overlay.render(painter, QPoint(), QRegion(), QWidget.RenderFlag.DrawChildren)
    painter.end()
    return image


def copy_image(panel: ProteinPanel) -> bool:
    """Put the view on the clipboard at twice its size; False when it cannot be drawn."""
    if panel.view is None:
        return False
    width, height = panel.view.width() * COPY_SCALE, panel.view.height() * COPY_SCALE
    image = picture(panel, width, height)
    if image.isNull():
        panel.status.setText("The view could not be drawn")
        return False
    QGuiApplication.clipboard().setImage(image)
    panel.status.setText(f"Copied the view ({width:,} × {height:,} px)")
    return True


def save_image(image: QImage, path: str, dpi: int) -> bool:
    """Write *image* with its resolution recorded (JPEG has no transparency: white behind)."""
    dots_per_metre = round(dpi / CM_PER_INCH * 100)
    image.setDotsPerMeterX(dots_per_metre)
    image.setDotsPerMeterY(dots_per_metre)
    if Path(path).suffix.lower() in (".jpg", ".jpeg") and image.hasAlphaChannel():
        flat = QImage(image.size(), QImage.Format.Format_RGB32)
        flat.fill(Qt.GlobalColor.white)
        painter = QPainter(flat)
        painter.drawImage(0, 0, image)
        painter.end()
        flat.setDotsPerMeterX(dots_per_metre)
        flat.setDotsPerMeterY(dots_per_metre)
        image = flat
    return image.save(path)


class ExportImageDialog(ToolDialog):
    def __init__(self, panel: ProteinPanel) -> None:
        super().__init__(panel.window(), "Export Image", guide="proteins")
        self.panel = panel
        view = panel.view
        self.aspect = view.width() / max(view.height(), 1) if view is not None else 4 / 3
        self.path = ""
        size = self.section("Size")
        self.units = Segmented(UNITS, PIXELS)
        size.row("Units", self.units, help="Size the picture in pixels, or in print units")
        self.width = _size_field()
        self.height = _size_field()
        size.row("Width", self.width)
        size.row("Height", self.height, help="Follows the width: the picture keeps the view's "
                 "shape, so labels and measurements stay where they are")  # fmt: skip
        self.dpi = QSpinBox()
        self.dpi.setRange(72, 2400)
        self.dpi.setValue(300)
        self.dpi.setSuffix(" dpi")
        size.row("Resolution", self.dpi, help="Dots per inch: turns cm or inches into pixels, "
                 "and is written into the file for page layout")  # fmt: skip
        output = self.section("Output")
        self.transparent = QCheckBox("See-through background")
        output.row("Background", self.transparent,
                   help="Leave the background out (PNG and TIFF), to lay the picture on a "
                   "coloured slide")  # fmt: skip
        self._linking = False
        start = view.width() * COPY_SCALE if view is not None else 1600
        self.width.setValue(start)
        self._follow_width()
        self.width.valueChanged.connect(lambda _v: self._follow_width())
        self.height.valueChanged.connect(lambda _v: self._follow_height())
        self.units.changed.connect(self._convert)
        self.dpi.valueChanged.connect(lambda _v: self.refresh())
        self._unit = PIXELS
        self.refresh()

    # sizes
    def _to_pixels(self, value: float, unit: str | None = None) -> int:
        unit = unit or self.units.value()
        if unit == CENTIMETRES:
            return round(value / CM_PER_INCH * self.dpi.value())
        if unit == INCHES:
            return round(value * self.dpi.value())
        return round(value)

    def _from_pixels(self, pixels: float, unit: str) -> float:
        if unit == CENTIMETRES:
            return pixels / self.dpi.value() * CM_PER_INCH
        if unit == INCHES:
            return pixels / self.dpi.value()
        return pixels

    def pixel_size(self) -> tuple[int, int]:
        width = max(self._to_pixels(self.width.value()), 1)
        return width, max(round(width / self.aspect), 1)

    def _follow_width(self) -> None:
        if self._linking:
            return
        self._linking = True
        self.height.setValue(self.width.value() / self.aspect)
        self._linking = False
        self.refresh()

    def _follow_height(self) -> None:
        if self._linking:
            return
        self._linking = True
        self.width.setValue(self.height.value() * self.aspect)
        self._linking = False
        self.refresh()

    def _convert(self, unit: str) -> None:
        """Keep the picture's size when the units change."""
        pixels = self._to_pixels(self.width.value(), self._unit)
        self._unit = unit
        decimals = 0 if unit == PIXELS else 2
        for field in (self.width, self.height):
            field.setDecimals(decimals)
        self.width.setValue(self._from_pixels(pixels, unit))
        self.refresh()

    # the kit
    def problem(self) -> str:
        if self.panel.view is None:
            return "3D drawing is not available on this computer"
        width, height = self.pixel_size()
        if max(width, height) > MAX_PIXELS:
            return f"At most {MAX_PIXELS:,} pixels a side ({width:,} × {height:,} asked for)"
        return ""

    def verb(self) -> str:
        width, height = self.pixel_size()
        return f"Export {width:,} × {height:,} px"

    def run(self):
        path = self.path
        if not path:
            path, _ = QFileDialog.getSaveFileName(self, "Export Image", "protein.png", FORMATS)
        if not path:
            return None
        width, height = self.pixel_size()
        image = picture(self.panel, width, height, self.transparent.isChecked())
        if image.isNull() or not save_image(image, path, self.dpi.value()):
            self.panel.status.setText(f"Could not export {Path(path).name}")
            return None
        self.panel.status.setText(f"Exported {Path(path).name} ({width:,} × {height:,} px)")
        return None


def _size_field() -> QDoubleSpinBox:
    field = QDoubleSpinBox()
    field.setRange(1, 100_000)
    field.setDecimals(0)
    return field
