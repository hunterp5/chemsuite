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
"""A read-only, sortable table of analysis results, on a model: filling it with a hundred
thousand rows costs nothing on the GUI thread (cells are formatted only when drawn).

A cell can be a :class:`Drawing` (a structure): drawn only when it scrolls into view, and
sorted and copied by its text."""

from __future__ import annotations

import math
from collections import OrderedDict
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    QPersistentModelIndex,
    QRect,
    QSize,
    Qt,
    Signal,
)
from PySide6.QtGui import QGuiApplication, QImage, QPainter
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTableView,
)

_NUMBER = (int, float, np.integer, np.floating)
_RIGHT = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
#: rows measured to size the columns (Qt measures every row by default)
_MEASURED_ROWS = 100
#: the role that hands a cell's :class:`Drawing` to the delegate
DRAWING_ROLE = Qt.ItemDataRole.UserRole + 1
#: drawings kept, at the sizes they were last drawn
_CACHED_IMAGES = 400


@dataclass(frozen=True, eq=False)
class Drawing:
    """A cell drawn as a picture: *draw* ``(width, height, scale, dark)`` gives PNG bytes in
    device pixels (None when it cannot be drawn: the text is shown instead). *text* is what
    the cell sorts and copies by, and its tooltip; *caption* is written under the picture."""

    text: str
    draw: Callable[[int, int, float, bool], bytes | None]
    caption: str = ""

    def __str__(self) -> str:
        return self.text


def _shown(value: object) -> object:
    """What a cell shows: numbers as numbers (floats to 3 places; not-a-number blank)."""
    if isinstance(value, _NUMBER):
        number = float(value)
        if not math.isfinite(number):
            return None
        return round(number, 3) if isinstance(value, float | np.floating) else int(value)
    return str(value)


def _text(value: object) -> str:
    shown = _shown(value)
    return "" if shown is None else str(shown)


class ResultModel(QAbstractTableModel):
    """Rows of raw cell values in a chosen order; the rows themselves never move."""

    def __init__(self, headers: Sequence[str]) -> None:
        super().__init__()
        self.headers = list(headers)
        self.rows: list[Sequence[object]] = []
        self.order: list[int] = []

    def set_rows(self, rows: Sequence[Sequence[object]]) -> None:
        self.beginResetModel()
        self.rows = list(rows)
        self.order = list(range(len(self.rows)))
        self.endResetModel()

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:  # noqa: B008
        return 0 if parent.isValid() else len(self.order)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:  # noqa: B008
        return 0 if parent.isValid() else len(self.headers)

    def value(self, row: int, column: int) -> object:
        cells = self.rows[self.order[row]]
        return cells[column] if column < len(cells) else ""

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        value = self.value(index.row(), index.column())
        if role == Qt.ItemDataRole.DisplayRole:
            return _shown(value)
        if role == Qt.ItemDataRole.TextAlignmentRole and isinstance(value, _NUMBER):
            return _RIGHT
        if isinstance(value, Drawing):
            if role == DRAWING_ROLE:
                return value
            if role == Qt.ItemDataRole.ToolTipRole:
                return value.text
        return None

    def headerData(self, section: int, orientation, role: int = Qt.ItemDataRole.DisplayRole):
        if (
            role == Qt.ItemDataRole.DisplayRole
            and orientation == Qt.Orientation.Horizontal
            and 0 <= section < len(self.headers)
        ):
            return self.headers[section]
        return None

    def set_headers(self, headers: Sequence[str]) -> None:
        self.headers = list(headers)
        self.headerDataChanged.emit(Qt.Orientation.Horizontal, 0, len(self.headers) - 1)

    def sort(self, column: int, order: Qt.SortOrder = Qt.SortOrder.AscendingOrder) -> None:
        """Sort by *column* (numbers as numbers, blanks last); -1 is the given order."""
        self.layoutAboutToBeChanged.emit()
        if column < 0 or column >= len(self.headers):
            self.order = list(range(len(self.rows)))
        else:
            keyed, blank = [], []
            for i, cells in enumerate(self.rows):
                shown = _shown(cells[column]) if column < len(cells) else None
                (blank if shown is None or shown == "" else keyed).append((shown, i))
            numeric = all(isinstance(v, int | float) for v, _i in keyed)
            keyed.sort(
                key=(lambda vi: vi[0]) if numeric else (lambda vi: str(vi[0]).casefold()),
                reverse=order == Qt.SortOrder.DescendingOrder,
            )
            self.order = [i for _v, i in keyed] + [i for _v, i in blank]
        self.layoutChanged.emit()


class DrawingDelegate(QStyledItemDelegate):
    """Paints :class:`Drawing` cells (other cells as usual), keeping recent pictures."""

    def __init__(self, width: int, parent=None) -> None:
        super().__init__(parent)
        self.width = width
        self.images: OrderedDict[tuple, QImage] = OrderedDict()
        self.drawn = 0  # pictures made (not taken from the cache)

    def sizeHint(self, option: QStyleOptionViewItem, index: QModelIndex | QPersistentModelIndex):
        hint = super().sizeHint(option, index)
        if isinstance(index.data(DRAWING_ROLE), Drawing):
            return QSize(max(self.width, hint.width()), hint.height())
        return hint

    def paint(
        self,
        painter: QPainter,
        option: QStyleOptionViewItem,
        index: QModelIndex | QPersistentModelIndex,
    ) -> None:
        drawing = index.data(DRAWING_ROLE)
        if not isinstance(drawing, Drawing):
            super().paint(painter, option, index)
            return
        item = QStyleOptionViewItem(option)
        self.initStyleOption(item, index)
        item.text = ""
        widget = item.widget
        style = widget.style() if widget is not None else QApplication.style()
        style.drawControl(QStyle.ControlElement.CE_ItemViewItem, item, painter, widget)
        box = option.rect.adjusted(2, 2, -2, -2)
        if drawing.caption:
            line = option.fontMetrics.height()
            caption = QRect(box.left(), box.bottom() - line + 1, box.width(), line)
            box.setBottom(caption.top() - 1)
            elided = option.fontMetrics.elidedText(
                drawing.caption, Qt.TextElideMode.ElideRight, caption.width()
            )
            painter.drawText(caption, Qt.AlignmentFlag.AlignCenter, elided)
        scale = widget.devicePixelRatioF() if widget is not None else 1.0
        dark = option.palette.base().color().lightness() < 128
        image = self.image(drawing, box.width(), box.height(), scale, dark)
        if image is None:
            painter.drawText(box, Qt.AlignmentFlag.AlignCenter, drawing.text)
            return
        painter.drawImage(box.topLeft(), image)

    def image(
        self, drawing: Drawing, width: int, height: int, scale: float, dark: bool
    ) -> QImage | None:
        if width < 8 or height < 8:
            return None
        key = (drawing, width, height, scale, dark)
        if key in self.images:
            self.images.move_to_end(key)
            return self.images[key]
        png = drawing.draw(round(width * scale), round(height * scale), scale, dark)
        self.drawn += 1
        if not png:
            return None
        image = QImage.fromData(png, "PNG")
        if image.isNull():
            return None
        image.setDevicePixelRatio(scale)
        self.images[key] = image
        while len(self.images) > _CACHED_IMAGES:
            self.images.popitem(last=False)
        return image


class ResultTable(QTableView):
    """A read-only, sortable table of results. With :class:`Drawing` cells, give a
    *row_height* and *drawing_width* to suit the pictures."""

    #: the selected rows changed
    itemSelectionChanged = Signal()
    #: a cell was double-clicked: its row (as shown) and column
    cellDoubleClicked = Signal(int, int)

    def __init__(
        self, headers: Sequence[str], *, row_height: int = 22, drawing_width: int = 160
    ) -> None:
        super().__init__()
        self.result_model = ResultModel(headers)
        self.setModel(self.result_model)
        self.drawings = DrawingDelegate(drawing_width, self)
        self.setItemDelegate(self.drawings)
        self.selectionModel().selectionChanged.connect(lambda *_: self.itemSelectionChanged.emit())
        self.doubleClicked.connect(lambda i: self.cellDoubleClicked.emit(i.row(), i.column()))
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.verticalHeader().hide()
        self.verticalHeader().setDefaultSectionSize(row_height)
        if row_height > 40:  # tall rows: scroll smoothly, not a whole picture at a time
            self.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        header = self.horizontalHeader()
        header.setStretchLastSection(False)
        header.setResizeContentsPrecision(_MEASURED_ROWS)
        self._data: list | None = None

    def fill(self, rows: Sequence[Sequence[object]], data: Sequence[object] | None = None) -> None:
        """One list of cells per row; numbers align right and sort as numbers.

        *data* attaches a value per row (read back with :meth:`row_data`), so sorting by a
        column does not lose track of which result a row shows.
        """
        self.setSortingEnabled(False)
        self.result_model.set_rows(rows)
        self._data = list(data) if data is not None else None
        self.resizeColumnsToContents()
        for c in range(self.result_model.columnCount()):  # room for the sort arrow
            self.setColumnWidth(c, min(self.columnWidth(c) + 18, 360))
        # Keep the given order (e.g. best first) until the user clicks a header.
        self.horizontalHeader().setSortIndicator(-1, Qt.SortOrder.AscendingOrder)
        self.setSortingEnabled(True)

    def setHorizontalHeaderLabels(self, headers: Sequence[str]) -> None:
        self.result_model.set_headers(headers)

    def rowCount(self) -> int:
        return self.result_model.rowCount()

    def header_text(self, column: int) -> str:
        return self.result_model.headers[column]

    def row_data(self, row: int):
        if self._data is None or not 0 <= row < len(self.result_model.order):
            return None
        return self._data[self.result_model.order[row]]

    def selected_data(self) -> list:
        rows = sorted(i.row() for i in self.selectionModel().selectedRows())
        return [self.row_data(r) for r in rows]

    def copy_text(self) -> str:
        model = self.result_model
        lines = ["\t".join(model.headers)]
        for r in range(model.rowCount()):
            lines.append("\t".join(_text(model.value(r, c)) for c in range(model.columnCount())))
        text = "\n".join(lines)
        QGuiApplication.clipboard().setText(text)
        return text
