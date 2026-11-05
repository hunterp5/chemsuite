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

"""The compound table view: structure drawing, selection sync, copy/paste and cell editing."""

from __future__ import annotations

import csv
import io
from collections import OrderedDict

import numpy as np
from PySide6.QtCore import QItemSelection, QItemSelectionModel, QModelIndex, QRect, QSize, Qt
from PySide6.QtGui import QColor, QGuiApplication, QImage, QPainter, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QHeaderView,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTableView,
    QWidget,
)

from chemsuite.core.depict import render_png
from chemsuite.core.selection import Selection
from chemsuite.core.table.columns import STRUCTURE, StructureColumn

from .table_model import KIND_ROLE, STATUS_ROLE, STRUCTURE_ROLE, TableModel

STRUCTURE_ROW_HEIGHT = 110
STRUCTURE_COLUMN_WIDTH = 200
MAX_AUTO_COLUMN_WIDTH = 320


class StructureImages:
    """Rendered structure pictures, keyed by molecule and size, within a byte budget."""

    def __init__(self, max_bytes: int = 192 * 1024**2) -> None:
        self._cache: OrderedDict[tuple, QImage] = OrderedDict()
        self._bytes = 0
        self.max_bytes = max_bytes

    def image(self, column: StructureColumn, slot: int, size: QSize, dpr: float, dark: bool):
        mol = column.mol(slot)  # parses a pending cell on first use
        if mol is None:
            return None
        blob = column.blobs(np.array([slot]))[0]
        key = (blob, size.width(), size.height(), round(dpr, 2), dark)
        image = self._cache.get(key)
        if image is not None:
            self._cache.move_to_end(key)
            return image
        png = render_png(
            mol, int(size.width() * dpr), int(size.height() * dpr), scale=dpr, dark=dark
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

    def paint(self, painter: QPainter, option: QStyleOptionViewItem, index: QModelIndex) -> None:
        if index.data(KIND_ROLE) != STRUCTURE:
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
            image = self._images.image(column, slot, rect.size(), dpr, dark)
            if image is not None:
                painter.drawImage(rect.topLeft(), image)
                return
        self._paint_invalid(painter, rect, column.display(slot))

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


class TableView(QTableView):
    """A table view bound to one :class:`TableModel` and the table's :class:`Selection`."""

    def __init__(
        self,
        model: TableModel,
        selection: Selection,
        images: StructureImages,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._origin = f"table:{id(self)}"
        self._selection = selection
        self._syncing = False
        self._anchor_oid: int | None = None
        self.setModel(model)
        self.setItemDelegate(CellDelegate(images, self))
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
            | QAbstractItemView.EditTrigger.AnyKeyPressed
        )
        self.setWordWrap(False)
        self.setAlternatingRowColors(True)
        self.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        vertical = self.verticalHeader()
        vertical.setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        horizontal = self.horizontalHeader()
        horizontal.setResizeContentsPrecision(200)  # never measure all 200k rows
        horizontal.setSortIndicator(-1, Qt.SortOrder.AscendingOrder)
        self.setSortingEnabled(True)

        model.modelAboutToBeReset.connect(self._remember_anchor)
        model.modelReset.connect(self._after_reset)
        model.layoutChanged.connect(self._restore_selection)
        self.selectionModel().selectionChanged.connect(self._on_view_selection)
        self._unsubscribe = selection.subscribe(self._on_store_selection)
        self._fit_layout()

    @property
    def table_model(self) -> TableModel:
        return self.model()  # type: ignore[return-value]

    def close_view(self) -> None:
        self._unsubscribe()

    # -- layout -------------------------------------------------------------------------------
    def _fit_layout(self) -> None:
        model = self.table_model
        has_structures = any(
            model.table.kind(name) == STRUCTURE for name in model.table.column_names
        )
        self.verticalHeader().setDefaultSectionSize(
            STRUCTURE_ROW_HEIGHT if has_structures else self.fontMetrics().height() + 8
        )
        self.resizeColumnsToContents()
        for col in range(model.columnCount()):
            name = model.column_name(col)
            if model.table.kind(name) == STRUCTURE:
                self.setColumnWidth(col, STRUCTURE_COLUMN_WIDTH)
            elif self.columnWidth(col) > MAX_AUTO_COLUMN_WIDTH:
                self.setColumnWidth(col, MAX_AUTO_COLUMN_WIDTH)

    # -- selection: view <-> store ------------------------------------------------------------
    def selected_rows(self) -> np.ndarray:
        """Display rows touched by the selection, ascending (vectorised over ranges)."""
        ranges = self.selectionModel().selection()
        parts = [np.arange(r.top(), r.bottom() + 1, dtype=np.int64) for r in ranges]
        return np.unique(np.concatenate(parts)) if parts else np.empty(0, dtype=np.int64)

    def selected_columns(self) -> list[int]:
        columns: set[int] = set()
        for r in self.selectionModel().selection():
            columns.update(range(r.left(), r.right() + 1))
        return sorted(columns)

    def selected_oids(self) -> np.ndarray:
        return self.table_model.visible_oids[self.selected_rows()]

    def _on_view_selection(self, *_args) -> None:
        if not self._syncing:
            self._selection.set(self.selected_oids().tolist(), self._origin)

    def _on_store_selection(self, oids: frozenset[int], origin: str) -> None:
        if origin != self._origin:
            self._apply_rows(self._rows_for(oids))

    def _rows_for(self, oids: frozenset[int]) -> np.ndarray:
        if not oids:
            return np.empty(0, dtype=np.int64)
        rows = self.table_model.rows_of_oids(np.fromiter(oids, dtype=np.int64, count=len(oids)))
        return np.sort(rows[rows >= 0])

    def _apply_rows(self, rows: np.ndarray) -> None:
        model = self.table_model
        selection = QItemSelection()
        if len(rows) and model.columnCount():
            breaks = np.flatnonzero(np.diff(rows) != 1)
            starts = np.concatenate([[rows[0]], rows[breaks + 1]])
            ends = np.concatenate([rows[breaks], [rows[-1]]])
            last = model.columnCount() - 1
            for start, end in zip(starts.tolist(), ends.tolist(), strict=True):
                selection.select(model.index(start, 0), model.index(end, last))
        self._syncing = True
        try:
            self.selectionModel().select(
                selection, QItemSelectionModel.SelectionFlag.ClearAndSelect
            )
        finally:
            self._syncing = False

    def _restore_selection(self) -> None:
        self._apply_rows(self._rows_for(self._selection.oids))

    def _remember_anchor(self) -> None:
        index = self.indexAt(self.viewport().rect().topLeft())
        self._anchor_oid = self.table_model.oid_at(index.row()) if index.isValid() else None

    def _after_reset(self) -> None:
        self._fit_layout()
        self._restore_selection()
        if self._anchor_oid is not None:
            rows = self.table_model.rows_of_oids([self._anchor_oid])
            if len(rows) and rows[0] >= 0:
                self.scrollTo(
                    self.table_model.index(int(rows[0]), 0),
                    QAbstractItemView.ScrollHint.PositionAtTop,
                )

    # -- clipboard ----------------------------------------------------------------------------
    def copy_selection(self) -> None:
        """Copy selected cells as tab-separated text (with a header row), Excel-compatible."""
        model = self.table_model
        rows = self.selected_rows()
        columns = self.selected_columns()
        if not len(rows):
            current = self.currentIndex()
            if not current.isValid():
                return
            rows, columns = np.array([current.row()]), [current.column()]
        names = [model.column_name(c) for c in columns]
        slots = model.table.slots_of(model.visible_oids[rows])
        cells = [model.table.column(name).displays(slots) for name in names]
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(names)
        writer.writerows(zip(*cells, strict=True))
        QGuiApplication.clipboard().setText(buffer.getvalue())

    def paste(self) -> str:
        """Paste tab-separated text at the current cell as one undoable step.

        A single copied value fills every selected cell. A first row equal to the target column
        names is treated as a header and skipped. Returns a status message.
        """
        model = self.table_model
        text = QGuiApplication.clipboard().text()
        if not text or not model.rowCount() or not model.columnCount():
            return "Nothing to paste"
        grid = [row for row in csv.reader(io.StringIO(text), delimiter="\t")]
        if not grid:
            return "Nothing to paste"
        if len(grid) == 1 and len(grid[0]) == 1 and len(self.selected_rows()) > 1:
            return self._fill_selection(grid[0][0])
        current = self.currentIndex()
        top = current.row() if current.isValid() else 0
        left = current.column() if current.isValid() else 0
        width = min(max(len(r) for r in grid), model.columnCount() - left)
        target_names = [model.column_name(left + c) for c in range(width)]
        if len(grid) > 1 and grid[0][:width] == target_names:
            grid = grid[1:]  # a header row copied along with the cells
        height = min(len(grid), model.rowCount() - top)
        oids = model.visible_oids[top : top + height]
        updates = {}
        for c in range(width):
            values = [grid[r][c] if c < len(grid[r]) else "" for r in range(height)]
            updates[model.column_name(left + c)] = (oids, values)
        model.table.set_block(updates, label="Paste")
        skipped = len(grid) - height
        return f"Pasted {height} row(s)" + (
            f"; {skipped} row(s) past the end skipped" if skipped else ""
        )

    def selected_cells_by_column(self) -> dict[str, np.ndarray]:
        """``{column name: oids}`` for every selected cell."""
        model = self.table_model
        cells: dict[str, np.ndarray] = {}
        for r in self.selectionModel().selection():
            oids = model.visible_oids[r.top() : r.bottom() + 1]
            for col in range(r.left(), r.right() + 1):
                name = model.column_name(col)
                cells[name] = np.union1d(cells.get(name, np.empty(0, dtype=np.int64)), oids)
        return cells

    def _fill_selection(self, value: str) -> str:
        cells = self.selected_cells_by_column()
        self.table_model.table.set_block(
            {name: (oids, [value] * len(oids)) for name, oids in cells.items()}, label="Paste"
        )
        return f"Filled {sum(len(oids) for oids in cells.values())} cell(s)"

    def clear_selected_cells(self) -> None:
        cells = self.selected_cells_by_column()
        if cells:
            self.table_model.table.set_block(
                {name: (oids, [""] * len(oids)) for name, oids in cells.items()},
                label="Clear cells",
            )

    def keyPressEvent(self, event) -> None:
        if self.state() != QAbstractItemView.State.EditingState and event.key() in (
            Qt.Key.Key_Delete,
            Qt.Key.Key_Backspace,
        ):
            self.clear_selected_cells()
            return
        super().keyPressEvent(event)
