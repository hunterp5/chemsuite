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
from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import (
    QElapsedTimer,
    QEvent,
    QItemSelection,
    QItemSelectionModel,
    QModelIndex,
    QPoint,
    QPointF,
    QRect,
    QSize,
    Qt,
    QTimer,
    Signal,
)
from PySide6.QtGui import (
    QColor,
    QFont,
    QGuiApplication,
    QImage,
    QMouseEvent,
    QPainter,
    QPalette,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QFrame,
    QHeaderView,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTableView,
    QWidget,
)
from rdkit import Chem

from chemsuite.analysis.structure_search import match_atoms
from chemsuite.core.depict import depiction_mol, render_png, render_reaction_png
from chemsuite.core.mols import looks_like_molblock
from chemsuite.core.reactions import looks_like_rxn_block, parse_reaction
from chemsuite.core.selection import Selection
from chemsuite.core.structure_input import read_structure
from chemsuite.core.table.columns import REACTION, STRUCTURE, ParsedColumn, ReactionColumn

from .structure_clipboard import structure_text_from_mime
from .table_model import KIND_ROLE, STATUS_ROLE, STRUCTURE_ROLE, TableModel

STRUCTURE_ROW_HEIGHT = 110
STRUCTURE_COLUMN_WIDTH = 200
REACTION_COLUMN_WIDTH = 440
#: Column kinds whose cells are drawn as pictures.
DRAWN = (STRUCTURE, REACTION)
MAX_AUTO_COLUMN_WIDTH = 320
#: keys that leave a column selection made on the header in place
_KEEP_COLUMN_MODE = frozenset(
    {
        Qt.Key.Key_Delete,
        Qt.Key.Key_Backspace,
        Qt.Key.Key_Shift,
        Qt.Key.Key_Control,
        Qt.Key.Key_Meta,
        Qt.Key.Key_Alt,
    }
)


@dataclass(frozen=True)
class Highlight:
    """Atoms to highlight in a structure column: those matching *query* (a substructure)."""

    key: str  # identifies the query in the image cache
    query: Chem.Mol
    stereo: bool = False


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
        highlight: Highlight | None = None,
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
            if highlight is not None:
                mol = depiction_mol(mol)
                atoms, bonds = match_atoms(mol, highlight.query, highlight.stereo)
            png = render_png(
                mol,
                width,
                height,
                scale=dpr,
                dark=dark,
                highlight_atoms=atoms,
                highlight_bonds=bonds,
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
        #: column name -> atoms to highlight in that column's drawings
        self.highlights: dict[str, Highlight] = {}

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
            highlight = None
            if self.highlights:
                highlight = self.highlights.get(index.model().column_name(index.column()))
            image = self._images.image(column, slot, rect.size(), dpr, dark, highlight)
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


#: Within this many pixels of the table's edge, a dragged column scrolls the table.
EDGE_ZONE = 56
#: Scrolling speed (pixels per second) entering the edge zone, and the fastest it gets.
EDGE_MIN_SPEED, EDGE_MAX_SPEED = 150.0, 2400.0


def edge_speed(x: int, width: int) -> float:
    """Sideways scroll speed for a column dragged at *x* in a header *width* wide.

    Zero away from the edges. Inside the edge zone it grows the closer the pointer gets and
    keeps growing past the edge (up to two zones beyond), so a short way is precise and a long
    table is crossed quickly. Negative scrolls left.
    """
    zone = min(EDGE_ZONE, width // 4)
    if zone <= 0:
        return 0.0
    if x < zone:
        depth, sign = (zone - x) / zone, -1.0
    elif x > width - zone:
        depth, sign = (x - (width - zone)) / zone, 1.0
    else:
        return 0.0
    t = min(depth, 3.0) / 3.0  # 0 at the zone's inner edge, 1/3 at the edge, 1 far beyond
    return sign * (EDGE_MIN_SPEED + (EDGE_MAX_SPEED - EDGE_MIN_SPEED) * t * t)


class TableView(QTableView):
    """A table view bound to one :class:`TableModel` and the table's :class:`Selection`.

    Column widths are remembered by column name, so a model reset (filtering, adding rows)
    keeps widths the user set; only new columns, or columns whose kind changed, are fitted.
    Dragging a header moves the column in the table (one undoable step) rather than only in
    the view. Clicking a header selects the column (Shift and Ctrl/⌘ extend it); sorting is in
    the header's menu. Whole columns selected from the header are not a compound selection, so
    they leave plots and other views alone, and Delete removes them.
    """

    #: (column index, global position) for a right-click on the horizontal header
    headerMenuRequested = Signal(int, QPoint)
    #: (cell, global position) for a right-click on a cell
    cellMenuRequested = Signal(QModelIndex, QPoint)
    #: (cell, starting text or "" to keep the cell's) when a structure cell is to be edited
    structureEditRequested = Signal(QModelIndex, str)
    #: a structure cell was double-clicked: open the sketcher on it
    structureDrawRequested = Signal(QModelIndex)
    #: column names to delete (Delete pressed with whole columns selected)
    deleteColumnsRequested = Signal(list)

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
        self._anchor_x = 0  # horizontal scroll kept across a model reset
        self._anchor_columns: list[str] = []  # whole columns selected before a reset
        self._reveal: str | None = None  # a moved column to keep in view after the reset
        # Dragging a column header near the table's left or right edge scrolls the table.
        self._drag_press: QPoint | None = None
        self._drag_pos: QPoint | None = None
        self._drag_speed = 0.0
        self._drag_clock = QElapsedTimer()
        self._drag_timer = QTimer(self)
        self._drag_timer.setInterval(16)
        self._drag_timer.timeout.connect(self._drag_scroll)
        self._column_mode = False  # the selection was made on the column header
        self._widths: dict[str, int] = {}
        self._kinds: dict[str, str] = {}
        self._structure_row_height = STRUCTURE_ROW_HEIGHT
        self._frozen = None  # the FrozenView over the left edge, made when first needed
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
        self.setFrameShape(QFrame.Shape.NoFrame)  # flush with the window and its panes
        vertical = self.verticalHeader()
        vertical.setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        vertical.setDefaultAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        # Row numbers: flat, muted, no per-row separators; one hairline against the cells.
        vertical.setStyleSheet(
            "QHeaderView::section { background: palette(window); color: palette(placeholder-text);"
            " border: none; border-right: 1px solid palette(mid); padding: 0 6px 0 8px; }"
        )
        horizontal = self.horizontalHeader()
        horizontal.setResizeContentsPrecision(200)  # never measure all 200k rows
        horizontal.setSortIndicator(-1, Qt.SortOrder.AscendingOrder)
        horizontal.setSortIndicatorShown(True)
        horizontal.setSectionsClickable(True)  # a click selects the column; sorting is a menu
        # Note a header press before Qt selects the column, so that selection is already known
        # to be a column selection when it reaches the compound store.
        horizontal.viewport().installEventFilter(self)
        self.verticalHeader().viewport().installEventFilter(self)
        horizontal.setSectionsMovable(True)
        horizontal.sectionMoved.connect(self._on_section_moved)
        horizontal.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        horizontal.customContextMenuRequested.connect(self._on_header_menu)
        # A header click flips the indicator but sorts nothing: put it back.
        horizontal.sortIndicatorChanged.connect(
            lambda *_: QTimer.singleShot(0, self._sync_sort_indicator)
        )

        model.modelAboutToBeReset.connect(self._remember_anchor)
        model.modelReset.connect(self._after_reset)
        model.layoutChanged.connect(self._restore_selection)
        model.layoutChanged.connect(self._sync_sort_indicator)
        model.headerDataChanged.connect(self._on_header_data_changed)
        self.selectionModel().selectionChanged.connect(self._on_view_selection)
        self._unsubscribe = selection.subscribe(self._on_store_selection)
        self._fit_layout()

    def apply_appearance(self, appearance) -> None:
        """Table font, text alignment and structure row height from Settings."""
        font = QFont(QApplication.font())
        if appearance.table_font_pt:
            font.setPointSize(appearance.table_font_pt)
        changed = (
            font != self.font() or appearance.structure_row_height != self._structure_row_height
        )
        self.setFont(font)
        self._structure_row_height = appearance.structure_row_height
        self.table_model.set_alignment(appearance.table_align_h, appearance.table_align_v)
        if changed:
            self._capture_widths()
            self._fit_layout()

    def set_highlight(self, column: str, highlight: Highlight | None) -> None:
        """Highlight atoms matching a query in *column*'s drawings (None: stop)."""
        highlights = self.itemDelegate().highlights
        if highlight is None:
            if highlights.pop(column, None) is None:
                return
        else:
            highlights[column] = highlight
        self.viewport().update()

    def highlights(self) -> dict[str, Highlight]:
        return dict(self.itemDelegate().highlights)

    @property
    def table_model(self) -> TableModel:
        return self.model()  # type: ignore[return-value]

    def close_view(self) -> None:
        self._unsubscribe()

    # -- layout -------------------------------------------------------------------------------
    def _capture_widths(self) -> None:
        """Remember current widths by name (call while the header still matches the model)."""
        header = self.horizontalHeader()
        model = self.table_model
        for col in range(min(header.count(), model.columnCount())):
            self._widths[model.column_name(col)] = header.sectionSize(col)

    def column_widths(self) -> dict[str, int]:
        """Current column widths by name (saved with the document)."""
        self._capture_widths()
        names = set(self.table_model.table.column_names)
        return {name: width for name, width in self._widths.items() if name in names}

    def set_column_widths(self, widths: dict[str, int]) -> None:
        self._capture_widths()
        self._widths.update(widths)
        self._fit_layout()

    def _fit_layout(self) -> None:
        model = self.table_model
        table = model.table
        has_structures = any(table.kind(name) in DRAWN for name in table.column_names)
        self.verticalHeader().setDefaultSectionSize(
            self._structure_row_height if has_structures else self.fontMetrics().height() + 8
        )
        for col in range(model.columnCount()):
            name = model.column_name(col)
            kind = table.kind(name)
            known = self._widths.get(name)
            if known is not None and self._kinds.get(name) == kind:
                self.setColumnWidth(col, known)
            elif kind == STRUCTURE:
                self.setColumnWidth(col, STRUCTURE_COLUMN_WIDTH)
            elif kind == REACTION:
                self.setColumnWidth(col, REACTION_COLUMN_WIDTH)
            else:
                self.resizeColumnToContents(col)
                self.setColumnWidth(col, min(self.columnWidth(col), MAX_AUTO_COLUMN_WIDTH))
        self._kinds = {name: table.kind(name) for name in table.column_names}
        self._sync_sort_indicator()
        if self._frozen is not None:
            self._frozen.refresh()

    def _on_header_data_changed(self, orientation: Qt.Orientation, *_args) -> None:
        if orientation == Qt.Orientation.Horizontal:  # e.g. a column changed kind
            self._capture_widths()
            self._fit_layout()

    def _sync_sort_indicator(self) -> None:
        """Show the table's primary sort key without triggering another sort."""
        model = self.table_model
        header = self.horizontalHeader()
        keys = model.table.sort_keys
        blocked = header.blockSignals(True)
        try:
            if keys and keys[0][0] in model.table:
                order = Qt.SortOrder.AscendingOrder if keys[0][1] else Qt.SortOrder.DescendingOrder
                header.setSortIndicator(model.column_index(keys[0][0]), order)
            else:
                header.setSortIndicator(-1, Qt.SortOrder.AscendingOrder)
        finally:
            header.blockSignals(blocked)

    def _on_section_moved(self, logical: int, old_visual: int, new_visual: int) -> None:
        """Turn a header drag into a table column move (undoable, saved with the document)."""
        header = self.horizontalHeader()
        blocked = header.blockSignals(True)
        try:
            header.moveSection(new_visual, old_visual)  # the view itself stays unmoved
        finally:
            header.blockSignals(blocked)
        name = self.table_model.column_name(logical)
        self._reveal = name
        QTimer.singleShot(0, lambda: self.table_model.table.move_column(name, new_visual))

    def _on_header_menu(self, pos: QPoint) -> None:
        column = self.horizontalHeader().logicalIndexAt(pos)
        if column >= 0:
            self.headerMenuRequested.emit(column, self.horizontalHeader().mapToGlobal(pos))

    def eventFilter(self, watched, event) -> bool:
        if event.type() == QEvent.Type.MouseButtonPress:
            self._column_mode = watched is self.horizontalHeader().viewport()
        if watched is self.horizontalHeader().viewport():
            self._follow_header_drag(event)
        return super().eventFilter(watched, event)

    # -- dragging columns past the edge ---------------------------------------------------------
    def _follow_header_drag(self, event) -> None:
        """While a column is dragged, scroll when the pointer nears or passes an edge."""
        kind = event.type()
        if kind == QEvent.Type.MouseButtonPress and event.button() == Qt.MouseButton.LeftButton:
            self._drag_press = event.position().toPoint()
        elif kind == QEvent.Type.MouseButtonRelease:
            self._stop_drag_scroll()
            self._drag_press = None
        elif (
            kind == QEvent.Type.MouseMove
            and self._drag_press is not None
            and event.buttons() & Qt.MouseButton.LeftButton
        ):
            pos = event.position().toPoint()
            if (pos - self._drag_press).manhattanLength() < QApplication.startDragDistance():
                return
            self._drag_pos = pos
            self._drag_speed = edge_speed(pos.x(), self.horizontalHeader().viewport().width())
            if self._drag_speed and not self._drag_timer.isActive():
                self._drag_clock.start()
                self._drag_timer.start()
            elif not self._drag_speed:
                self._stop_drag_scroll()

    def _stop_drag_scroll(self) -> None:
        self._drag_timer.stop()
        self._drag_speed = 0.0

    def _drag_scroll(self) -> None:
        bar = self.horizontalScrollBar()
        before = bar.value()
        # pixels per second, scaled by the time since the last step: smooth at any frame rate
        elapsed = max(1, self._drag_clock.restart())
        step = round(self._drag_speed * elapsed / 1000) or (1 if self._drag_speed > 0 else -1)
        bar.setValue(before + step)
        if bar.value() == before or self._drag_pos is None:
            return
        # Tell the header the pointer is still there so its drop marker follows the scroll.
        viewport = self.horizontalHeader().viewport()
        pos = QPointF(self._drag_pos)
        move = QMouseEvent(
            QEvent.Type.MouseMove,
            pos,
            viewport.mapToGlobal(pos),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            QGuiApplication.keyboardModifiers(),
        )
        QApplication.sendEvent(viewport, move)

    def whole_columns(self) -> list[str]:
        """The selected column names when whole columns were picked on the header, else []."""
        model = self.table_model
        rows = model.rowCount()
        ranges = self.selectionModel().selection()
        if not self._column_mode or not rows or not len(ranges):
            return []
        if any(r.top() != 0 or r.bottom() != rows - 1 for r in ranges):
            return []
        return [model.column_name(c) for c in self.selected_columns()]

    def selected_column_names(self) -> list[str]:
        """Names of the columns that are selected top to bottom."""
        model = self.table_model
        selection = self.selectionModel()
        return [
            model.column_name(col)
            for col in range(model.columnCount())
            if selection.isColumnSelected(col, QModelIndex())
        ]

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
        """Oids of the selected rows; none while whole columns are selected on the header."""
        if self.whole_columns():
            return np.empty(0, dtype=np.int64)
        return self.table_model.visible_oids[self.selected_rows()]

    def _on_view_selection(self, *_args) -> None:
        if self._syncing:
            return
        if self.whole_columns():
            self._selection.set([], self._origin)  # columns are not compounds
        else:
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
        self._column_mode = False
        try:
            self.selectionModel().select(
                selection, QItemSelectionModel.SelectionFlag.ClearAndSelect
            )
        finally:
            self._syncing = False

    def _restore_selection(self) -> None:
        if not self._selection.oids and self.whole_columns():
            return  # a column selection survives re-sorting
        self._apply_rows(self._rows_for(self._selection.oids))

    def _remember_anchor(self) -> None:
        self._capture_widths()
        index = self.indexAt(self.viewport().rect().topLeft())
        self._anchor_oid = self.table_model.oid_at(index.row()) if index.isValid() else None
        self._anchor_x = self.horizontalScrollBar().value()
        self._anchor_columns = self.whole_columns()

    def _after_reset(self) -> None:
        self._fit_layout()
        columns = [n for n in self._anchor_columns if n in self.table_model.table]
        self._anchor_columns = []
        if columns:
            self._select_columns(columns)
        else:
            self._restore_selection()
        if self._anchor_oid is not None:
            rows = self.table_model.rows_of_oids([self._anchor_oid])
            if len(rows) and rows[0] >= 0:
                self.scrollTo(
                    self.table_model.index(int(rows[0]), 0),
                    QAbstractItemView.ScrollHint.PositionAtTop,
                )
        # stay where the user was sideways (scrollTo above went back to the first column)
        self.horizontalScrollBar().setValue(self._anchor_x)
        reveal, self._reveal = self._reveal, None
        if reveal is not None and reveal in self.table_model.table:
            self.reveal_column(reveal)

    def _select_columns(self, names: list[str]) -> None:
        """Select whole columns by name, as a header click does."""
        model = self.table_model
        rows = model.rowCount()
        if not rows:
            return
        selection = QItemSelection()
        for name in names:
            col = model.column_index(name)
            selection.select(model.index(0, col), model.index(rows - 1, col))
        self._column_mode = True
        self.selectionModel().select(selection, QItemSelectionModel.SelectionFlag.ClearAndSelect)

    def reveal_column(self, name: str) -> None:
        """Scroll sideways just enough to show column *name* whole (if it fits)."""
        header = self.horizontalHeader()
        col = self.table_model.column_index(name)
        left = header.sectionPosition(col)
        right = left + header.sectionSize(col)
        bar = self.horizontalScrollBar()
        width = self.viewport().width()
        if left < bar.value():
            bar.setValue(left)
        elif right > bar.value() + width:
            bar.setValue(min(left, right - width))

    def mousePressEvent(self, event) -> None:
        self._column_mode = False  # a cell click starts an ordinary selection
        super().mousePressEvent(event)

    def start_column_selection(self) -> None:
        """A header press (here or on the frozen columns' header): whole columns follow."""
        self._column_mode = True

    def start_cell_selection(self) -> None:
        self._column_mode = False

    # -- frozen columns -----------------------------------------------------------------------
    @property
    def frozen_count(self) -> int:
        """How many columns, from the left, stay put while the rest scroll sideways."""
        return self._frozen.count if self._frozen is not None else 0

    def set_frozen(self, count: int) -> None:
        """Freeze the first *count* columns (0 unfreezes)."""
        count = max(0, min(count, self.table_model.columnCount()))
        if count == self.frozen_count:
            return
        if self._frozen is None:
            from .frozen_columns import FrozenView

            self._frozen = FrozenView(self)
        self._frozen.set_count(count)
        self.horizontalScrollBar().setValue(0)

    def updateGeometries(self) -> None:
        super().updateGeometries()
        if self._frozen is not None:
            self._frozen.place()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if self._frozen is not None:
            self._frozen.place()

    def scrollTo(self, index, hint=QAbstractItemView.ScrollHint.EnsureVisible) -> None:
        """Scroll as usual, then make sure a cell right of the frozen part is not under it.

        A frozen cell is always in view sideways, so going to one only scrolls up or down.
        """
        frozen = self.frozen_count
        sideways = self.horizontalScrollBar().value()
        super().scrollTo(index, hint)
        if not frozen or not index.isValid():
            return
        if index.column() < frozen:
            self.horizontalScrollBar().setValue(sideways)
            return
        covered = self._frozen.frozen_width()
        left = self.columnViewportPosition(index.column())
        if left < covered:
            bar = self.horizontalScrollBar()
            bar.setValue(bar.value() - (covered - left))

    def contextMenuEvent(self, event) -> None:
        index = self.indexAt(event.pos())
        if index.isValid():
            self.cellMenuRequested.emit(index, event.globalPos())

    # -- editing ------------------------------------------------------------------------------
    def edit(self, index, trigger=None, event=None):  # type: ignore[override]
        """Structure and reaction cells open the sketcher on double-click and F2. Typing in a
        structure cell opens the sketcher too; typing in a reaction cell edits its text in
        place."""
        if trigger is None:
            return super().edit(index)
        kind = index.data(KIND_ROLE) if index.isValid() else None
        if kind not in DRAWN or self.state() == QAbstractItemView.State.EditingState:
            return super().edit(index, trigger, event)
        Trigger = QAbstractItemView.EditTrigger
        if trigger == Trigger.DoubleClicked:
            self.structureDrawRequested.emit(index)
            return False
        if trigger == Trigger.EditKeyPressed or (
            kind == STRUCTURE and trigger == Trigger.AnyKeyPressed
        ):
            typed = event.text() if trigger == Trigger.AnyKeyPressed else ""
            self.structureEditRequested.emit(index, typed)
            return False
        return super().edit(index, trigger, event)

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
        if not model.rowCount() or not model.columnCount():
            return "Nothing to paste"
        mime = QGuiApplication.clipboard().mimeData()
        rxn = _rxn_text_from_mime(mime)
        if rxn:
            return self._paste_rxn(rxn)
        structure = structure_text_from_mime(mime)
        if looks_like_molblock(structure):
            return self._paste_molfile(structure)
        text = QGuiApplication.clipboard().text()
        if not text:
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

    def _target_cells(self, kind: str) -> dict[str, np.ndarray]:
        """Selected cells in columns of *kind*, else the current cell if it is one."""
        table = self.table_model.table
        cells = {
            name: oids
            for name, oids in self.selected_cells_by_column().items()
            if table.kind(name) == kind
        }
        current = self.currentIndex()
        if not cells and current.isValid() and current.data(KIND_ROLE) == kind:
            name = self.table_model.column_name(current.column())
            cells = {name: np.array([self.table_model.oid_at(current.row())])}
        return cells

    def _paste_rxn(self, block: str) -> str:
        """An RXN file goes whole into the selected reaction cells (or the current one)."""
        cells = self._target_cells(REACTION)
        if not cells:
            return "Select a reaction cell to paste the RXN file into"
        blob, shown = parse_reaction(block)
        if not blob:
            return "The pasted RXN file is not a valid reaction"
        self.table_model.table.set_block(
            {
                name: (oids, ReactionColumn.from_parsed([blob] * len(oids), [shown] * len(oids)))
                for name, oids in cells.items()
            },
            label="Paste reaction",
        )
        return f"Pasted reaction into {sum(len(o) for o in cells.values())} cell(s)"

    def _paste_molfile(self, block: str) -> str:
        """A molfile goes whole into the selected structure cells (or the current one)."""
        table = self.table_model.table
        cells = self._target_cells(STRUCTURE)
        if not cells:
            return "Select a structure cell to paste the molfile into"
        parsed = read_structure(block)
        if parsed.mol is None:
            return f"The pasted molfile is not valid: {parsed.error}"
        one = parsed.as_column()
        table.set_block(
            {
                name: (oids, one.take(np.zeros(len(oids), dtype=np.int64)))
                for name, oids in cells.items()
            },
            label="Paste structure",
        )
        return f"Pasted structure into {sum(len(o) for o in cells.values())} cell(s)"

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
        if event.key() not in _KEEP_COLUMN_MODE:
            self._column_mode = False  # Ctrl+A, arrows …: an ordinary cell selection
        if self.state() != QAbstractItemView.State.EditingState and event.key() in (
            Qt.Key.Key_Delete,
            Qt.Key.Key_Backspace,
        ):
            columns = self.whole_columns()
            if columns:
                self.deleteColumnsRequested.emit(columns)
            else:
                self.clear_selected_cells()
            return
        super().keyPressEvent(event)


def _rxn_text_from_mime(mime) -> str:
    """An RXN file on the clipboard (as ``chemical/x-mdl-rxnfile`` or as text), else ""."""
    if mime is None:
        return ""
    if mime.hasFormat("chemical/x-mdl-rxnfile"):
        text = bytes(mime.data("chemical/x-mdl-rxnfile").data()).decode("utf-8", "replace")
        if looks_like_rxn_block(text):
            return text
    text = mime.text() if mime.hasText() else ""
    return text if looks_like_rxn_block(text) else ""
