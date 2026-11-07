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

"""The selection browser: step through rows one at a time with a large drawing.

It browses either the selected rows or every shown row, in the table's display order, and
keeps the table's current row in step (without changing the selection) so the two views
always point at the same compound. Its tick box adds or removes the row from the selection.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QEvent, QItemSelectionModel, QModelIndex, Qt, QTimer
from PySide6.QtGui import QKeyEvent, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.depict import render_png
from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed
from chemsuite.core.table.columns import STRUCTURE

from .table_view import TableView

ORIGIN = "browser"
SCOPES = {"selected": "Selected rows", "shown": "All shown rows"}
PREVIEW_SCALE = 1.6  # bonds drawn larger than in table cells (shrunk to fit if needed)


class _Preview(QLabel):
    """A drawing that re-renders at its own size after resizes settle."""

    def __init__(self, owner: SelectionBrowser) -> None:
        super().__init__()
        self._owner = owner
        self.setObjectName("BrowserPreview")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumSize(160, 140)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(60)
        self._timer.timeout.connect(owner._draw)

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._timer.start()


class SelectionBrowser(QWidget):
    def __init__(self, view: TableView, selection: Selection, parent: QWidget | None = None):
        super().__init__(parent)
        self.view = view
        self.model = view.table_model
        self.table = self.model.table
        self.selection = selection
        self._oids = np.empty(0, dtype=np.int64)  # rows in scope, display order
        self._index = -1
        self._oid: int | None = None
        self._following = False  # moving the table cursor ourselves
        self._stale = False

        self.scope = QComboBox()
        for key, title in SCOPES.items():
            self.scope.addItem(title, key)
        self.scope.setToolTip("Which rows to step through")
        self.structure = QComboBox()
        self.structure.setToolTip("The structure column to draw")
        self.preview = _Preview(self)
        self._glyphs: dict[QToolButton, str] = {}
        self.first_button = self._button("first", "First (Home)", self.first)
        self.prev_button = self._button("previous", "Previous (Left)", self.previous)
        self.next_button = self._button("next", "Next (Right)", self.next)
        self.last_button = self._button("last", "Last (End)", self.last)
        self.position = QSpinBox()
        self.position.setKeyboardTracking(False)
        self.position.setMinimum(1)
        self.position.setButtonSymbols(QSpinBox.ButtonSymbols.NoButtons)
        self.position.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.position.setToolTip("Row number in this list; type one to jump to it")
        self.count_label = QLabel("of 0")
        self.count_label.setEnabled(False)  # the muted colour of a hint
        self.select_box = QCheckBox("Selected")
        self.select_box.setToolTip("Add or remove this row from the selection (Space)")
        self.properties = QTableWidget(0, 2)
        self.properties.setObjectName("BrowserProperties")
        self.properties.horizontalHeader().hide()
        self.properties.verticalHeader().hide()
        self.properties.verticalHeader().setDefaultSectionSize(self.fontMetrics().height() + 8)
        self.properties.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.properties.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.properties.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self.properties.horizontalHeader().setStretchLastSection(True)
        self.properties.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.properties.setWordWrap(False)
        self.properties.setAlternatingRowColors(True)
        self.properties.setShowGrid(False)
        self.properties.horizontalHeader().setDefaultAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )

        top = QHBoxLayout()
        top.setSpacing(8)
        top.addWidget(self.scope, 1)
        top.addWidget(self.structure, 1)
        nav = QHBoxLayout()
        nav.setSpacing(2)
        nav.addWidget(self.first_button)
        nav.addWidget(self.prev_button)
        nav.addSpacing(4)
        nav.addWidget(self.position)
        nav.addSpacing(4)
        nav.addWidget(self.count_label)
        nav.addSpacing(4)
        nav.addWidget(self.next_button)
        nav.addWidget(self.last_button)
        nav.addStretch(1)
        nav.addWidget(self.select_box)
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.setHandleWidth(6)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self.preview)
        splitter.addWidget(self.properties)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 12)
        layout.setSpacing(8)
        layout.addLayout(top)
        layout.addLayout(nav)
        layout.addWidget(splitter, 1)
        self.setMinimumWidth(280)
        self._restyle()
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.setInterval(40)
        self._refresh_timer.timeout.connect(self.refresh)
        self.scope.currentIndexChanged.connect(lambda _i: self.refresh(restart=True))
        self.structure.currentIndexChanged.connect(lambda _i: self._draw())
        self.position.valueChanged.connect(lambda value: self.go(value - 1))
        self.select_box.clicked.connect(self.set_selected)
        self.view.selectionModel().currentChanged.connect(self._on_table_current)
        self._unsubscribe_table = self.table.subscribe(self._on_table_event)
        self._unsubscribe_selection = selection.subscribe(lambda *_: self._schedule())
        self.refresh(restart=True)

    def _button(self, name: str, tip: str, slot) -> QToolButton:
        button = QToolButton()
        button.setAutoRaise(True)
        button.setFixedSize(24, 24)
        button.setToolTip(tip)
        button.clicked.connect(lambda _c=False: slot())
        self._glyphs[button] = name
        return button

    def _restyle(self) -> None:
        """Glyphs, frames and muted text in the palette's colours (again on theme changes)."""
        from .glyphs import flat_buttons, glyph

        palette = self.palette()
        ink = palette.text().color()
        for button, name in self._glyphs.items():
            button.setIcon(glyph(name, ink))
        base = palette.base().color().name()
        line = palette.mid().color().name()
        faint = palette.placeholderText().color()
        muted = f"rgba({faint.red()}, {faint.green()}, {faint.blue()}, {faint.alphaF():.2f})"
        frame = f"background: {base}; border: 1px solid {line}; border-radius: 6px;"
        self.setStyleSheet(
            flat_buttons(palette)
            + f" QLabel#BrowserPreview {{ {frame} color: {muted}; }}"
            + f" QTableWidget#BrowserProperties {{ {frame} padding: 2px; }}"
        )

    def close_browser(self) -> None:
        self._unsubscribe_table()
        self._unsubscribe_selection()
        self._refresh_timer.stop()

    # -- state ----------------------------------------------------------------------------------
    @property
    def oids(self) -> np.ndarray:
        """The browsed rows, in display order."""
        return self._oids

    @property
    def index(self) -> int:
        return self._index

    @property
    def current_oid(self) -> int | None:
        return self._oid

    def scope_key(self) -> str:
        return self.scope.currentData()

    def set_scope(self, key: str) -> None:
        self.scope.setCurrentIndex(self.scope.findData(key))

    def _scope_oids(self) -> np.ndarray:
        shown = self.model.visible_oids
        if self.scope_key() == "shown":
            return shown
        return shown[np.isin(shown, self.selection.as_array())]

    def _schedule(self) -> None:
        if self.isVisible():
            self._refresh_timer.start()
        else:
            self._stale = True  # catch up when shown

    def showEvent(self, event) -> None:
        super().showEvent(event)
        if self._stale:
            self.refresh()

    def _on_table_event(self, event: object) -> None:
        if isinstance(event, Changed):
            self._schedule()

    def refresh(self, *, restart: bool = False) -> None:
        """Recompute the browsed rows, staying on the same row when it is still there."""
        self._refresh_timer.stop()
        self._stale = False
        self._fill_structure_columns()
        self._oids = self._scope_oids()
        n = len(self._oids)
        if restart:
            index = self._start_index()
        elif self._oid is not None and n:
            hit = np.flatnonzero(self._oids == self._oid)
            # A row that left (e.g. deselected while browsing selected rows): stay at its place.
            index = int(hit[0]) if len(hit) else min(max(self._index, 0), n - 1)
        else:
            index = 0 if n else -1
        self._set_index(index if n else -1, move_table=restart)

    def _start_index(self) -> int:
        """Begin at the table's current row when it is in scope, else at the first row."""
        if not len(self._oids):
            return -1
        current = self.view.currentIndex()
        if current.isValid():
            hit = np.flatnonzero(self._oids == self.model.oid_at(current.row()))
            if len(hit):
                return int(hit[0])
        return 0

    def _fill_structure_columns(self) -> None:
        names = [n for n in self.table.column_names if self.table.kind(n) == STRUCTURE]
        current = self.structure.currentText()
        if names != [self.structure.itemText(i) for i in range(self.structure.count())]:
            self.structure.blockSignals(True)
            self.structure.clear()
            self.structure.addItems(names)
            if current in names:
                self.structure.setCurrentText(current)
            self.structure.blockSignals(False)
        self.structure.setVisible(len(names) > 1)

    # -- moving ---------------------------------------------------------------------------------
    def go(self, index: int) -> None:
        if len(self._oids):
            self._set_index(min(max(index, 0), len(self._oids) - 1), move_table=True)

    def first(self) -> None:
        self.go(0)

    def last(self) -> None:
        self.go(len(self._oids) - 1)

    def next(self) -> None:
        self.go(self._index + 1)

    def previous(self) -> None:
        self.go(self._index - 1)

    def _set_index(self, index: int, *, move_table: bool) -> None:
        self._index = index
        self._oid = int(self._oids[index]) if index >= 0 else None
        n = len(self._oids)
        self.position.blockSignals(True)
        self.position.setMaximum(max(n, 1))
        self.position.setValue(index + 1 if index >= 0 else 1)
        self.position.blockSignals(False)
        self.position.setEnabled(n > 0)
        digits = len(f"{max(n, 1):,}")
        self.position.setFixedWidth(self.fontMetrics().horizontalAdvance("0" * digits) + 20)
        self.count_label.setText(f"of {n:,}")
        self.first_button.setEnabled(index > 0)
        self.prev_button.setEnabled(index > 0)
        self.next_button.setEnabled(0 <= index < n - 1)
        self.last_button.setEnabled(0 <= index < n - 1)
        self.select_box.setEnabled(self._oid is not None)
        self.select_box.setChecked(self._oid is not None and self._oid in self.selection)
        if move_table and self._oid is not None:
            self._move_table_cursor(self._oid)
        self._draw()
        self._fill_properties()

    def _move_table_cursor(self, oid: int) -> None:
        row = int(self.model.rows_of_oids([oid])[0])
        if row < 0:
            return
        current = self.view.currentIndex()
        column = current.column() if current.isValid() else 0
        target = self.model.index(row, column)
        self._following = True
        try:
            self.view.selectionModel().setCurrentIndex(
                target, QItemSelectionModel.SelectionFlag.NoUpdate
            )
        finally:
            self._following = False
        self.view.scrollTo(target, QAbstractItemView.ScrollHint.EnsureVisible)

    def _on_table_current(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if self._following or not current.isValid() or not len(self._oids):
            return
        hit = np.flatnonzero(self._oids == self.model.oid_at(current.row()))
        if len(hit) and int(hit[0]) != self._index:
            self._set_index(int(hit[0]), move_table=False)

    def set_selected(self, on: bool) -> None:
        if self._oid is None:
            return
        if on:
            self.selection.add([self._oid], ORIGIN)
        else:
            self.selection.remove([self._oid], ORIGIN)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        key = event.key()
        moves = {
            Qt.Key.Key_Home: self.first,
            Qt.Key.Key_End: self.last,
            Qt.Key.Key_Left: self.previous,
            Qt.Key.Key_Up: self.previous,
            Qt.Key.Key_Right: self.next,
            Qt.Key.Key_Down: self.next,
        }
        if key in moves:
            moves[key]()
        elif key == Qt.Key.Key_Space and self._oid is not None:
            self.set_selected(self._oid not in self.selection)
        else:
            super().keyPressEvent(event)

    # -- showing ----------------------------------------------------------------------------------
    def _empty_text(self) -> str:
        if not len(self.table):
            return "The table is empty"
        return "No rows selected" if self.scope_key() == "selected" else "No rows shown"

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "preview"):
            self._restyle()
            self._draw()  # drawings are coloured for a light or dark background
            self._fill_properties()  # names are drawn in the muted colour

    def _draw(self) -> None:
        self.preview.setPixmap(QPixmap())
        if self._oid is None:
            self.preview.setText(self._empty_text())
            return
        name = self.structure.currentText()
        if not name or name not in self.table:
            self.preview.setText("No structure column")
            return
        column = self.table.column(name)
        mol = column.mol(int(self.table.slots_of([self._oid])[0]))
        if mol is None:
            self.preview.setText("No structure")
            return
        dpr = self.devicePixelRatioF()
        size = self.preview.size()
        dark = self.palette().window().color().lightness() < 128
        png = render_png(
            mol,
            int(size.width() * dpr),
            int(size.height() * dpr),
            scale=dpr * PREVIEW_SCALE,
            dark=dark,
        )
        pixmap = QPixmap()
        pixmap.loadFromData(png, "PNG")
        pixmap.setDevicePixelRatio(dpr)
        self.preview.setText("")
        self.preview.setPixmap(pixmap)

    def _fill_properties(self) -> None:
        if self._oid is None:
            self.properties.setRowCount(0)
            return
        names = [n for n in self.table.column_names if self.table.kind(n) != STRUCTURE]
        self.properties.setRowCount(len(names))
        muted = self.palette().placeholderText()
        for row, name in enumerate(names):
            value = self.table.cell(name, self._oid)
            for col, text in enumerate((name, value)):
                item = QTableWidgetItem(text)
                item.setToolTip(text)
                if col == 0:
                    item.setForeground(muted)
                self.properties.setItem(row, col, item)

    def property_values(self) -> dict[str, str]:
        """What the property table shows (for tests and copying)."""
        return {
            self.properties.item(r, 0).text(): self.properties.item(r, 1).text()
            for r in range(self.properties.rowCount())
        }
