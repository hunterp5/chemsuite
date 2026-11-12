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
"""The Log window's Messages tab: the session's log records, filtered by level, text and job."""

from __future__ import annotations

import logging
import time
from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QPersistentModelIndex, Qt, QTimer
from PySide6.QtGui import QColor, QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QStyle,
    QTableView,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.logbook import Entry, Logbook, logbook

from .error_report import show_log_folder

POLL_MS = 500
LEVEL_CHOICES = (
    ("Debug and above", logging.DEBUG),
    ("Info and above", logging.INFO),
    ("Warnings and errors", logging.WARNING),
    ("Errors only", logging.ERROR),
)
COLOURS = {
    logging.DEBUG: QColor(130, 130, 130),
    logging.WARNING: QColor(205, 120, 0),
    logging.ERROR: QColor(210, 40, 40),
    logging.CRITICAL: QColor(210, 40, 40),
}
Index = QModelIndex | QPersistentModelIndex


def _level_colour(level: int) -> QColor | None:
    if level >= logging.ERROR:
        return COLOURS[logging.ERROR]
    if level >= logging.WARNING:
        return COLOURS[logging.WARNING]
    if level < logging.INFO:
        return COLOURS[logging.DEBUG]
    return None


class LogModel(QAbstractTableModel):
    """The logbook's entries that pass the level and text filters, oldest first."""

    COLUMNS = ("Time", "Level", "Source", "Message")

    def __init__(self, book: Logbook) -> None:
        super().__init__()
        self.book = book
        self.min_level = logging.INFO
        self.needle = ""
        self.job: int | None = None  # only this job's entries
        self.rows: list[Entry] = []
        self._last = 0  # newest seq looked at
        self._cleared = book.cleared_through

    # -- filtering ---------------------------------------------------------------------------------
    def passes(self, entry: Entry) -> bool:
        if entry.level < self.min_level:
            return False
        if self.job is not None and entry.job != self.job:
            return False
        if not self.needle:
            return True
        text = f"{entry.source} {entry.message} {entry.details}".casefold()
        return self.needle in text

    def set_filter(self, min_level: int, needle: str) -> None:
        self.min_level, self.needle = min_level, needle.strip().casefold()
        self.beginResetModel()
        entries = self.book.entries()
        self.rows = [e for e in entries if self.passes(e)]
        self._last = entries[-1].seq if entries else self._last
        self.endResetModel()

    def poll(self) -> int:
        """Take in new entries (and drop ones the logbook let go); returns how many arrived."""
        if self.book.cleared_through != self._cleared:  # Clear, from this or another window
            self._cleared = self.book.cleared_through
            self.set_filter(self.min_level, self.needle)
            return 0
        first = self.book.first_seq
        stale = next((i for i, e in enumerate(self.rows) if e.seq >= first), len(self.rows))
        if stale:
            self.beginRemoveRows(QModelIndex(), 0, stale - 1)
            del self.rows[:stale]
            self.endRemoveRows()
        new = self.book.entries(after=self._last)
        if not new:
            return 0
        self._last = new[-1].seq
        shown = [e for e in new if self.passes(e)]
        if shown:
            start = len(self.rows)
            self.beginInsertRows(QModelIndex(), start, start + len(shown) - 1)
            self.rows.extend(shown)
            self.endInsertRows()
        return len(shown)

    # -- model -------------------------------------------------------------------------------------
    def rowCount(self, parent: Index = QModelIndex()) -> int:  # noqa: B008
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent: Index = QModelIndex()) -> int:  # noqa: B008
        return 0 if parent.isValid() else len(self.COLUMNS)

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = 0):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            return self.COLUMNS[section]
        return None

    def data(self, index: Index, role: int = 0):
        if not index.isValid():
            return None
        entry = self.rows[index.row()]
        column = index.column()
        if role == Qt.ItemDataRole.DisplayRole:
            if column == 0:
                return time.strftime("%H:%M:%S", time.localtime(entry.time))
            if column == 1:
                return entry.level_name.title()
            if column == 2:
                return entry.source + (f" · job {entry.job}" if entry.job is not None else "")
            first = entry.message.splitlines()[0] if entry.message else ""
            return first + (" …" if entry.details else "")
        if role == Qt.ItemDataRole.ForegroundRole:
            return _level_colour(entry.level)
        if role == Qt.ItemDataRole.ToolTipRole and column == 0:
            return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(entry.time))
        return None

    def entry(self, row: int) -> Entry:
        return self.rows[row]


class LogPanel(QWidget):
    """The log, live (the Log window's Messages tab). *show_job(job_id)* shows a job's row
    in the Jobs tab."""

    def __init__(
        self,
        book: Logbook | None = None,
        show_job: Callable[[int], None] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.book = book or logbook()
        self._show_job = show_job
        self.model = LogModel(self.book)
        self.level = QComboBox()
        for title, level in LEVEL_CHOICES:
            self.level.addItem(title, level)
        self.level.setCurrentIndex(1)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Filter messages")
        self.search.setClearButtonEnabled(True)
        self.follow = QCheckBox("Follow")
        self.follow.setChecked(True)
        self.follow.setToolTip("Keep the newest entry in view")
        self.job_filter = QToolButton()
        self.job_filter.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.job_filter.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogCloseButton)
        )
        self.job_filter.setToolTip("Show every job's messages")
        self.job_filter.clicked.connect(lambda: self.show_job_entries(None))
        self.job_filter.hide()
        self.view = QTableView()
        self.view.setModel(self.model)
        self.view.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.view.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.view.verticalHeader().hide()
        self.view.verticalHeader().setDefaultSectionSize(self.fontMetrics().height() + 8)
        self.view.setWordWrap(False)
        self.view.setAlternatingRowColors(True)
        self.view.setShowGrid(False)
        header = self.view.horizontalHeader()
        for column in range(3):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        header.setStretchLastSection(True)
        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setPlaceholderText("Select an entry to see all of it.")
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(self.view)
        splitter.addWidget(self.details)
        splitter.setSizes([400, 140])
        copy = QPushButton("Copy")
        copy.setToolTip("Copy the selected entries (all shown entries if none are selected)")
        copy.clicked.connect(self.copy)
        save = QPushButton("Save As…")
        save.clicked.connect(self.save_as)
        clear = QPushButton("Clear")
        clear.clicked.connect(self.book.clear)
        folder = QPushButton("Show Log Folder")
        folder.clicked.connect(show_log_folder)
        self.count = QLabel()
        top = QHBoxLayout()
        top.addWidget(self.level)
        top.addWidget(self.job_filter)
        top.addWidget(self.search, 1)
        top.addWidget(self.follow)
        bottom = QHBoxLayout()
        bottom.addWidget(self.count, 1)
        for button in (copy, save, clear):
            bottom.addWidget(button)
        bottom.addSpacing(12)
        bottom.addWidget(folder)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        layout.addLayout(top)
        layout.addWidget(splitter, 1)
        layout.addLayout(bottom)
        self.level.currentIndexChanged.connect(self.refilter)
        self.search.textChanged.connect(self.refilter)
        self.view.selectionModel().selectionChanged.connect(self._show_details)
        self.view.doubleClicked.connect(self._open_job)
        self._timer = QTimer(self)
        self._timer.setInterval(POLL_MS)
        self._timer.timeout.connect(self.poll)
        self.refilter()

    def refilter(self, *_args) -> None:
        self.model.set_filter(self.level.currentData(), self.search.text())
        self._after_change()

    def show_job_entries(self, job: int | None, title: str = "") -> None:
        """Show only *job*'s entries (None: every entry's)."""
        self.model.job = job
        self.job_filter.setText(title or f"Job {job}")
        self.job_filter.setVisible(job is not None)
        if job is not None and self.level.currentData() > logging.DEBUG:
            self.level.setCurrentIndex(0)  # a job's story includes its debug lines
        self.refilter()

    def poll(self) -> None:
        if self.model.poll():
            self._after_change()

    def _after_change(self) -> None:
        self.count.setText(f"{self.model.rowCount():,} of {len(self.book):,} entries")
        if self.follow.isChecked():
            self.view.scrollToBottom()

    def selected_entries(self) -> list[Entry]:
        rows = sorted({i.row() for i in self.view.selectionModel().selectedRows()})
        return [self.model.entry(r) for r in rows]

    def _show_details(self, *_args) -> None:
        entries = self.selected_entries()
        self.details.setPlainText(entries[-1].as_text() if entries else "")

    def _open_job(self, index: QModelIndex) -> None:
        entry = self.model.entry(index.row())
        if entry.job is not None and self._show_job is not None:
            self._show_job(entry.job)

    def text(self, entries: list[Entry] | None = None) -> str:
        entries = self.model.rows if entries is None else entries
        return "\n".join(e.as_text() for e in entries)

    def copy(self) -> str:
        text = self.text(self.selected_entries() or None)
        QGuiApplication.clipboard().setText(text)
        return text

    def save_as(self) -> None:
        stamp = time.strftime("%Y%m%d-%H%M%S")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Log", f"chemsuite-log-{stamp}.txt", "Text (*.txt)"
        )
        if path:
            Path(path).write_text(self.text() + "\n", encoding="utf-8")

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self.poll()
        self._timer.start()

    def hideEvent(self, event) -> None:
        super().hideEvent(event)
        self._timer.stop()
