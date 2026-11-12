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
"""The jobs indicator in each window's status bar, and the Log window's Jobs tab."""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QItemSelectionModel, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QGuiApplication, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QHeaderView,
    QMenu,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .jobs import CANCELLED, FAILED, FINISHED, QUEUE, QUEUED, RUNNING, Job, JobRunner

STATE_TEXT = {
    QUEUED: "Queued",
    RUNNING: "Running",
    FINISHED: "Done",
    FAILED: "Failed",
    CANCELLED: "Cancelled",
}
FAILED_COLOUR = QColor(200, 40, 40)
TICK_MS = 250


def duration(seconds: float) -> str:
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds} s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}:{seconds:02d}"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}:{minutes:02d}:{seconds:02d}"


def fraction(job: Job) -> float | None:
    done, total = job.progress
    return min(done, total) / total if total else None


#: text, not a raised button: it sits in the slim status bar like the other labels
BUTTON_STYLE = "QToolButton { border: none; padding: 0px 2px; }"


def progress_text(job: Job) -> str:
    """What a job is doing, e.g. ``Generating conformers… (120/500)``."""
    name = job.name
    if job.cancelling:
        return f"Cancelling {name[:1].lower() + name[1:]}…"
    done, total = job.progress
    text = f"{name}…"
    return f"{text} ({min(done, total):,}/{total:,})" if total else text


class JobsIndicator(QWidget):
    """The running job's progress as a sentence, plus what else is running, queued or failed:
    ``Generating conformers… (120/500) · 3 queued``. Click to open the Log window's jobs.

    Failures stay flagged until the jobs have been seen.
    """

    clicked = Signal()
    #: after each refresh (the Log button follows it)
    refreshed = Signal()

    def __init__(self, runner: JobRunner, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.runner = runner
        self._seen_failures: set[int] = set()
        self.button = QToolButton()
        self.button.setAutoRaise(True)
        self.button.setStyleSheet(BUTTON_STYLE)
        self.button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.button.setToolTip("Show the jobs")
        self.button.clicked.connect(self.clicked)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.button)
        self._timer = QTimer(self)
        self._timer.setInterval(TICK_MS)
        self._timer.timeout.connect(self.refresh)
        runner.changed.connect(self.refresh)
        self.refresh()

    def acknowledge(self) -> None:
        """The user has seen the failures (the jobs were shown)."""
        self._seen_failures |= {j.id for j in self.runner.jobs if j.state == FAILED}
        self.refresh()

    def unseen_failures(self) -> int:
        return sum(j.state == FAILED and j.id not in self._seen_failures for j in self.runner.jobs)

    def current_job(self) -> Job | None:
        """The job the status bar describes: the running queue job, else any running job."""
        running = [j for j in self.runner.jobs if j.state == RUNNING]
        return next((j for j in running if j.lane == QUEUE), running[0] if running else None)

    def refresh(self) -> None:
        jobs = self.runner.jobs
        running = [j for j in jobs if j.state == RUNNING]
        queued = sum(j.state == QUEUED for j in jobs)
        failed = self.unseen_failures()
        current = self.current_job()
        parts = []
        if current is not None:
            parts.append(progress_text(current))
            if len(running) > 1:
                parts.append(f"{len(running) - 1} more running")
        if queued:
            parts.append(f"{queued} queued")
        if failed:
            parts.append(f"{failed} failed")
        self.button.setText(" · ".join(parts))
        self.button.setToolTip(
            f"{current.title} for {duration(current.elapsed)}. Click to show the jobs."
            if current is not None
            else "Show the jobs"
        )
        colour = FAILED_COLOUR.name() if failed and current is None else ""
        self.button.setStyleSheet(
            BUTTON_STYLE + (f" QToolButton {{ color: {colour}; }}" if colour else "")
        )
        self.setVisible(bool(parts))
        if running or queued:
            self._timer.start()
        else:
            self._timer.stop()
        self.refreshed.emit()


class JobsPanel(QWidget):
    """Every job of the session: running, queued and ended, with controls (the Log window's
    Jobs tab). *show_messages(job)* shows a job's log messages."""

    COLUMNS = ("Job", "Document", "Status", "Progress", "Time", "Result")

    def __init__(
        self,
        runner: JobRunner,
        show_messages: Callable[[Job], None] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.runner = runner
        self._show_messages = show_messages
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(list(self.COLUMNS))
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._menu)
        self.table.verticalHeader().hide()
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        self.table.setColumnWidth(0, 260)
        self.table.itemSelectionChanged.connect(self._update_buttons)
        self.table.itemDoubleClicked.connect(lambda _item: self._with_current(self.show_messages))
        self.table.setAlternatingRowColors(True)
        self.table.setShowGrid(False)
        self.table.verticalHeader().setDefaultSectionSize(self.fontMetrics().height() + 8)
        self.cancel_button = QPushButton("Cancel Job")
        self.up_button = QPushButton("Move Up")
        self.down_button = QPushButton("Move Down")
        self.copy_button = QPushButton("Copy Details")
        self.messages_button = QPushButton("Show Messages")
        self.messages_button.setToolTip("Show the selected job's log messages")
        self.messages_button.setVisible(show_messages is not None)
        self.messages_button.clicked.connect(lambda: self._with_current(self.show_messages))
        self.remove_button = QPushButton("Remove")
        self.remove_button.setToolTip("Remove the selected ended jobs from the list (Delete)")
        clear = QPushButton("Clear Finished")
        self.cancel_button.clicked.connect(self.cancel_selected)
        self.remove_button.clicked.connect(self.remove_selected)
        for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
            shortcut = QShortcut(key, self.table)
            shortcut.setContext(Qt.ShortcutContext.WidgetShortcut)
            shortcut.activated.connect(self.remove_selected)
        self.up_button.clicked.connect(lambda: self._with_current(lambda j: runner.move(j, -1)))
        self.down_button.clicked.connect(lambda: self._with_current(lambda j: runner.move(j, 1)))
        self.copy_button.clicked.connect(lambda: self._with_current(self.copy_details))
        clear.clicked.connect(runner.clear_finished)
        self.table.setToolTip(
            "Heavy jobs run one at a time, in order; quick jobs (opening files, panels) run "
            "alongside them."
        )
        buttons = QHBoxLayout()
        for button in (self.cancel_button, self.up_button, self.down_button):
            buttons.addWidget(button)
        buttons.addSpacing(12)
        buttons.addWidget(self.copy_button)
        buttons.addWidget(self.messages_button)
        buttons.addStretch(1)
        buttons.addWidget(self.remove_button)
        buttons.addWidget(clear)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)
        self._timer = QTimer(self)
        self._timer.setInterval(TICK_MS * 2)
        self._timer.timeout.connect(self.refresh)
        runner.changed.connect(self.refresh)
        self.refresh()

    # -- contents ----------------------------------------------------------------------------------
    def ordered(self) -> list[Job]:
        """Running first, then the queue in order, then ended jobs newest first."""
        jobs = self.runner.jobs
        running = [j for j in jobs if j.state == RUNNING]
        queued = self.runner.queued_jobs
        queued += [j for j in jobs if j.state == QUEUED and j not in queued]  # quick, waiting
        ended = [j for j in reversed(jobs) if not j.running]
        return running + queued + ended

    def refresh(self) -> None:
        if not self.isVisible():
            self._timer.stop()
            return
        chosen = {job.id for job in self.selected()}
        jobs = self.ordered()
        self.table.setRowCount(len(jobs))
        for row, job in enumerate(jobs):
            share = fraction(job)
            progress = f"{share:.0%}" if share is not None and job.running else ""
            if job.state == RUNNING and share is None:
                progress = "…"
            elapsed = duration(job.elapsed) if job.started_at is not None else ""
            result = job.summary or (job.error_text.strip().splitlines() or [""])[-1]
            values = (job.title, job.document, STATE_TEXT[job.state], progress, elapsed, result)
            for col, value in enumerate(values):
                item = self.table.item(row, col)
                if item is None:
                    item = QTableWidgetItem()
                    self.table.setItem(row, col, item)
                item.setText(value)
                item.setData(Qt.ItemDataRole.UserRole, job.id)
                item.setForeground(FAILED_COLOUR if job.state == FAILED else self.palette().text())
                item.setToolTip(job.error_text if job.state == FAILED else value)
        selection = self.table.selectionModel()
        flags = QItemSelectionModel.SelectionFlag
        selection.clearSelection()
        for row, job in enumerate(jobs):  # rows move as jobs start and end: follow the jobs
            if job.id in chosen:
                selection.select(self.table.model().index(row, 0), flags.Select | flags.Rows)
        if any(j.running for j in jobs):
            self._timer.start()
        else:
            self._timer.stop()
        self._update_buttons()

    def selected(self) -> list[Job]:
        """The selected jobs, top to bottom."""
        rows = sorted(index.row() for index in self.table.selectionModel().selectedRows())
        ids = [self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) for row in rows]
        by_id = {job.id: job for job in self.runner.jobs}
        return [by_id[i] for i in ids if i in by_id]

    def current(self) -> Job | None:
        """The selected job when exactly one is selected."""
        chosen = self.selected()
        return chosen[0] if len(chosen) == 1 else None

    def cancel_selected(self) -> None:
        for job in self.selected():
            if job.running:
                job.cancel()
        self.refresh()

    def remove_selected(self) -> int:
        """Remove the selected ended jobs from the list; running and queued ones stay."""
        ended = [job for job in self.selected() if not job.running]
        for job in ended:
            self.runner.remove(job)
        self.refresh()
        return len(ended)

    def _menu(self, pos) -> None:
        row = self.table.rowAt(pos.y())
        if row >= 0 and row not in {i.row() for i in self.table.selectionModel().selectedRows()}:
            self.table.selectRow(row)
        jobs = self.selected()
        if not jobs:
            return
        menu = QMenu(self)
        cancel = menu.addAction("Cancel Job")
        cancel.setEnabled(any(job.running for job in jobs))
        remove = menu.addAction("Remove")
        remove.setEnabled(any(not job.running for job in jobs))
        copy = menu.addAction("Copy Details")
        copy.setEnabled(len(jobs) == 1)
        messages = menu.addAction("Show Messages")
        messages.setEnabled(len(jobs) == 1)
        messages.setVisible(self._show_messages is not None)
        chosen = menu.exec(self.table.viewport().mapToGlobal(pos))
        if chosen is cancel:
            self.cancel_selected()
        elif chosen is remove:
            self.remove_selected()
        elif chosen is copy:
            self.copy_details(jobs[0])
        elif chosen is messages:
            self.show_messages(jobs[0])

    def show_messages(self, job: Job) -> None:
        if self._show_messages is not None:
            self._show_messages(job)

    def select(self, job: Job) -> None:
        for row in range(self.table.rowCount()):
            if self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) == job.id:
                self.table.selectRow(row)
                return

    def _with_current(self, action) -> None:
        job = self.current()
        if job is not None:
            action(job)
            self.refresh()

    def _update_buttons(self) -> None:
        job = self.current()
        chosen = self.selected()
        queue = self.runner.queued_jobs
        self.cancel_button.setEnabled(any(j.running for j in chosen))
        self.remove_button.setEnabled(any(not j.running for j in chosen))
        self.up_button.setEnabled(job in queue and queue.index(job) > 0)
        self.down_button.setEnabled(job in queue and queue.index(job) < len(queue) - 1)
        self.copy_button.setEnabled(job is not None)
        self.messages_button.setEnabled(job is not None)

    def copy_details(self, job: Job) -> str:
        lines = [
            job.title,
            f"Document: {job.document or '-'}",
            f"Status: {STATE_TEXT[job.state]}",
        ]
        if job.started_at is not None:
            lines.append(f"Time: {duration(job.elapsed)}")
        if job.summary:
            lines.append(f"Result: {job.summary}")
        if job.error_text:
            lines += ["", job.error_text.rstrip()]
        text = "\n".join(lines)
        QGuiApplication.clipboard().setText(text)
        return text

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self.refresh()
