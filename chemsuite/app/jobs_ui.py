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
"""The jobs indicator in each window's status bar, and the Jobs window (View ▸ Jobs)."""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QHeaderView,
    QLabel,
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
    ``Generating conformers… (120/500) · 3 queued``. Click to open the Jobs window.

    Failures stay flagged until the Jobs window has been seen.
    """

    clicked = Signal()

    def __init__(self, runner: JobRunner, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.runner = runner
        self._seen_failures: set[int] = set()
        self.button = QToolButton()
        self.button.setAutoRaise(True)
        self.button.setStyleSheet(BUTTON_STYLE)
        self.button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.button.setToolTip("Show the Jobs window")
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
        """The user has seen the failures (the Jobs window opened)."""
        self._seen_failures |= {j.id for j in self.runner.jobs if j.state == FAILED}
        self.refresh()

    def current_job(self) -> Job | None:
        """The job the status bar describes: the running queue job, else any running job."""
        running = [j for j in self.runner.jobs if j.state == RUNNING]
        return next((j for j in running if j.lane == QUEUE), running[0] if running else None)

    def refresh(self) -> None:
        jobs = self.runner.jobs
        running = [j for j in jobs if j.state == RUNNING]
        queued = sum(j.state == QUEUED for j in jobs)
        failed = sum(j.state == FAILED and j.id not in self._seen_failures for j in jobs)
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
            f"{current.title} for {duration(current.elapsed)}. Click for the Jobs window."
            if current is not None
            else "Show the Jobs window"
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


class JobsWindow(QWidget):
    """Every job of the session: running, queued and ended, with controls."""

    COLUMNS = ("Job", "Document", "Status", "Progress", "Time", "Result")

    def __init__(self, runner: JobRunner, parent: QWidget | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Jobs")
        self.runner = runner
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(list(self.COLUMNS))
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.verticalHeader().hide()
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        self.table.setColumnWidth(0, 260)
        self.table.itemSelectionChanged.connect(self._update_buttons)
        self.cancel_button = QPushButton("Cancel Job")
        self.up_button = QPushButton("Move Up")
        self.down_button = QPushButton("Move Down")
        self.copy_button = QPushButton("Copy Details")
        clear = QPushButton("Clear Finished")
        self.cancel_button.clicked.connect(lambda: self._with_current(lambda j: j.cancel()))
        self.up_button.clicked.connect(lambda: self._with_current(lambda j: runner.move(j, -1)))
        self.down_button.clicked.connect(lambda: self._with_current(lambda j: runner.move(j, 1)))
        self.copy_button.clicked.connect(lambda: self._with_current(self.copy_details))
        clear.clicked.connect(runner.clear_finished)
        self.note = QLabel(
            "Heavy jobs run one at a time, in order; quick jobs (opening files, panels) run "
            "alongside them."
        )
        self.note.setWordWrap(True)
        buttons = QHBoxLayout()
        for button in (self.cancel_button, self.up_button, self.down_button, self.copy_button):
            buttons.addWidget(button)
        buttons.addStretch(1)
        buttons.addWidget(clear)
        layout = QVBoxLayout(self)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)
        layout.addWidget(self.note)
        self._timer = QTimer(self)
        self._timer.setInterval(TICK_MS * 2)
        self._timer.timeout.connect(self.refresh)
        runner.changed.connect(self.refresh)
        self.resize(820, 360)
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
        current = self.current()
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
            if current is not None and job is current:
                self.table.selectRow(row)
        if any(j.running for j in jobs):
            self._timer.start()
        else:
            self._timer.stop()
        self._update_buttons()

    def current(self) -> Job | None:
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return None
        job_id = self.table.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        return next((j for j in self.runner.jobs if j.id == job_id), None)

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
        queue = self.runner.queued_jobs
        self.cancel_button.setEnabled(job is not None and job.running)
        self.up_button.setEnabled(job in queue and queue.index(job) > 0)
        self.down_button.setEnabled(job in queue and queue.index(job) < len(queue) - 1)
        self.copy_button.setEnabled(job is not None)

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
