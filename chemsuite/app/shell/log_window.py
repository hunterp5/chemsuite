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
"""The Log window: the session's jobs and its log messages, as two tabs.

It opens from the **Log** button at the right of the menu bar (Window ▸ Log on macOS) and from
the jobs indicator in the status bar. A job's row shows its messages; a message from a job
shows that job's row.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import QTabWidget, QVBoxLayout, QWidget

from chemsuite.core.logbook import Logbook

from .jobs import FAILED, QUEUED, RUNNING, Job, JobRunner
from .jobs_ui import FAILED_COLOUR, JobsPanel
from .log_ui import LogPanel

JOBS, MESSAGES = 0, 1


class LogWindow(QWidget):
    """One per application: the Jobs and Messages tabs."""

    def __init__(
        self, runner: JobRunner, book: Logbook | None = None, parent: QWidget | None = None
    ) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Log")
        self.runner = runner
        self.jobs = JobsPanel(runner, show_messages=self.show_messages)
        self.messages = LogPanel(book, show_job=self.show_job)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.jobs, "Jobs")
        self.tabs.addTab(self.messages, "Messages")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.addWidget(self.tabs)
        for key in (QKeySequence.StandardKey.Close, QKeySequence(Qt.Key.Key_Escape)):
            QShortcut(key, self, self.close)
        runner.changed.connect(self._label_jobs)
        self._label_jobs()
        self.resize(960, 600)

    def show_tab(self, index: int) -> None:
        self.tabs.setCurrentIndex(index)

    def show_job(self, job_id: int) -> None:
        """The Jobs tab, with *job_id*'s row selected."""
        self.show_tab(JOBS)
        job = next((j for j in self.runner.jobs if j.id == job_id), None)
        if job is not None:
            self.jobs.refresh()
            self.jobs.select(job)

    def show_messages(self, job: Job) -> None:
        """The Messages tab, showing only *job*'s messages."""
        self.messages.show_job_entries(job.id, job.title)
        self.show_tab(MESSAGES)

    def _label_jobs(self) -> None:
        jobs = self.runner.jobs
        active = sum(j.state in (RUNNING, QUEUED) for j in jobs)
        failed = any(j.state == FAILED for j in jobs)
        self.tabs.setTabText(JOBS, f"Jobs ({active})" if active else "Jobs")
        bar = self.tabs.tabBar()
        bar.setTabTextColor(
            JOBS, FAILED_COLOUR if failed and not active else self.palette().text().color()
        )
