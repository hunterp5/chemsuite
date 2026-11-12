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
"""Running a window's work as jobs, the same way in every kind of window.

Work runs in the background (``JobRunner``); its progress shows in the status bar and the Jobs
window. A failure says so without blocking (details are one click away); long work that ends
while the window is not in front alerts it and posts a notification; cancelled per-row work
can keep what it finished. A *guarded* job's result is dropped when the window has moved on to
another document meanwhile.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any

import shiboken6
from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication, QProgressDialog, QWidget

from .error_report import details_of, show_error
from .jobs import QUEUE, QUICK, Job, JobRunner
from .jobs_ui import duration
from .settings_hub import settings_hub

log = logging.getLogger(__name__)

NOTIFY_AFTER_S = 30  # a background job this long gets a "done" notice

#: run(label, work, on_done, ...) as a callable: what panels and searchers are handed
Run = Callable[..., Job]


class JobHost:
    """The jobs of one window.

    *current* names what results apply to (the window's document): a guarded job whose
    *current* changed meanwhile is not applied. *name* is the document's name, shown with each
    job; *message* shows a line in the window's status bar.
    """

    def __init__(
        self,
        window: QWidget,
        runner: JobRunner,
        *,
        current: Callable[[], object],
        name: Callable[[], str],
        message: Callable[[str, int], None],
    ) -> None:
        self.window = window
        self.runner = runner
        self._current = current
        self._name = name
        self.message = message
        #: the last failure's box (kept so it is not deleted while shown)
        self.failure_box: QWidget | None = None

    def run(
        self,
        label: str,
        work: Callable[[Job], Any],
        on_done: Callable[[Any], None],
        *,
        wait: bool = False,
        lane: str = QUICK,
        guard: bool = False,
        keep_partial: bool = False,
    ) -> Job:
        """Run *work* in the background; *on_done* gets its result on the GUI thread.

        *wait* blocks (with a progress dialog) until the job ends, for callers that need the
        outcome now, such as Save from the close prompt. With *guard*, the result is dropped
        if the window has moved on to another document meanwhile. With *keep_partial*, a
        cancelled job's ``Cancelled(partial)`` still reaches *on_done*.
        """
        current = self._current()
        job = self.runner.submit(label, work, lane=lane, owner=self.window, document=self._name())

        def finished(result: Any) -> None:
            if not self.alive():
                log.info("dropping the result of %r: its window is closed", label)
                job.set_summary("Not applied: the window was closed")
                return
            if guard and self._current() is not current:
                log.info("dropping the result of %r: its document is gone", label)
                job.set_summary("Not applied: the document was closed")
                return
            on_done(result)
            if job.lane == QUEUE and job.elapsed > NOTIFY_AFTER_S:
                self.message(f"{label.rstrip('…')}: done", 8000)
                self.notify_unseen(job, job.summary or "Done")

        job.finished.connect(finished)
        if lane == QUEUE:
            job.failed.connect(lambda err: self.alive() and self.report_failure(job, err))
        else:
            job.failed.connect(lambda err: self.alive() and self.show_error(label, err))

        def cancelled() -> None:
            if not self.alive():  # cancelled because its window closed
                return
            name = label.rstrip("…")
            if keep_partial and job.partial_result is not None:
                if guard and self._current() is not current:
                    job.set_summary("Cancelled; partial results not applied: document closed")
                    return
                on_done(job.partial_result)  # writes what was finished before the cancel
                kept = job.summary or "kept the finished rows"
                job.set_summary(f"Cancelled, partial results kept: {kept}")
                self.message(f"Cancelled: {name}; partial results kept", 6000)
            else:
                self.message(f"Cancelled: {name}", 3000)

        job.cancelled.connect(cancelled)
        if wait:
            self.wait_for(job)
        return job

    def alive(self) -> bool:
        """Whether the window is still there: a job's outcome can arrive after it closed (the
        close cancels the window's jobs, and the cancel is delivered later)."""
        return shiboken6.isValid(self.window)

    def wait_for(self, job: Job) -> None:
        """Block, with a progress dialog that can cancel, until *job* ends."""
        dialog = QProgressDialog(job.title, "Cancel", 0, 0, self.window)
        dialog.setWindowTitle("chemsuite")
        dialog.setMinimumDuration(400)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)
        dialog.setModal(True)
        timer = QTimer(dialog)

        def poll() -> None:
            done, total = job.progress
            if total:
                dialog.setMaximum(1000)
                dialog.setValue(int(1000 * min(done, total) / total))

        timer.timeout.connect(poll)
        timer.start(100)
        dialog.canceled.connect(job.cancel)
        loop = QEventLoop()
        job.done.connect(timer.stop)
        job.done.connect(dialog.close)
        job.done.connect(loop.quit)
        if job.running:  # `done` is delivered on this thread, so it cannot be missed
            loop.exec()
        dialog.deleteLater()

    def notify_unseen(self, job: Job, outcome: str) -> None:
        """A long job ended while the window was not in front: alert the dock or taskbar and,
        unless Settings ▸ General turns it off, post a desktop notification."""
        if self.window.isActiveWindow():
            return
        QApplication.alert(self.window)
        if settings_hub().settings.jobs.notify:
            from .notify import notify

            where = f" ({job.document})" if job.document else ""
            notify(f"{job.name}{where}", f"{outcome} after {duration(job.elapsed)}")

    def report_failure(self, job: Job, error: BaseException) -> None:
        """A background job failed: say so without blocking; details are one click away."""
        if job.elapsed > NOTIFY_AFTER_S:
            self.notify_unseen(job, f"Failed: {error}" if str(error) else "Failed")
        self.message(f"{job.title.rstrip('…')} failed: {error}", 10000)
        self.failure_box = show_error(
            self.window,
            f"{job.title.rstrip('…')} failed.",
            str(error) or type(error).__name__,
            job.error_text,
        )

    def show_error(self, label: str, error: BaseException) -> None:
        self.failure_box = show_error(
            self.window,
            f"{label.rstrip('…')} failed.",
            str(error) or type(error).__name__,
            details_of(error),
        )
