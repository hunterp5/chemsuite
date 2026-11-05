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

"""Running core jobs off the GUI thread.

A job function receives its :class:`Job` and uses ``job.report(done, total)``,
``job.cancel_token`` and ``job.emit_partial(value)``. It runs on a worker thread (and may fan
out to the shared process pool); its result, partial results and errors are delivered back on
the GUI thread as Qt signals. Progress is polled rather than signalled, so a busy job can never
flood the event loop.
"""

from __future__ import annotations

import logging
import threading
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from typing import Any

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal, Slot

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import process_pool

log = logging.getLogger(__name__)


class Job(QObject):
    """Public signals are always emitted on the GUI thread (the thread that owns the job).

    The worker emits the private ``_worker_*`` signals; they arrive at the job's own slots by
    queued connection and are re-emitted from there. Connecting plain callables directly to a
    signal emitted on a worker thread would run them on that thread in PySide6.
    """

    finished = Signal(object)
    failed = Signal(object)  # the exception
    cancelled = Signal()
    partial = Signal(object)
    done = Signal()  # after any of finished / failed / cancelled

    _worker_partial = Signal(object)
    _worker_end = Signal(str, object)

    def __init__(self, title: str, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.title = title
        self.cancel_token = CancelToken()
        self._lock = threading.Lock()
        self._progress = (0, 0)
        self.running = True
        self._worker_partial.connect(self._on_partial)
        self._worker_end.connect(self._on_end)

    # called from the worker thread ----------------------------------------------------------
    def report(self, done: int, total: int) -> None:
        with self._lock:
            self._progress = (int(done), int(total))

    def emit_partial(self, value: Any) -> None:
        self._worker_partial.emit(value)

    # relayed onto the GUI thread ------------------------------------------------------------
    @Slot(object)
    def _on_partial(self, value: Any) -> None:
        if self.running:
            self.partial.emit(value)

    @Slot(str, object)
    def _on_end(self, outcome: str, value: Any) -> None:
        self.running = False
        if outcome == "finished":
            self.finished.emit(value)
        elif outcome == "failed":
            self.failed.emit(value)
        else:
            self.cancelled.emit()
        self.done.emit()

    # called from the GUI thread -------------------------------------------------------------
    @property
    def progress(self) -> tuple[int, int]:
        with self._lock:
            return self._progress

    def cancel(self) -> None:
        self.cancel_token.cancel()


class _Runnable(QRunnable):
    def __init__(self, job: Job, fn: Callable[[Job], Any]) -> None:
        super().__init__()
        self.setAutoDelete(True)
        self._job = job
        self._fn = fn

    def run(self) -> None:
        job = self._job
        try:
            result = self._fn(job)
        except Cancelled:
            job._worker_end.emit("cancelled", None)
        except BaseException as err:  # top-level job handler: logged here, shown by the GUI
            log.exception("job %r failed", job.title)
            job._worker_end.emit("failed", err)
        else:
            job._worker_end.emit("finished", result)


class JobRunner(QObject):
    """Owns the worker threads and the (lazily started) process pool for one application."""

    def __init__(self, parent: QObject | None = None, *, threads: int = 3) -> None:
        super().__init__(parent)
        self._threads = QThreadPool(self)
        self._threads.setMaxThreadCount(threads)
        self._pool: ProcessPoolExecutor | None = None
        self._jobs: set[Job] = set()

    def process_pool(self) -> ProcessPoolExecutor:
        if self._pool is None:
            self._pool = process_pool()
        return self._pool

    def submit(self, title: str, fn: Callable[[Job], Any]) -> Job:
        job = Job(title, self)
        self._jobs.add(job)
        job.done.connect(lambda: self._jobs.discard(job))
        self._threads.start(_Runnable(job, fn))
        return job

    @property
    def active_jobs(self) -> list[Job]:
        return [job for job in self._jobs if job.running]

    def shutdown(self, timeout_ms: int = 5000) -> None:
        for job in list(self._jobs):
            job.cancel()
        self._threads.waitForDone(timeout_ms)
        if self._pool is not None:
            self._pool.shutdown(wait=False, cancel_futures=True)
            self._pool = None
