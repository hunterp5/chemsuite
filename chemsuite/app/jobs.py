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
"""Running core jobs off the GUI thread, in two lanes.

A job function receives its :class:`Job` and uses ``job.report(done, total)``,
``job.cancel_token`` and ``job.emit_partial(value)``. It runs on a worker thread (and may fan
out to the shared process pool); its result, partial results and errors are delivered back on
the GUI thread as Qt signals. Progress is polled rather than signalled, so a busy job can never
flood the event loop.

Lanes:

- ``QUEUE``: heavy work (descriptors, conformers, MMP, maps…) runs one job at a time, in the
  order queued. Queued jobs can be cancelled or moved before they start, and the computer is
  kept awake while any is waiting or running.
- ``QUICK``: short or interactive work (opening files, a panel's statistics) starts at once on
  a small thread pool, never waiting behind the queue.
"""

from __future__ import annotations

import itertools
import logging
import threading
import time
import traceback
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from typing import Any

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal, Slot

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.awake import KeepAwake
from chemsuite.core.jobs.parallel import process_pool

log = logging.getLogger(__name__)

QUEUE, QUICK = "queue", "quick"
QUEUED, RUNNING, FINISHED, FAILED, CANCELLED = (
    "queued",
    "running",
    "finished",
    "failed",
    "cancelled",
)
ENDED = (FINISHED, FAILED, CANCELLED)
HISTORY_LIMIT = 200  # ended jobs kept for the Jobs window

_ids = itertools.count(1)


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
    changed = Signal()  # state, summary or position changed

    _worker_started = Signal()
    _worker_partial = Signal(object)
    _worker_end = Signal(str, object, object)  # outcome, value, the work function (carried)

    def __init__(
        self,
        title: str,
        parent: QObject | None = None,
        *,
        lane: str = QUICK,
        owner: object = None,
        document: str = "",
    ) -> None:
        super().__init__(parent)
        self.id = next(_ids)
        self.title = title
        self.lane = lane
        self.owner = owner  # what the job belongs to (a window); its jobs end with it
        self.document = document
        self.state = QUEUED
        self.queued_at = time.monotonic()
        self.started_at: float | None = None
        self.ended_at: float | None = None
        self.summary = ""  # one line about the result ("Added 3 columns")
        self.error_text = ""  # the traceback of a failure
        self.cancel_token = CancelToken()
        self._lock = threading.Lock()
        self._progress = (0, 0)
        self._worker_started.connect(self._on_started)
        self._worker_partial.connect(self._on_partial)
        self._worker_end.connect(self._on_end)

    # called from the worker thread ----------------------------------------------------------
    def report(self, done: int, total: int) -> None:
        with self._lock:
            self._progress = (int(done), int(total))

    def emit_partial(self, value: Any) -> None:
        self._worker_partial.emit(value)

    # relayed onto the GUI thread ------------------------------------------------------------
    @Slot()
    def _on_started(self) -> None:
        if self.state == QUEUED:
            self.state = RUNNING
            self.started_at = time.monotonic()
            log.info("job %d started: %s", self.id, self.title)
            self.changed.emit()

    @Slot(object)
    def _on_partial(self, value: Any) -> None:
        if self.running:
            self.partial.emit(value)

    @Slot(str, object, object)
    def _on_end(self, outcome: str, value: Any, carried: list | None = None) -> None:
        if carried:
            carried.clear()  # release the work function here, on the GUI thread (see _Runnable)
        if self.state in ENDED:
            return
        self.state = outcome
        self.ended_at = time.monotonic()
        if outcome == FINISHED:
            log.info("job %d finished in %.1f s: %s", self.id, self.elapsed, self.title)
            self.finished.emit(value)
        elif outcome == FAILED:
            log.error("job %d failed: %s\n%s", self.id, self.title, self.error_text)
            self.failed.emit(value)
        else:
            log.info("job %d cancelled: %s", self.id, self.title)
            self.cancelled.emit()
        self.done.emit()
        self.changed.emit()

    # called from the GUI thread -------------------------------------------------------------
    @property
    def running(self) -> bool:
        """Not ended yet: queued or running."""
        return self.state not in ENDED

    @property
    def progress(self) -> tuple[int, int]:
        with self._lock:
            return self._progress

    @property
    def elapsed(self) -> float:
        """Seconds spent running (0 while queued)."""
        if self.started_at is None:
            return 0.0
        return (self.ended_at or time.monotonic()) - self.started_at

    def set_summary(self, text: str) -> None:
        self.summary = text
        self.changed.emit()

    def cancel(self) -> None:
        self.cancel_token.cancel()
        runner = self.parent()
        if self.state == QUEUED and isinstance(runner, JobRunner):
            runner._drop_queued(self)


class _Runnable(QRunnable):
    def __init__(self, job: Job, fn: Callable[[Job], Any]) -> None:
        super().__init__()
        self.setAutoDelete(True)
        self._job = job
        self._fn = fn

    def run(self) -> None:
        # The work function's closure may hold the last reference to a window or dialog. Python
        # would destroy it on this thread when the function is dropped, and Qt objects must die
        # on the GUI thread. So the function travels back in a list the GUI thread empties.
        job = self._job
        carried = [self._fn]
        self._fn = None
        outcome, value = CANCELLED, None
        if not job.cancel_token.cancelled:
            job._worker_started.emit()
            try:
                value = carried[0](job)
                outcome = FINISHED
            except Cancelled:
                value = None
            except BaseException as err:  # noqa: BLE001 - top-level job handler; logged on end
                job.error_text = "".join(traceback.format_exception(err))
                outcome, value = FAILED, err
        job._worker_end.emit(outcome, value, carried)


class JobRunner(QObject):
    """Owns the worker threads, the job queue and the (lazily started) process pool."""

    changed = Signal()  # a job was added, removed, moved or changed state

    def __init__(
        self, parent: QObject | None = None, *, threads: int = 3, keep_awake: bool = False
    ) -> None:
        super().__init__(parent)
        self._threads = QThreadPool(self)
        self._threads.setMaxThreadCount(threads)
        self._heavy = QThreadPool(self)  # the queue's one worker
        self._heavy.setMaxThreadCount(1)
        self._pool: ProcessPoolExecutor | None = None
        self._jobs: list[Job] = []  # every job not cleared, in the order submitted
        self._waiting: list[tuple[Job, Callable[[Job], Any]]] = []
        self._current: Job | None = None
        self._awake = KeepAwake() if keep_awake else None

    def process_pool(self) -> ProcessPoolExecutor:
        if self._pool is None:
            from .settings_hub import settings_hub

            # Settings ▸ General ▸ Processor cores (0: automatic); read once, at first use.
            self._pool = process_pool(settings_hub().settings.performance.max_workers or None)
        return self._pool

    # -- submitting ------------------------------------------------------------------------------
    def submit(
        self,
        title: str,
        fn: Callable[[Job], Any],
        *,
        lane: str = QUICK,
        owner: object = None,
        document: str = "",
    ) -> Job:
        job = Job(title, self, lane=lane, owner=owner, document=document)
        self._jobs.append(job)
        job.changed.connect(self.changed)
        job.done.connect(lambda: self._ended(job))
        log.info("job %d queued (%s): %s", job.id, lane, title)
        if lane == QUICK:
            self._threads.start(_Runnable(job, fn))
        else:
            self._waiting.append((job, fn))
            self._pump()
        self._trim()
        self._update_awake()
        self.changed.emit()
        return job

    def _pump(self) -> None:
        if self._current is None and self._waiting:
            job, fn = self._waiting.pop(0)
            self._current = job
            self._heavy.start(_Runnable(job, fn))

    def _ended(self, job: Job) -> None:
        if job is self._current:
            self._current = None
            self._pump()
        self._update_awake()
        self.changed.emit()

    def _drop_queued(self, job: Job) -> None:
        """Cancel a job that has not started: it ends now, without running."""
        entry = next((e for e in self._waiting if e[0] is job), None)
        if entry is not None:
            self._waiting.remove(entry)
            job._on_end(CANCELLED, None, None)

    def _update_awake(self) -> None:
        if self._awake is not None:
            self._awake.hold(self._current is not None or bool(self._waiting))

    def _trim(self) -> None:
        ended = [j for j in self._jobs if not j.running]
        for job in ended[: max(0, len(ended) - HISTORY_LIMIT)]:
            self._jobs.remove(job)

    # -- looking and managing --------------------------------------------------------------------
    @property
    def jobs(self) -> list[Job]:
        """Every job still listed: running, queued (in queue order), then ended ones."""
        return list(self._jobs)

    @property
    def queued_jobs(self) -> list[Job]:
        return [job for job, _fn in self._waiting]

    @property
    def active_jobs(self) -> list[Job]:
        return [job for job in self._jobs if job.running]

    def jobs_of(self, owner: object) -> list[Job]:
        return [job for job in self.active_jobs if job.owner is owner]

    def move(self, job: Job, delta: int) -> bool:
        """Move a queued job earlier (negative) or later in the queue."""
        index = next((i for i, e in enumerate(self._waiting) if e[0] is job), None)
        if index is None:
            return False
        target = max(0, min(len(self._waiting) - 1, index + delta))
        if target == index:
            return False
        self._waiting.insert(target, self._waiting.pop(index))
        self.changed.emit()
        return True

    def remove(self, job: Job) -> None:
        """Forget an ended job (the Jobs window's list)."""
        if not job.running and job in self._jobs:
            self._jobs.remove(job)
            self.changed.emit()

    def clear_finished(self) -> None:
        self._jobs = [job for job in self._jobs if job.running]
        self.changed.emit()

    def cancel_all(self, owner: object = None) -> None:
        """Cancel every job (or every job of *owner*); queued ones end at once."""
        for job in list(self.active_jobs):
            if owner is None or job.owner is owner:
                job.cancel()

    def shutdown(self, timeout_ms: int = 5000) -> None:
        self.cancel_all()
        self._threads.waitForDone(timeout_ms)
        self._heavy.waitForDone(timeout_ms)
        if self._awake is not None:
            self._awake.hold(False)
        if self._pool is not None:
            self._pool.shutdown(wait=False, cancel_futures=True)
            self._pool = None
