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

from PySide6.QtCore import QObject, QRunnable, QThread, QThreadPool, Signal, Slot

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.awake import KeepAwake
from chemsuite.core.jobs.parallel import default_workers, process_pool

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
HISTORY_LIMIT = 200  # ended jobs kept for the Log window's Jobs tab

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
        #: work nobody asked for (made ahead for a view): closing the document cancels it
        #: without asking
        self.background = False
        self.state = QUEUED
        self.queued_at = time.monotonic()
        self.started_at: float | None = None
        self.ended_at: float | None = None
        self.summary = ""  # one line about the result ("Added 3 columns")
        self.error_text = ""  # the traceback of a failure
        #: what a cancelled job finished before stopping (see ``Cancelled.partial``), or None
        self.partial_result: Any = None
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
            log.info("started: %s", self.name, extra={"job": self.id})
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
            log.info("finished in %.1f s: %s", self.elapsed, self.name, extra={"job": self.id})
            self.finished.emit(value)
        elif outcome == FAILED:
            log.error(  # the traceback goes with the record (the Log window's details)
                "failed: %s: %s",
                self.name,
                value,
                exc_info=value if isinstance(value, BaseException) else None,
                extra={"job": self.id},
            )
            self.failed.emit(value)
        else:
            self.partial_result = value
            log.info(
                "cancelled%s: %s",
                " with partial results" if value is not None else "",
                self.name,
                extra={"job": self.id},
            )
            self.cancelled.emit()
        self.done.emit()
        self.changed.emit()

    # called from the GUI thread -------------------------------------------------------------
    @property
    def name(self) -> str:
        """The title without a trailing ellipsis, for messages."""
        return self.title.rstrip("…").rstrip()

    @property
    def running(self) -> bool:
        """Not ended yet: queued or running."""
        return self.state not in ENDED

    @property
    def cancelling(self) -> bool:
        """Cancel was asked for but the work has not stopped yet."""
        return self.state == RUNNING and self.cancel_token.cancelled

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
        else:
            self.changed.emit()  # show "Cancelling…" at once


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
            except Cancelled as stopped:
                value = stopped.partial  # what was finished before the cancel, if anything
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
        self._wide: ProcessPoolExecutor | None = None  # one worker, many threads (wide_process)
        self._pool_lock = threading.Lock()  # jobs on two threads may both want the pool first
        #: Settings ▸ General ▸ Processor cores, read on the GUI thread (see process_pool)
        self._pool_workers: int | None = None
        self._jobs: list[Job] = []  # every job not cleared, in the order submitted
        self._waiting: list[tuple[Job, Callable[[Job], Any]]] = []
        self._current: Job | None = None
        self._awake = KeepAwake() if keep_awake else None

    def process_pool(self) -> ProcessPoolExecutor:
        """The shared worker processes, started on first use (usually from a job's thread).

        The worker count comes from Settings, read on the GUI thread when a job is submitted:
        reading it here could create the settings hub on a worker thread, and the hub applies
        the theme to Qt widgets, which must only happen on the GUI thread.
        """
        with self._pool_lock:
            if self._pool is None:
                workers = self._pool_workers
                if workers is None and QThread.currentThread() is self.thread():
                    workers = self._configured_workers()
                self._pool = process_pool(workers or None)  # 0 or unknown: automatic
            return self._pool

    def wide_process(self) -> ProcessPoolExecutor:
        """One worker process whose numeric libraries (numba, BLAS) use as many threads as
        the shared pool has workers: for one multi-threaded task that holds Python's GIL
        (UMAP, t-SNE), which would freeze the window from a thread. Started on first use and
        kept, so what it compiles is compiled once."""
        with self._pool_lock:
            if self._wide is None:
                workers = self._pool_workers
                if workers is None and QThread.currentThread() is self.thread():
                    workers = self._configured_workers()
                self._wide = process_pool(1, threads=workers or default_workers())
            return self._wide

    @staticmethod
    def _configured_workers() -> int:
        from .settings_hub import settings_hub

        return settings_hub().settings.performance.max_workers

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
        if self._pool is None and self._pool_workers is None:
            self._pool_workers = self._configured_workers()  # on the GUI thread, see process_pool
        job = Job(title, self, lane=lane, owner=owner, document=document)
        self._jobs.append(job)
        job.changed.connect(self.changed)
        job.done.connect(lambda: self._ended(job))
        log.info("queued (%s lane): %s", lane, job.name, extra={"job": job.id})
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
        """Forget an ended job (the Log window's Jobs tab)."""
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
        for pool in (self._pool, self._wide):
            if pool is not None:
                pool.shutdown(wait=False, cancel_futures=True)
        self._pool = self._wide = None
