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

"""Process-pool helpers that bound memory and keep results in order."""

from __future__ import annotations

import multiprocessing
import multiprocessing.connection
import os
import sys
import threading
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Executor, Future, ProcessPoolExecutor, wait
from typing import TypeVar

from . import CancelToken, check

T = TypeVar("T")
R = TypeVar("R")
#: how often a wait for a worker's result looks at the cancel token (seconds)
CANCEL_POLL_S = 0.05


def default_workers() -> int:
    """Leave one core for the GUI; at least one worker."""
    return max(1, (os.cpu_count() or 2) - 1)


#: environment variables that size the thread pools of OpenMP, BLAS, numexpr and numba
THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMBA_NUM_THREADS",
)


def limit_app_threads(threads: int | None = None) -> None:
    """Size the app process's own numeric thread pools to *threads* (default: one per core but
    one), unless the user set them. Call before numpy, scikit-learn or numba is imported."""
    count = str(threads or default_workers())
    for name in THREAD_VARS:
        os.environ.setdefault(name, count)


def process_pool(workers: int | None = None) -> ProcessPoolExecutor:
    """A pool using the *spawn* start method on every OS.

    Spawn behaves identically on Windows, macOS and Linux and never forks a process that holds
    Qt or open SQLite handles. Workers import only what the submitted function needs, and start
    with :func:`worker_start`.
    """
    return ProcessPoolExecutor(
        max_workers=workers or default_workers(),
        mp_context=multiprocessing.get_context("spawn"),
        initializer=worker_start,
    )


def worker_start() -> None:
    """Set up a pool worker before its first task.

    * One numeric thread each: the pool already runs a worker per core, and a BLAS, OpenMP or
      numba pool per worker on top would put cores x cores busy threads on the machine.
    * Below-normal priority, so the app window and the rest of the system stay responsive
      while the workers keep every core busy.
    * Exit when the app does, even if it is killed or crashes, instead of lingering.
    """
    for name in THREAD_VARS:
        os.environ[name] = "1"  # libraries loaded from now on
    try:  # libraries already loaded
        from threadpoolctl import threadpool_limits

        threadpool_limits(1)
    except Exception:  # noqa: BLE001, S110 - optional; never stop a worker from starting
        pass
    numba = sys.modules.get("numba")
    if numba is not None:
        numba.set_num_threads(1)
    lower_priority()
    exit_with_parent()


#: Windows BELOW_NORMAL_PRIORITY_CLASS
_BELOW_NORMAL = 0x4000
#: Unix nice increment for background work
_NICE = 10


def lower_priority() -> None:
    """Run this process below normal priority (best effort)."""
    try:
        if sys.platform == "win32":
            import ctypes

            kernel = ctypes.windll.kernel32  # type: ignore[attr-defined]
            kernel.SetPriorityClass(kernel.GetCurrentProcess(), _BELOW_NORMAL)
        else:
            os.nice(_NICE)
    except Exception:  # noqa: BLE001, S110 - priority is a courtesy, not a requirement
        pass


def exit_with_parent() -> None:
    """End this (spawned) process as soon as its parent process ends."""
    parent = multiprocessing.parent_process()
    if parent is None:
        return

    def watch() -> None:
        multiprocessing.connection.wait([parent.sentinel])
        os._exit(0)

    threading.Thread(target=watch, name="exit-with-parent", daemon=True).start()


def ordered_map(
    fn: Callable[[T], R],
    items: Iterable[T],
    executor: Executor | None = None,
    *,
    window: int | None = None,
    cancel: CancelToken | None = None,
) -> Iterator[R]:
    """Like ``map(fn, items)``, run on *executor*, yielding results in input order.

    At most *window* items are in flight, so a large input (a 1 GB SDF split into blocks) is
    never all in memory at once. Without an executor, runs inline.

    With *cancel*, raises :class:`~chemsuite.core.jobs.Cancelled` soon after it is set: before
    each inline item, and every few milliseconds while waiting for a worker, rather than only
    when a slow chunk comes back. Work not started is cancelled; a worker already busy with a
    chunk finishes that chunk on its own (a process cannot be interrupted safely).
    """
    if executor is None:
        for item in items:
            check(cancel)
            yield fn(item)
        return
    window = window or 2 * (getattr(executor, "_max_workers", None) or default_workers())
    pending: deque[Future[R]] = deque()

    def next_result() -> R:
        head = pending.popleft()
        if cancel is not None:
            while not head.done():
                check(cancel)
                wait([head], timeout=CANCEL_POLL_S)
        return head.result()

    try:
        for item in items:
            check(cancel)
            pending.append(executor.submit(fn, item))
            if len(pending) >= window:
                yield next_result()
        while pending:
            yield next_result()
    finally:
        for future in pending:
            future.cancel()
