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
import time
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Executor, Future, ProcessPoolExecutor, wait
from typing import TypeVar

from . import CancelToken, check, machine

T = TypeVar("T")
R = TypeVar("R")
#: how often a wait for a worker's result looks at the cancel token (seconds)
CANCEL_POLL_S = 0.05


#: a worker process's own memory at its busiest (MB): pKa's measured 70–120 MB with its
#: network shared (ionisation/unipka/mapped.py), with room for larger molecules
WORKER_MB = 250
#: the most of the machine's memory the workers may take together, and of what is free when
#: they start
MEMORY_SHARE, FREE_SHARE = 0.25, 0.5


def default_workers() -> int:
    """Worker processes that leave the machine usable: one per physical core but one (the
    window and other programs keep a core, and every core's second hardware thread), and no
    more than :data:`MEMORY_SHARE` of the memory, or :data:`FREE_SHARE` of what is free now,
    holds at :data:`WORKER_MB` each. At least one."""
    count = machine.physical_cores() - 1
    found = machine.memory()
    if found is not None:
        budget = found.total_mb * MEMORY_SHARE
        if found.available_mb is not None:
            budget = min(budget, found.available_mb * FREE_SHARE)
        count = min(count, int(budget // WORKER_MB))
    return max(1, count)


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


def process_pool(workers: int | None = None, threads: int = 1) -> ProcessPoolExecutor:
    """A pool using the *spawn* start method on every OS.

    Spawn behaves identically on Windows, macOS and Linux and never forks a process that holds
    Qt or open SQLite handles. Workers import only what the submitted function needs, and start
    with :func:`worker_start` (*threads* numeric threads each).
    """
    return ProcessPoolExecutor(
        max_workers=workers or default_workers(),
        mp_context=multiprocessing.get_context("spawn"),
        initializer=worker_start,
        initargs=(threads,),
    )


def worker_start(threads: int = 1) -> None:
    """Set up a pool worker before its first task.

    * *threads* numeric threads: one in the shared pool, which already runs a worker per core
      (a BLAS, OpenMP or numba pool per worker on top would put cores x cores busy threads on
      the machine); more in a pool of one worker that runs a single multi-threaded task.
    * Below-normal priority, so the app window and the rest of the system stay responsive
      while the workers keep every core busy.
    * Exit when the app does, even if it is killed or crashes, instead of lingering.
    """
    for name in THREAD_VARS:
        os.environ[name] = str(threads)  # libraries loaded from now on
    try:  # libraries already loaded
        from threadpoolctl import threadpool_limits

        threadpool_limits(threads)
    except Exception:  # noqa: BLE001, S110 - optional; never stop a worker from starting
        pass
    numba = sys.modules.get("numba")
    if numba is not None:
        numba.set_num_threads(min(threads, numba.config.NUMBA_NUM_THREADS))
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


def stop_pool(pool: ProcessPoolExecutor, grace: float = 2.0) -> bool:
    """Shut *pool* down now, ending workers that are still busy (on the way out).

    ``shutdown(wait=False)`` alone leaves a busy worker running its task, and the interpreter
    then waits for it at exit. Each worker gets *grace* seconds after a terminate request
    before it is killed. Returns whether every worker and the pool's own thread ended.
    """
    workers = list((pool._processes or {}).values())
    manager = pool._executor_manager_thread
    pool.shutdown(wait=False, cancel_futures=True)
    for worker in workers:
        if worker.is_alive():
            worker.terminate()
    deadline = time.monotonic() + grace
    for worker in workers:
        worker.join(max(0.0, deadline - time.monotonic()))
        if worker.is_alive():
            worker.kill()
            worker.join(grace)
    if manager is not None:
        manager.join(grace)
    stuck = any(w.is_alive() for w in workers) or (manager is not None and manager.is_alive())
    return not stuck


def ordered_map(
    fn: Callable[[T], R],
    items: Iterable[T],
    executor: Executor | None = None,
    *,
    window: int | None = None,
    cancel: CancelToken | None = None,
    finished: Callable[[R], None] | None = None,
) -> Iterator[R]:
    """Like ``map(fn, items)``, run on *executor*, yielding results in input order.

    *finished* is called with each result as soon as it is ready, in whatever order the
    workers finish (one call at a time, from a pool thread): for progress that counts every
    item done, not only those whose turn has come.

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
            result = fn(item)
            if finished is not None:
                finished(result)
            yield result
        return
    window = window or 2 * (getattr(executor, "_max_workers", None) or default_workers())
    pending: deque[Future[R]] = deque()
    lock = threading.Lock()

    def report(future: Future[R]) -> None:
        if finished is not None and not future.cancelled() and future.exception() is None:
            with lock:
                finished(future.result())

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
            future = executor.submit(fn, item)
            if finished is not None:
                future.add_done_callback(report)
            pending.append(future)
            if len(pending) >= window:
                yield next_result()
        while pending:
            yield next_result()
    finally:
        for future in pending:
            future.cancel()
