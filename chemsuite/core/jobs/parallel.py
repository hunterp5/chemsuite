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
import os
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Executor, Future, ProcessPoolExecutor
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def default_workers() -> int:
    """Leave one core for the GUI; at least one worker."""
    return max(1, (os.cpu_count() or 2) - 1)


def process_pool(workers: int | None = None) -> ProcessPoolExecutor:
    """A pool using the *spawn* start method on every OS.

    Spawn behaves identically on Windows, macOS and Linux and never forks a process that holds
    Qt or open SQLite handles. Workers import only what the submitted function needs.
    """
    return ProcessPoolExecutor(
        max_workers=workers or default_workers(),
        mp_context=multiprocessing.get_context("spawn"),
    )


def ordered_map(
    fn: Callable[[T], R],
    items: Iterable[T],
    executor: Executor | None = None,
    *,
    window: int | None = None,
) -> Iterator[R]:
    """Like ``map(fn, items)``, run on *executor*, yielding results in input order.

    At most *window* items are in flight, so a large input (a 1 GB SDF split into blocks) is
    never all in memory at once. Without an executor, runs inline. If the consumer stops early
    (for example on cancel), queued work is cancelled.
    """
    if executor is None:
        for item in items:
            yield fn(item)
        return
    window = window or 2 * (getattr(executor, "_max_workers", None) or default_workers())
    pending: deque[Future[R]] = deque()
    try:
        for item in items:
            pending.append(executor.submit(fn, item))
            if len(pending) >= window:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()
    finally:
        for future in pending:
            future.cancel()
