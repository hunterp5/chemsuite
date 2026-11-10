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
"""ordered_map: results in order, and a cancel that does not wait for the remaining work.
Pool workers: one numeric thread, below-normal priority, gone when the app is."""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import THREAD_VARS, ordered_map, process_pool


def slow_square(x: int) -> int:
    time.sleep(0.2 if x else 0)
    return x * x


@pytest.mark.parametrize("pooled", [False, True])
def test_results_in_order(pooled) -> None:
    executor = ThreadPoolExecutor(3) if pooled else None
    try:
        assert list(ordered_map(lambda x: x * x, range(10), executor)) == [x * x for x in range(10)]
    finally:
        if executor:
            executor.shutdown()


@pytest.mark.parametrize("pooled", [False, True])
def test_cancel_stops_without_waiting_for_running_work(pooled) -> None:
    executor = ThreadPoolExecutor(2) if pooled else None
    token = CancelToken()
    seen: list[int] = []
    start = time.monotonic()
    try:
        with pytest.raises(Cancelled):
            for value in ordered_map(slow_square, range(40), executor, cancel=token):
                seen.append(value)
                token.cancel()  # after the first (instant) result
        elapsed = time.monotonic() - start
    finally:
        if executor:
            executor.shutdown(cancel_futures=True)
    assert seen == [0]
    assert elapsed < 0.15  # the cancel did not wait for a 0.2 s item to finish


def worker_state() -> dict:
    import numpy  # noqa: F401 - loaded after start-up, as a real task would
    from threadpoolctl import threadpool_info

    return {
        "vars": {name: os.environ.get(name) for name in THREAD_VARS},
        "threads": {info["num_threads"] for info in threadpool_info()},
        "nice": None if sys.platform == "win32" else os.getpriority(os.PRIO_PROCESS, 0),
    }


def test_pool_workers_use_one_thread_at_lower_priority() -> None:
    with process_pool(1) as pool:
        state = pool.submit(worker_state).result(timeout=60)
    assert set(state["vars"].values()) == {"1"}
    assert state["threads"] <= {1}
    if state["nice"] is not None:
        assert state["nice"] > os.getpriority(os.PRIO_PROCESS, 0)


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)  # signal 0: existence check only
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


@pytest.mark.skipif(sys.platform == "win32", reason="signal-0 liveness check is POSIX")
def test_pool_workers_exit_when_the_app_dies(tmp_path) -> None:
    """A crashed or killed app (no shutdown) does not leave its workers behind."""
    script = (
        "import os, sys\n"
        "from chemsuite.core.jobs.parallel import process_pool\n"
        "pool = process_pool(1)\n"
        "print(pool.submit(os.getpid).result(timeout=60), flush=True)\n"
        "os._exit(1)\n"
    )
    out = tmp_path / "worker.pid"
    with out.open("w") as sink:  # a file, not a pipe: a lingering worker would hold a pipe open
        subprocess.run(
            [sys.executable, "-c", script], stdout=sink, stdin=subprocess.DEVNULL, timeout=120
        )
    worker = int(out.read_text().strip())
    deadline = time.monotonic() + 10
    while _alive(worker) and time.monotonic() < deadline:
        time.sleep(0.05)
    lingering = _alive(worker)
    if lingering:
        os.kill(worker, signal.SIGKILL)  # do not leave it behind when the test fails
    assert not lingering
