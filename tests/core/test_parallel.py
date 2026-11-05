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
"""ordered_map: results in order, and a cancel that does not wait for the remaining work."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map


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
