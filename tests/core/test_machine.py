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

"""Worker processes sized to the machine: its physical cores and its memory."""

from __future__ import annotations

import pytest

from chemsuite.core.jobs import machine, parallel
from chemsuite.core.jobs.machine import Memory


def test_the_machine_is_read() -> None:
    assert 1 <= machine.physical_cores() <= machine.logical_cores()
    found = machine.memory()
    assert found is not None and found.total_mb > 0
    assert found.available_mb is None or 0 <= found.available_mb <= found.total_mb


@pytest.mark.parametrize(
    ("cores", "memory", "workers"),
    [
        (8, Memory(16_384, 10_000), 7),  # 8 cores, 16 hyperthreads, 16 GB: one core kept
        (16, Memory(16_384, 10_000), 15),  # 16 physical cores: 4 GB holds 16 workers
        (16, Memory(16_384, 3_000), 6),  # little free memory: half of it
        (8, Memory(4_096, 2_000), 4),  # 4 GB: a quarter is 1 GB
        (8, Memory(2_048, 200), 1),  # never none
        (1, Memory(16_384, 10_000), 1),
        (8, None, 7),  # memory unknown: the cores decide
        (8, Memory(16_384, None), 7),
    ],
)
def test_workers_leave_a_core_and_most_memory_free(monkeypatch, cores, memory, workers) -> None:
    monkeypatch.setattr(machine, "physical_cores", lambda: cores)
    monkeypatch.setattr(machine, "memory", lambda: memory)
    assert parallel.default_workers() == workers
