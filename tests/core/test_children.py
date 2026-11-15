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

"""External programs end with everything they started, when stopped or when chemsuite quits."""

from __future__ import annotations

import subprocess
import sys

from chemsuite.core.jobs import children
from tests.processes import alive, gone

# A program that starts a sleeping program of its own, prints that one's pid and sleeps.
PARENT = (
    "import subprocess, sys, time\n"
    "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(600)'])\n"
    "print(child.pid, flush=True)\n"
    "time.sleep(600)\n"
)


def start_parent() -> tuple[subprocess.Popen, int]:
    process = children.start([sys.executable, "-c", PARENT], stdout=subprocess.PIPE, text=True)
    grandchild = int(process.stdout.readline())
    assert alive(process.pid) and alive(grandchild)
    return process, grandchild


def test_stop_ends_the_program_and_what_it_started() -> None:
    process, grandchild = start_parent()
    assert process in children.running()
    children.stop(process, grace=1.0)
    assert process.poll() is not None
    assert gone(grandchild)
    assert process not in children.running()
    process.stdout.close()


def test_stop_all_ends_every_program() -> None:
    started = [start_parent() for _ in range(2)]
    children.stop_all(grace=1.0)
    for process, grandchild in started:
        assert process.poll() is not None
        assert gone(grandchild)
        process.stdout.close()
    assert children.running() == []


def test_stopping_a_finished_program_is_harmless() -> None:
    process = children.start([sys.executable, "-c", "pass"])
    process.wait(timeout=30)
    children.stop(process)
    children.stop(process)
    assert process.returncode == 0
