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

from __future__ import annotations

import threading

import pytest

from chemsuite.app.jobs import JobRunner
from chemsuite.core.jobs import check


@pytest.fixture
def runner(qapp):
    r = JobRunner(threads=2)
    yield r
    r.shutdown()


def test_signals_arrive_on_the_gui_thread(qtbot, runner: JobRunner) -> None:
    gui = threading.get_ident()
    seen: list[tuple[str, object, bool]] = []

    def work(job):
        for i in range(3):
            job.emit_partial(i)
            job.report(i + 1, 3)
        return "ok"

    job = runner.submit("t", work)
    job.partial.connect(lambda v: seen.append(("partial", v, threading.get_ident() == gui)))
    job.finished.connect(lambda v: seen.append(("finished", v, threading.get_ident() == gui)))
    with qtbot.waitSignal(job.done, timeout=5000):
        pass
    assert seen == [
        ("partial", 0, True),
        ("partial", 1, True),
        ("partial", 2, True),
        ("finished", "ok", True),
    ]
    assert job.progress == (3, 3) and not job.running
    assert runner.active_jobs == []


def test_failure_and_cancel(qtbot, runner: JobRunner) -> None:
    errors: list[BaseException] = []
    job = runner.submit("bad", lambda job: 1 / 0)
    job.failed.connect(errors.append)
    with qtbot.waitSignal(job.done, timeout=5000):
        pass
    assert isinstance(errors[0], ZeroDivisionError)

    started = threading.Event()

    def spin(job):
        started.set()
        while True:
            check(job.cancel_token)

    job = runner.submit("spin", spin)
    started.wait(5)
    with qtbot.waitSignal(job.cancelled, timeout=5000):
        job.cancel()
