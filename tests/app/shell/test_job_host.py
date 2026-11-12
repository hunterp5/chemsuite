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

"""A window's jobs: their outcomes in its status bar, and none once it is gone."""

from __future__ import annotations

import threading

import pytest
from PySide6.QtWidgets import QWidget

from chemsuite.app.shell.job_host import JobHost
from chemsuite.app.shell.jobs import JobRunner
from chemsuite.core.jobs import check
from tests.gui import destroy_now


@pytest.fixture
def runner(qapp):
    r = JobRunner(threads=2)
    yield r
    r.shutdown()


def _spin(started: threading.Event):
    def work(job):
        started.set()
        while True:
            check(job.cancel_token)

    return work


def _host(runner: JobRunner):
    window = QWidget()
    shown: list[str] = []

    def message(text: str, ms: int) -> None:
        window.setWindowTitle(text)  # raises once the window is deleted, as a status bar does
        shown.append(text)

    host = JobHost(window, runner, current=lambda: None, name=lambda: "", message=message)
    return window, host, shown


def test_a_cancelled_job_says_so(qtbot, runner: JobRunner) -> None:
    window, host, shown = _host(runner)
    started = threading.Event()
    job = host.run("Spin…", _spin(started), lambda result: None)
    started.wait(5)
    with qtbot.waitSignal(job.done, timeout=5000):
        job.cancel()
    assert shown == ["Cancelled: Spin"]
    destroy_now(window)


def test_a_job_that_ends_after_its_window_closed_reports_nothing(qtbot, runner) -> None:
    window, host, shown = _host(runner)
    started = threading.Event()
    applied: list[object] = []
    cancelled = host.run("Spin…", _spin(started), applied.append)
    started.wait(5)
    finished = host.run("Quick…", lambda job: 1, applied.append)
    failed = host.run("Bad…", lambda job: 1 / 0, applied.append)
    destroy_now(window)  # before any outcome is delivered
    cancelled.cancel()
    for job in (cancelled, finished, failed):
        if job.running:
            qtbot.waitSignal(job.done, timeout=5000).wait()
    qtbot.wait(50)
    assert shown == [] and applied == []
