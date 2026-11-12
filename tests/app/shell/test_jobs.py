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

from chemsuite.app.shell.jobs import JobRunner
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


def _gate():
    """A job function that waits until released, and the event that releases it."""
    release = threading.Event()

    def work(job):
        while not release.wait(0.01):
            check(job.cancel_token)
        return job.title

    return work, release


def test_queue_runs_one_at_a_time_in_order(qtbot, runner: JobRunner) -> None:
    from chemsuite.app.shell.jobs import QUEUE, QUEUED, RUNNING

    order: list[str] = []
    work, release = _gate()
    first = runner.submit("first", work, lane=QUEUE)
    second = runner.submit("second", lambda job: order.append("second"), lane=QUEUE)
    third = runner.submit("third", lambda job: order.append("third"), lane=QUEUE)
    qtbot.waitUntil(lambda: first.state == RUNNING, timeout=5000)
    assert second.state == QUEUED and third.state == QUEUED
    assert runner.queued_jobs == [second, third]
    assert runner.move(third, -1)
    assert runner.queued_jobs == [third, second]
    assert not runner.move(third, -1)  # already first
    release.set()
    qtbot.waitUntil(lambda: not second.running, timeout=5000)
    assert order == ["third", "second"]
    assert runner.active_jobs == []


def test_cancel_queued_job_ends_it_without_running(qtbot, runner: JobRunner) -> None:
    from chemsuite.app.shell.jobs import CANCELLED, QUEUE, RUNNING

    ran: list[bool] = []
    work, release = _gate()
    first = runner.submit("first", work, lane=QUEUE)
    second = runner.submit("second", lambda job: ran.append(True), lane=QUEUE)
    qtbot.waitUntil(lambda: first.state == RUNNING, timeout=5000)
    with qtbot.waitSignal(second.cancelled, timeout=1000):
        second.cancel()
    assert second.state == CANCELLED and runner.queued_jobs == []
    release.set()
    qtbot.waitUntil(lambda: not first.running, timeout=5000)
    assert ran == []


def test_quick_lane_does_not_wait_for_the_queue(qtbot, runner: JobRunner) -> None:
    from chemsuite.app.shell.jobs import FINISHED, QUEUE, QUICK

    work, release = _gate()
    heavy = runner.submit("heavy", work, lane=QUEUE)
    quick = runner.submit("quick", lambda job: 42, lane=QUICK)
    qtbot.waitUntil(lambda: quick.state == FINISHED, timeout=5000)
    assert heavy.running
    release.set()
    qtbot.waitUntil(lambda: not heavy.running, timeout=5000)


def test_failure_keeps_traceback_and_cancel_all_by_owner(qtbot, runner: JobRunner) -> None:
    from chemsuite.app.shell.jobs import FAILED, QUEUE

    def broken(job):
        raise ValueError("no structures")

    job = runner.submit("broken", broken, lane=QUEUE)
    qtbot.waitUntil(lambda: job.state == FAILED, timeout=5000)
    assert "ValueError: no structures" in job.error_text
    assert job.elapsed >= 0 and job.started_at is not None

    a, b = object(), object()
    work, release = _gate()
    mine = runner.submit("mine", work, lane=QUEUE, owner=a)
    other = runner.submit("other", lambda job: None, lane=QUEUE, owner=b)
    assert runner.jobs_of(a) == [mine]
    runner.cancel_all(owner=a)
    qtbot.waitUntil(lambda: not other.running, timeout=5000)
    assert mine.state == "cancelled" and other.state == "finished"
    runner.clear_finished()
    assert runner.jobs == []


def test_keep_awake_follows_the_queue(qtbot, qapp, monkeypatch) -> None:
    from chemsuite.app.shell.jobs import QUEUE
    from chemsuite.core.jobs import awake

    calls: list[bool] = []
    monkeypatch.setattr(awake.KeepAwake, "_start", lambda self: calls.append(True))
    monkeypatch.setattr(awake.KeepAwake, "_stop", lambda self: calls.append(False))
    runner = JobRunner(keep_awake=True)
    try:
        work, release = _gate()
        job = runner.submit("long", work, lane=QUEUE)
        runner.submit("quick", lambda job: None)  # the quick lane does not hold it
        assert calls == [True]
        release.set()
        qtbot.waitUntil(lambda: not job.running, timeout=5000)
        assert calls == [True, False]
    finally:
        runner.shutdown()
