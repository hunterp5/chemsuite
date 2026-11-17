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

"""What closing a window asks about running work, and what quitting asks."""

from __future__ import annotations

import threading

import pytest
from PySide6.QtCore import QObject

from chemsuite.app.shell import quitting
from chemsuite.app.shell.jobs import JobRunner


@pytest.fixture
def runner(qapp, monkeypatch):
    monkeypatch.setattr(quitting.state, "with_last_window", True)
    release = threading.Event()
    r = JobRunner(threads=8)
    r.release = release
    yield r
    release.set()
    r.shutdown()


def submit(runner, owner=None, *, background=False):
    job = runner.submit("work", lambda job: runner.release.wait(30), owner=owner)
    job.background = background
    return job


def test_another_window_still_open_asks_only_about_this_document(runner) -> None:
    mine, other = QObject(), QObject()
    submit(runner, mine)
    submit(runner, other)
    submit(runner, mine, background=True)
    question = quitting.jobs_question(runner, mine, last=False)
    assert question.startswith("1 job(s) for this document")
    assert quitting.jobs_question(runner, QObject(), last=False) is None


def test_the_last_window_asks_about_every_job_the_user_started(runner) -> None:
    window = QObject()
    submit(runner, QObject())
    submit(runner)  # no owner
    submit(runner, background=True)
    question = quitting.jobs_question(runner, window, last=True)
    assert question.startswith("2 jobs are still running or queued (and 1 background task).")
    assert question.endswith("Cancel them and quit?")


def test_one_job_reads_naturally(runner) -> None:
    submit(runner)
    assert quitting.quit_question(runner) == (
        "1 job is still running or queued.\n\nCancel it and quit?"
    )


def test_background_work_alone_quits_without_asking(runner) -> None:
    submit(runner, background=True)
    submit(runner, QObject(), background=True)
    assert quitting.jobs_question(runner, QObject(), last=True) is None


def test_without_quit_on_last_window_closing_the_last_asks_as_before(runner, monkeypatch) -> None:
    monkeypatch.setattr(quitting.state, "with_last_window", False)
    window = QObject()
    submit(runner)
    assert quitting.jobs_question(runner, window, last=True) is None  # not this document's


def test_closing_the_last_window_closes_the_others(qtbot, monkeypatch) -> None:
    from PySide6.QtWidgets import QWidget

    monkeypatch.setattr(quitting.state, "with_last_window", True)
    exits: list[int] = []
    monkeypatch.setattr(quitting.QCoreApplication, "exit", lambda code=0: exits.append(code))
    log = QWidget()
    qtbot.addWidget(log)
    log.show()
    quitting.window_closed(lambda: False)  # another document window opened meanwhile
    qtbot.wait(50)
    assert log.isVisible() and exits == []
    quitting.window_closed(lambda: True)
    qtbot.waitUntil(lambda: exits == [0], timeout=2000)
    assert not log.isVisible()
