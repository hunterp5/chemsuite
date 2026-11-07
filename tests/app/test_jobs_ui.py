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
"""Background jobs in the window: the status-bar indicator, the Jobs window, closing."""

from __future__ import annotations

import threading

from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QMessageBox

from chemsuite.app.jobs import FAILED, QUEUE, RUNNING
from chemsuite.app.tools import ToolContext
from chemsuite.core.jobs import check

from .conftest import gui

pytestmark = gui


def _gate():
    release = threading.Event()

    def work(job):
        while not release.wait(0.01):
            check(job.cancel_token)
        return "result"

    return work, release


def test_tool_job_runs_in_background_and_records_its_result(window, qtbot) -> None:
    ctx = ToolContext(window)
    work, release = _gate()
    applied: list[str] = []

    def done(result):
        applied.append(result)
        ctx.show_message("Added 2 columns")

    job = ctx.run("Calculating…", work, done)
    assert job.lane == QUEUE and job.document == "Untitled"
    qtbot.waitUntil(lambda: job.state == RUNNING, timeout=5000)
    assert not window.jobs_indicator.isHidden()
    assert window.jobs_indicator.button.text() == "Calculating…"
    job.report(120, 500)
    window.jobs_indicator.refresh()
    assert window.jobs_indicator.button.text() == "Calculating… (120/500)"
    # Left-aligned with the other status text; no progress bar.
    status = window.statusBar()
    row = status.document.parentWidget().layout()
    order = [row.itemAt(i).widget() for i in range(row.count())]
    assert order == [window.jobs_indicator, status.message, status.document]
    window.show_message("Hello")
    assert status.currentMessage() == "Hello" and not window.jobs_indicator.isHidden()
    assert window.isEnabled()  # nothing modal: the window keeps working
    release.set()
    qtbot.waitUntil(lambda: not job.running, timeout=5000)
    assert applied == ["result"]
    assert job.summary == "Added 2 columns"
    assert window.jobs_indicator.isHidden()


def test_result_dropped_when_document_replaced(window, qtbot) -> None:
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    ctx = ToolContext(window)
    work, release = _gate()
    applied: list[str] = []
    job = ctx.run("Calculating…", work, applied.append)
    window._install(Workspace(Table()))
    release.set()
    qtbot.waitUntil(lambda: not job.running, timeout=5000)
    assert applied == []
    assert "Not applied" in job.summary


def test_failed_job_reports_without_blocking(window, qtbot) -> None:
    ctx = ToolContext(window)

    def broken(job):
        raise RuntimeError("bad input")

    job = ctx.run("Clustering…", broken, lambda r: None)
    qtbot.waitUntil(lambda: job.state == FAILED, timeout=5000)
    box = window._failure_box
    assert box.isVisible() and "bad input" in box.informativeText()
    assert "RuntimeError" in box.detailedText()
    assert "1 failed" in window.jobs_indicator.button.text()
    jobs = window.show_jobs()
    assert window.jobs_indicator.isHidden()  # seen
    jobs.select(job)
    text = jobs.copy_details(job)
    assert "Clustering" in text and "RuntimeError: bad input" in text
    assert QGuiApplication.clipboard().text() == text
    box.close()
    jobs.close()


def test_jobs_window_controls(window, qtbot) -> None:
    ctx = ToolContext(window)
    work, release = _gate()
    first = ctx.run("First…", work, lambda r: None)
    second = ctx.run("Second…", lambda job: None, lambda r: None)
    third = ctx.run("Third…", lambda job: None, lambda r: None)
    qtbot.waitUntil(lambda: first.state == RUNNING, timeout=5000)
    jobs = window.show_jobs()
    titles = [jobs.table.item(r, 0).text() for r in range(jobs.table.rowCount())]
    assert titles[:3] == ["First…", "Second…", "Third…"]
    jobs.select(third)
    assert jobs.up_button.isEnabled() and not jobs.down_button.isEnabled()
    jobs.up_button.click()
    assert window._runner.queued_jobs == [third, second]
    jobs.select(second)
    jobs.cancel_button.click()
    assert second.state == "cancelled"
    release.set()
    qtbot.waitUntil(lambda: not third.running, timeout=5000)
    jobs.refresh()
    assert {jobs.table.item(r, 2).text() for r in range(jobs.table.rowCount())} == {
        "Done",
        "Cancelled",
    }
    window._runner.clear_finished()
    assert jobs.table.rowCount() == 0
    jobs.close()


def test_close_asks_about_running_jobs(window, qtbot, monkeypatch) -> None:
    ctx = ToolContext(window)
    work, release = _gate()
    job = ctx.run("Long…", work, lambda r: None)
    answers = [QMessageBox.StandardButton.Cancel]
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: answers.pop(0))
    assert not window.close()
    assert job.running
    answers.append(QMessageBox.StandardButton.Yes)
    window.close()
    qtbot.waitUntil(lambda: not job.running, timeout=5000)
    assert job.state == "cancelled"
    release.set()


def test_cancel_shows_cancelling_then_keeps_partial_results(window, qtbot) -> None:
    from chemsuite.core.jobs import Cancelled

    ctx = ToolContext(window)
    started = threading.Event()

    def work(job):
        started.set()
        rows = []
        for i in range(1000):
            if job.cancel_token.cancelled:
                raise Cancelled(partial=rows)
            rows.append(i)
            threading.Event().wait(0.01)
        return rows

    applied: list[list[int]] = []

    def done(rows):
        applied.append(rows)
        ctx.show_message(f"Wrote {len(rows)} rows")

    job = ctx.run("Generating conformers…", work, done, keep_partial=True)
    qtbot.waitUntil(lambda: job.state == RUNNING and started.is_set(), timeout=5000)
    job.cancel()
    assert window.jobs_indicator.button.text().startswith("Cancelling generating conformers")
    qtbot.waitUntil(lambda: not job.running, timeout=5000)
    assert len(applied) == 1 and 0 < len(applied[0]) < 1000
    assert job.summary.startswith("Cancelled, partial results kept: Wrote")
    assert "partial results kept" in window.statusBar().currentMessage()


def test_cancel_without_keep_partial_discards(window, qtbot) -> None:
    from chemsuite.core.jobs import Cancelled

    ctx = ToolContext(window)

    def work(job):
        while not job.cancel_token.cancelled:
            threading.Event().wait(0.01)
        raise Cancelled(partial=[1, 2])

    applied: list = []
    job = ctx.run("Clustering…", work, applied.append)
    qtbot.waitUntil(lambda: job.state == RUNNING, timeout=5000)
    job.cancel()
    qtbot.waitUntil(lambda: not job.running, timeout=5000)
    assert applied == [] and window.statusBar().currentMessage() == "Cancelled: Clustering"


def test_jobs_window_removes_ended_jobs_one_by_one(window, qtbot) -> None:
    ctx = ToolContext(window)
    first = ctx.run("First…", lambda job: None, lambda r: None)
    second = ctx.run("Second…", lambda job: None, lambda r: None)
    qtbot.waitUntil(lambda: not first.running and not second.running, timeout=5000)
    work, release = _gate()
    running = ctx.run("Running…", work, lambda r: None)
    try:
        qtbot.waitUntil(lambda: running.state == RUNNING, timeout=5000)
        jobs = window.show_jobs()
        jobs.select(first)
        assert jobs.remove_button.isEnabled()
        jobs.remove_button.click()
        assert first not in window._runner.jobs and second in window._runner.jobs
        jobs.select(running)
        assert not jobs.remove_button.isEnabled()  # a running job is cancelled, not removed
        assert jobs.remove_selected() == 0 and running in window._runner.jobs
        jobs.table.selectAll()  # several at once: only the ended ones go
        assert jobs.remove_selected() == 1
        assert window._runner.jobs == [running]
    finally:
        release.set()
    qtbot.waitUntil(lambda: not running.running, timeout=5000)
    jobs.close()


def test_notification_commands_pass_text_as_arguments(monkeypatch) -> None:
    from chemsuite.app import notify

    title, message = 'Conformers "x" (a.csd)', 'Done" & do shell script "rm -rf ~'
    argv = notify.command(title, message, "darwin")
    assert argv[0] == "osascript" and argv[-2:] == [title, message]
    assert all(title not in part and message not in part for part in argv[:-2])
    monkeypatch.setattr(notify.shutil, "which", lambda name: "/usr/bin/notify-send")
    assert notify.command(title, message, "linux") == [
        "/usr/bin/notify-send",
        "--app-name",
        "chemsuite",
        title,
        message,
    ]
    monkeypatch.setattr(notify.shutil, "which", lambda name: None)
    assert notify.command(title, message, "linux") is None  # falls back to the tray
    assert notify.command(title, message, "win32") is None


def test_a_long_job_ending_unseen_notifies_unless_turned_off(window, monkeypatch) -> None:
    from chemsuite.app import notify
    from chemsuite.app.jobs import Job
    from chemsuite.app.settings_hub import settings_hub

    posted = []
    monkeypatch.setattr(notify, "notify", lambda title, text: posted.append((title, text)))
    monkeypatch.setattr(window, "isActiveWindow", lambda: False)
    job = Job("Generating conformers…", document="a.csd")
    job.started_at, job.ended_at = 0.0, 75.0
    window._notify_unseen(job, "500 conformers")
    assert posted == [("Generating conformers (a.csd)", "500 conformers after 1:15")]

    def off(settings) -> None:
        settings.jobs.notify = False

    settings_hub().update(off)
    window._notify_unseen(job, "again")
    assert len(posted) == 1
    monkeypatch.setattr(window, "isActiveWindow", lambda: True)
    settings_hub().update(lambda s: setattr(s.jobs, "notify", True))
    window._notify_unseen(job, "seen")
    assert len(posted) == 1  # the window is in front: nothing to tell
