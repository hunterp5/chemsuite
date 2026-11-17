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

"""Quitting: what closing a document window asks about running work, and what closes with
the last one.

- Closing a document window while others stay open asks only about that document's jobs.
- Closing the last one asks about every job the user started, from any window, and
  names the background tasks that stop with them. Background work alone stops without asking.
- Once the last document window has closed, every other window (the Log, a Sketcher, tool
  windows) closes too and the app quits; ``main.run_app`` then stops whatever still runs.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QCoreApplication, QObject, QTimer
from PySide6.QtWidgets import QApplication, QMessageBox, QWidget

from .jobs import JobRunner


class _State:
    #: Set by ``main.run_app`` only: closing the last document window quits the app. Tests
    #: open and close windows freely without quitting the test application.
    with_last_window = False


state = _State()


def _count(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


def jobs_question(runner: JobRunner, owner: QObject, *, last: bool) -> str | None:
    """What to ask before closing *owner*'s window, or None to close without asking.

    *last*: no other document window is open, so closing it quits the app.
    """
    if not (last and state.with_last_window):
        jobs = [job for job in runner.jobs_of(owner) if not job.background]
        if not jobs:
            return None
        return (
            f"{len(jobs)} job(s) for this document are still running or queued.\n\n"
            "Cancel them and close?"
        )
    return quit_question(runner)


def quit_question(runner: JobRunner) -> str | None:
    """What to ask before quitting, or None when only background work (or nothing) runs."""
    active = runner.active_jobs
    started = [job for job in active if not job.background]
    if not started:
        return None
    one = len(started) == 1
    background = len(active) - len(started)
    also = f" (and {_count(background, 'background task')})" if background else ""
    return (
        f"{_count(len(started), 'job')} {'is' if one else 'are'} still running or queued{also}."
        f"\n\nCancel {'it' if one else 'them'} and quit?"
    )


def confirm(parent: QWidget, question: str) -> bool:
    answer = QMessageBox.question(
        parent,
        "chemsuite",
        question,
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
    )
    return answer == QMessageBox.StandardButton.Yes


def window_closed(none_left: Callable[[], bool]) -> None:
    """A document window closed; if *none_left* still holds once the close is done, close
    everything else and quit."""
    if state.with_last_window:
        QTimer.singleShot(0, lambda: none_left() and close_everything())


def close_everything() -> None:
    """Close every remaining window and leave the event loop."""
    for widget in QApplication.topLevelWidgets():
        if widget.isVisible():
            widget.close()
    QCoreApplication.exit(0)
