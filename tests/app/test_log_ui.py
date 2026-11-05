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
"""The Log window, error reports and the Help menu."""

from __future__ import annotations

import logging

from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QMessageBox

from chemsuite.app.log_ui import LogWindow
from chemsuite.core.logbook import Logbook, LogbookHandler

from .conftest import gui

pytestmark = gui


def _filled() -> Logbook:
    book = Logbook()
    book.add(logging.DEBUG, "chemsuite.io", "reading header")
    book.add(logging.INFO, "chemsuite.app.jobs", "queued: Clustering", job=3)
    book.add(logging.WARNING, "chemsuite.io", "2 rows have no structure")
    book.add(
        logging.ERROR, "chemsuite.app.jobs", "failed: Hansch", details="Traceback…\nValueError"
    )
    return book


def test_filters_by_level_and_text(qtbot) -> None:
    window = LogWindow(_filled())
    qtbot.addWidget(window)
    assert window.model.rowCount() == 3  # Info and above
    window.level.setCurrentIndex(0)
    assert window.model.rowCount() == 4
    window.level.setCurrentIndex(2)
    assert [e.message for e in window.model.rows] == [
        "2 rows have no structure",
        "failed: Hansch",
    ]
    window.level.setCurrentIndex(0)
    window.search.setText("valueerror")  # details are searched too
    assert [e.message for e in window.model.rows] == ["failed: Hansch"]
    assert window.model.data(window.model.index(0, 3)) == "failed: Hansch …"


def test_new_entries_arrive_and_clear(qtbot) -> None:
    book = _filled()
    window = LogWindow(book)
    qtbot.addWidget(window)
    window.show()
    book.add(logging.INFO, "x", "late arrival")
    qtbot.waitUntil(lambda: window.model.rows[-1].message == "late arrival", timeout=3000)
    assert "4 of 5 entries" in window.count.text()
    book.clear()
    window.poll()
    assert window.model.rowCount() == 0


def test_copy_details_and_job_link(qtbot) -> None:
    opened: list[int] = []
    window = LogWindow(_filled(), show_job=opened.append)
    qtbot.addWidget(window)
    window.view.selectRow(2)  # the error
    assert "ValueError" in window.details.toPlainText()
    text = window.copy()
    assert "failed: Hansch" in text and "ValueError" in text
    assert QGuiApplication.clipboard().text() == text
    window.view.doubleClicked.emit(window.model.index(0, 0))  # queued: Clustering, job 3
    assert opened == [3]


def test_save_as(qtbot, tmp_path, monkeypatch) -> None:
    from PySide6.QtWidgets import QFileDialog

    window = LogWindow(_filled())
    qtbot.addWidget(window)
    target = tmp_path / "log.txt"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), ""))
    window.save_as()
    assert "2 rows have no structure" in target.read_text()


def test_window_log_menu_and_open_failure(window, qtbot, tmp_path) -> None:
    from chemsuite.core.logbook import logbook

    handler = LogbookHandler(logbook())
    logging.getLogger().addHandler(handler)
    try:
        log_window = window.show_log()
        assert window.show_log() is log_window  # one per application
        assert window.command_action("help.log_folder").text() == "Show &Log Folder"
        window.open_path(tmp_path / "missing.csv")
        qtbot.waitUntil(lambda: getattr(window, "_failure_box", None) is not None, timeout=5000)
        box = window._failure_box
        assert isinstance(box, QMessageBox) and "Opening missing.csv failed" in box.text()
        assert "Traceback" in box.detailedText()
        texts = [b.text() for b in box.buttons()]
        assert "Copy Details" in texts and "Show Log Folder" in texts
        box.close()
        log_window.close()
    finally:
        logging.getLogger().removeHandler(handler)


def test_crash_handler_shows_report_from_a_thread(qtbot, monkeypatch) -> None:
    import sys
    import threading

    from chemsuite.app import error_report

    shown: list[tuple[str, str]] = []
    monkeypatch.setattr(
        error_report, "show_error", lambda parent, summary, reason, details, **k: shown.append(
            (reason, details))
    )  # fmt: skip
    hooks = sys.excepthook, threading.excepthook
    try:
        error_report.install_crash_handler()

        def boom() -> None:
            raise RuntimeError("from a worker")

        thread = threading.Thread(target=boom)
        thread.start()
        thread.join()
        qtbot.waitUntil(lambda: bool(shown), timeout=3000)
        assert "RuntimeError: from a worker" in shown[0][0]
        assert "Traceback" in shown[0][1]
    finally:
        sys.excepthook, threading.excepthook = hooks
