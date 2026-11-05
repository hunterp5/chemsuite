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
"""The session logbook, log files, the audit log and the uncaught-exception hooks."""

from __future__ import annotations

import json
import logging
import sys
import threading

import pytest

from chemsuite.core import logbook as lb


@pytest.fixture
def restore_logging():
    root = logging.getLogger()
    handlers, level = list(root.handlers), root.level
    hooks = sys.excepthook, threading.excepthook
    yield
    for handler in list(lb._installed):
        root.removeHandler(handler)
        handler.close()
    lb._installed.clear()
    root.handlers[:] = handlers
    root.setLevel(level)
    sys.excepthook, threading.excepthook = hooks
    lb.setup_audit(None)


def test_logbook_is_bounded_and_reads_after_a_sequence() -> None:
    book = lb.Logbook(capacity=3)
    for i in range(5):
        book.add(logging.INFO, "t", f"m{i}")
    assert [e.message for e in book.entries()] == ["m2", "m3", "m4"]
    assert [e.message for e in book.entries(after=4)] == ["m4"]
    assert book.entries(after=5) == []
    assert [e.message for e in book.entries(after=1)] == ["m2", "m3", "m4"]  # 2 was dropped
    book.clear()
    assert len(book) == 0 and book.cleared_through == 5
    assert book.add(logging.INFO, "t", "next").seq == 6


def test_handler_keeps_job_ids_and_tracebacks() -> None:
    book = lb.Logbook()
    logger = logging.getLogger("test.logbook")
    logger.propagate = False
    handler = lb.LogbookHandler(book)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    try:
        logger.info("queued: %s", "Clustering", extra={"job": 7})
        try:
            raise ZeroDivisionError("division by zero")
        except ZeroDivisionError:
            logger.exception("it broke")
    finally:
        logger.removeHandler(handler)
    first, second = book.entries()
    assert first.message == "queued: Clustering" and first.job == 7
    assert second.level == logging.ERROR and "ZeroDivisionError" in second.details
    assert "it broke" in second.as_text() and "Traceback" in second.as_text()


def test_setup_logging_writes_a_file(tmp_path, restore_logging) -> None:
    path = lb.setup_logging(level=logging.DEBUG, folder=tmp_path, to_stderr=False)
    logging.getLogger("chemsuite.test").debug("hello file")
    before = lb.logbook().last_seq
    logging.getLogger("chemsuite.test").warning("hello book")
    for handler in lb._installed:
        handler.flush()
    text = path.read_text()
    assert "hello file" in text and "hello book" in text
    assert lb.logbook().entries(after=before)[-1].message == "hello book"
    # Again (after the policy is read): handlers are replaced, not doubled.
    lb.setup_logging(level=logging.WARNING, folder=tmp_path, to_stderr=False)
    assert sum(isinstance(h, lb.LogbookHandler) for h in logging.getLogger().handlers) == 1


def test_log_folder_and_levels(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv(lb.LOG_ENV, str(tmp_path))
    assert lb.log_folder() == tmp_path
    assert lb.level_from("debug") == logging.DEBUG
    assert lb.level_from(40) == logging.ERROR
    assert lb.level_from("loud") == logging.INFO
    assert lb.level_from(True) == logging.INFO


def test_audit_only_when_enabled(tmp_path, restore_logging) -> None:
    lb.setup_audit(None)
    lb.audit("open", path=tmp_path / "a.csv")
    assert not (tmp_path / lb.AUDIT_FILE).exists()
    target = tmp_path / "audit" / "trail.log"
    assert lb.setup_audit(str(target)) == target
    lb.audit("export", path=tmp_path / "b.csv", rows=3)
    for handler in logging.getLogger(lb.AUDIT_LOGGER).handlers:
        handler.flush()
    record = json.loads(target.read_text().strip())
    assert record["event"] == "export" and record["rows"] == 3
    assert record["path"].endswith("b.csv") and record["user"]
    # Audit records stay out of the ordinary log.
    assert all("export" not in e.message for e in lb.logbook().entries()[-3:])


def test_excepthooks_log_and_report(restore_logging) -> None:
    seen: list[tuple[str, str]] = []
    lb.setup_logging(to_file=False, to_stderr=False)
    lb.install_excepthooks(lambda summary, details: seen.append((summary, details)))
    try:
        raise ValueError("main thread")
    except ValueError:
        sys.excepthook(*sys.exc_info())

    def boom() -> None:
        raise KeyError("worker")

    thread = threading.Thread(target=boom, name="worker-1")
    thread.start()
    thread.join()
    assert seen[0][0] == "ValueError: main thread" and "Traceback" in seen[0][1]
    assert seen[1][0].startswith("KeyError")
    messages = [e.message for e in lb.logbook().entries()[-2:]]
    assert "uncaught exception in thread worker-1" in messages[1]
