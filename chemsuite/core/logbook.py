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
"""The session log: every log record in memory, a rotating log file, and an optional audit log.

- :class:`Logbook` keeps the last 10,000 records (thread-safe) for the Log window. Readers ask
  for entries after a sequence number, so polling is cheap and nothing is signalled per record.
- :func:`setup_logging` wires the root logger to the logbook, to stderr and to a rotating file
  (5 × 2 MB) in :func:`log_folder`. IT can set the level through the policy file
  (``"log_level": "DEBUG"``).
- :func:`audit` records file opens, saves, exports, network use and external tool runs as JSON
  lines, but only when the policy enables it (``"audit_log": true`` or a file path).
- :func:`install_excepthooks` logs uncaught exceptions on any thread in full.

No Qt here.
"""

from __future__ import annotations

import itertools
import json
import logging
import logging.handlers
import os
import sys
import threading
import time
import traceback
from collections import deque
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

CAPACITY = 10_000
LOG_FILE = "chemsuite.log"
AUDIT_FILE = "audit.log"
FILE_BYTES = 2 * 1024 * 1024
FILE_BACKUPS = 4  # plus the live file: five in all
LOG_ENV = "CHEMSUITE_LOG_DIR"
FORMAT = "%(asctime)s %(levelname)-8s %(name)s [%(threadName)s]: %(message)s"
AUDIT_LOGGER = "chemsuite.audit"
LEVELS = {"DEBUG": logging.DEBUG, "INFO": logging.INFO, "WARNING": logging.WARNING,
          "ERROR": logging.ERROR, "CRITICAL": logging.CRITICAL}  # fmt: skip

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Entry:
    seq: int
    time: float  # time.time()
    level: int
    source: str  # the logger's name
    message: str
    details: str = ""  # traceback text
    job: int | None = None  # the job id, when a job wrote it

    @property
    def level_name(self) -> str:
        return logging.getLevelName(self.level)

    def as_text(self) -> str:
        stamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(self.time))
        text = f"{stamp} {self.level_name:<8} {self.source}: {self.message}"
        return text + ("\n" + self.details.rstrip() if self.details else "")


class Logbook:
    """A bounded, thread-safe list of entries, oldest first."""

    def __init__(self, capacity: int = CAPACITY) -> None:
        self._entries: deque[Entry] = deque(maxlen=capacity)
        self._lock = threading.Lock()
        self._seq = itertools.count(1)
        self.cleared_through = 0  # entries up to this seq were cleared by the user

    def add(
        self,
        level: int,
        source: str,
        message: str,
        *,
        details: str = "",
        job: int | None = None,
        when: float | None = None,
    ) -> Entry:
        with self._lock:
            entry = Entry(
                next(self._seq), when or time.time(), level, source, message, details, job
            )
            self._entries.append(entry)
        return entry

    def entries(self, after: int = 0) -> list[Entry]:
        """Entries with ``seq > after``, oldest first."""
        with self._lock:
            if not self._entries or self._entries[-1].seq <= after:
                return []
            if after < self._entries[0].seq:
                return list(self._entries)
            start = len(self._entries) - (self._entries[-1].seq - after)
            return list(itertools.islice(self._entries, start, None))

    @property
    def first_seq(self) -> int:
        with self._lock:
            return self._entries[0].seq if self._entries else 0

    @property
    def last_seq(self) -> int:
        with self._lock:
            return self._entries[-1].seq if self._entries else 0

    def clear(self) -> None:
        with self._lock:
            self.cleared_through = self._entries[-1].seq if self._entries else 0
            self._entries.clear()

    def __len__(self) -> int:
        with self._lock:
            return len(self._entries)


class LogbookHandler(logging.Handler):
    """Feeds log records into a :class:`Logbook`."""

    def __init__(self, book: Logbook, level: int = logging.DEBUG) -> None:
        super().__init__(level)
        self.book = book

    def emit(self, record: logging.LogRecord) -> None:
        try:
            details = ""
            if record.exc_info and record.exc_info[1] is not None:
                details = "".join(traceback.format_exception(*record.exc_info))
            elif record.stack_info:
                details = record.stack_info
            self.book.add(
                record.levelno,
                record.name,
                record.getMessage(),
                details=details,
                job=getattr(record, "job", None),
                when=record.created,
            )
        except Exception:  # noqa: BLE001 - logging must never raise
            self.handleError(record)


_book = Logbook()


def logbook() -> Logbook:
    """The session's logbook."""
    return _book


# -- folders and setup ----------------------------------------------------------------
def log_folder() -> Path:
    """Where log files go: ``CHEMSUITE_LOG_DIR``, else the platform's per-user log folder."""
    if os.environ.get(LOG_ENV):
        return Path(os.environ[LOG_ENV])
    home = Path.home()
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or str(home / "AppData" / "Local")
        return Path(base) / "chemsuite" / "Logs"
    if sys.platform == "darwin":
        return home / "Library" / "Logs" / "chemsuite"
    state = os.environ.get("XDG_STATE_HOME") or str(home / ".local" / "state")
    return Path(state) / "chemsuite" / "log"


def level_from(value: Any, default: int = logging.INFO) -> int:
    """A level from a policy value such as ``"DEBUG"`` or ``10`` (bad values: *default*)."""
    if isinstance(value, int) and not isinstance(value, bool) and value in LEVELS.values():
        return value
    if isinstance(value, str) and value.strip().upper() in LEVELS:
        return LEVELS[value.strip().upper()]
    return default


_installed: list[logging.Handler] = []


def setup_logging(
    *,
    level: int = logging.INFO,
    folder: Path | None = None,
    to_file: bool = True,
    to_stderr: bool = True,
) -> Path | None:
    """Send the root logger to the logbook, stderr and a rotating file; returns the file path.

    Safe to call again (e.g. after the policy is read): earlier handlers are replaced.
    """
    root = logging.getLogger()
    for handler in _installed:
        root.removeHandler(handler)
        handler.close()
    _installed.clear()
    root.setLevel(min(level, logging.INFO))  # the logbook always sees INFO and above
    book = LogbookHandler(_book, logging.DEBUG if level <= logging.DEBUG else logging.INFO)
    _installed.append(book)
    if to_stderr:
        stderr = logging.StreamHandler()
        stderr.setLevel(max(level, logging.INFO))
        stderr.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
        _installed.append(stderr)
    path = None
    if to_file:
        folder = folder or log_folder()
        try:
            folder.mkdir(parents=True, exist_ok=True)
            path = folder / LOG_FILE
            file = logging.handlers.RotatingFileHandler(
                path, maxBytes=FILE_BYTES, backupCount=FILE_BACKUPS, encoding="utf-8", delay=True
            )
            file.setLevel(level)
            file.setFormatter(logging.Formatter(FORMAT))
            _installed.append(file)
        except OSError as err:  # a read-only profile: carry on without a file
            path = None
            log.warning("cannot write log files in %s: %s", folder, err)
    for handler in _installed:
        root.addHandler(handler)
    return path


# -- audit log ----------------------------------------------------------------
_audit = logging.getLogger(AUDIT_LOGGER)
_audit.propagate = False  # never in the ordinary log or on screen
_audit.addHandler(logging.NullHandler())


def setup_audit(setting: Any, folder: Path | None = None) -> Path | None:
    """Turn the audit log on per the policy value: ``True`` (log folder) or a file path."""
    for handler in list(_audit.handlers):
        _audit.removeHandler(handler)
        handler.close()
    _audit.addHandler(logging.NullHandler())
    if not setting:
        _audit.setLevel(logging.CRITICAL + 1)
        return None
    path = Path(setting) if isinstance(setting, str) else (folder or log_folder()) / AUDIT_FILE
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.handlers.RotatingFileHandler(
            path, maxBytes=10 * FILE_BYTES, backupCount=FILE_BACKUPS, encoding="utf-8"
        )
    except OSError as err:
        log.error("cannot write the audit log %s: %s", path, err)
        return None
    handler.setFormatter(logging.Formatter("%(message)s"))
    _audit.addHandler(handler)
    _audit.setLevel(logging.INFO)
    return path


def audit(event: str, **fields: Any) -> None:
    """One audit record (a JSON line), when the audit log is on."""
    if not _audit.isEnabledFor(logging.INFO):
        return
    record = {
        "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "user": _user(),
        "event": event,
        **{k: (str(v) if isinstance(v, Path) else v) for k, v in fields.items()},
    }
    _audit.info(json.dumps(record, default=str, ensure_ascii=False))


def _user() -> str:
    try:
        return os.getlogin()
    except OSError:
        return os.environ.get("USER") or os.environ.get("USERNAME") or "?"


# -- uncaught exceptions ----------------------------------------------------------------
def install_excepthooks(on_error: Callable[[str, str], None] | None = None) -> None:
    """Log uncaught exceptions on every thread; *on_error(summary, details)* is told too."""

    def report(kind: str, exc_type, exc, tb) -> None:
        details = "".join(traceback.format_exception(exc_type, exc, tb))
        logging.getLogger("chemsuite").critical(
            "uncaught exception%s: %s", kind, exc, exc_info=(exc_type, exc, tb)
        )
        if on_error is not None:
            try:
                on_error(f"{exc_type.__name__}: {exc}", details)
            except Exception:
                log.exception("the error reporter failed")

    def hook(exc_type, exc, tb) -> None:
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc, tb)
            return
        report("", exc_type, exc, tb)

    def thread_hook(args: threading.ExceptHookArgs) -> None:
        if args.exc_type is SystemExit:
            return
        name = args.thread.name if args.thread is not None else "?"
        report(f" in thread {name}", args.exc_type, args.exc_value, args.exc_traceback)

    sys.excepthook = hook
    threading.excepthook = thread_hook
