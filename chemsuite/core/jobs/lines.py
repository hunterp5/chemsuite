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

"""A worker in another interpreter, spoken to in JSON lines.

Programs that cannot share chemsuite's interpreter (AiZynthFinder, OpenMM) run a small script
in their own, kept warm between jobs. This is the client both use.

**Protocol.** One JSON object per line each way.

- The worker's first line is ``{"id": 0, "ok": true, "ready": true}``.
- A request is ``{"id": n, "op": "...", ...}``; its reply carries the same ``id`` and
  ``"ok": true`` or ``"ok": false, "error": "..."``. A failed request leaves the worker up.
- Before its reply, a request may send any number of **events**,
  ``{"id": n, "event": {...}}`` (progress, energies, frames).
- ``{"id": n, "op": "stop"}`` asks the running request *n* to end early, at a point where
  what it has done is whole; it still replies. A worker that cannot stop ignores it.
- ``{"id": -1, "op": "quit"}`` ends the worker, as does closing its stdin.

A cancelled request is stopped this way when the caller says it can be
(``request(..., stoppable=True)``), and the worker is killed if it does not reply within
:data:`STOP_GRACE_S`; otherwise cancelling kills the worker at once.
"""

from __future__ import annotations

import json
import logging
import os
import queue
import subprocess
import sys
import threading
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Generic, TypeVar

from . import Cancelled, CancelToken, children

log = logging.getLogger(__name__)

#: how often a blocking wait looks at the cancel token
POLL_S = 0.05
#: how long a stopped request has to reply before its worker is killed
STOP_GRACE_S = 30.0
#: tail of the worker's stderr kept for an error message
TAIL_CHARS = 4000

Event = Callable[[dict], None]


class WorkerError(RuntimeError):
    """The worker could not be started, died, or answered nothing usable."""


class LineWorker:
    """One worker process. Requests are taken one at a time (a lock serialises them)."""

    #: raised for a worker that is gone or silent; a subclass names its own
    error: type[WorkerError] = WorkerError

    def __init__(
        self,
        command: Sequence[str],
        *,
        name: str = "The worker",
        env: dict[str, str] | None = None,
        cwd: str | None = None,
        start_timeout: float = 180.0,
    ) -> None:
        self.command = [str(part) for part in command]
        self.name = name
        self.env = env
        self.cwd = cwd
        self.start_timeout = start_timeout
        self._process: subprocess.Popen[str] | None = None
        self._lines: queue.Queue[dict | None] = queue.Queue()
        self._errors: list[str] = []
        self._next_id = 1
        self._lock = threading.Lock()

    # -- the process --------------------------------------------------------------------------
    @property
    def alive(self) -> bool:
        return self._process is not None and self._process.poll() is None

    def start(self, cancel: CancelToken | None = None) -> None:
        """Start the worker and wait for its ready line. A live worker is left alone."""
        if self.alive:
            return
        self.close()
        try:
            self._process = children.start(
                self.command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                env=self.env,
                cwd=self.cwd,
            )
        except OSError as err:
            raise self.error(f"Could not start {self.command[0]}: {err}") from err
        # each process gets its own queue and error list, handed to its readers: a killed
        # worker's last lines (its end of stream above all) never reach its successor
        self._lines = queue.Queue()
        self._errors = []
        process = self._process
        threading.Thread(target=_read_lines, args=(process, self._lines), daemon=True).start()
        threading.Thread(target=_read_errors, args=(process, self._errors), daemon=True).start()
        ready = self._await(0, self.start_timeout, cancel, None, stoppable=False)
        if not ready.get("ready"):
            raise self.error(self._died_message())

    def kill(self) -> None:
        """End the worker at once (a cancelled request it cannot stop: it is not listening)."""
        process = self._process
        self._process = None
        if process is not None:
            children.stop(process)

    def close(self) -> None:
        """Ask the worker to quit; kill it if it will not."""
        process = self._process
        self._process = None
        if process is None:
            return
        try:
            if process.poll() is None and process.stdin is not None:
                process.stdin.write(json.dumps({"id": -1, "op": "quit"}) + "\n")
                process.stdin.flush()
                process.wait(timeout=5)
        except (OSError, ValueError, subprocess.TimeoutExpired):
            pass
        finally:
            children.stop(process)

    # -- requests -----------------------------------------------------------------------------
    def request(
        self,
        op: str,
        body: dict | None = None,
        *,
        timeout: float,
        cancel: CancelToken | None = None,
        events: Event | None = None,
        stoppable: bool = False,
    ) -> dict:
        """Send one request and return its reply (``ok: false`` replies too: the caller's to
        report). *timeout* is how long the worker may stay silent: every event restarts it.
        *events* is called with each event's body, on this thread.

        Cancelling raises :class:`Cancelled`. A *stoppable* request is asked to stop first and
        its reply is the exception's ``partial``; otherwise, or if it does not reply in
        :data:`STOP_GRACE_S`, the worker is killed and ``partial`` is None.
        """
        with self._lock:
            self.start(cancel)
            ident = self._next_id
            self._next_id += 1
            self._send({"id": ident, "op": op, **(body or {})})
            return self._await(ident, timeout, cancel, events, stoppable=stoppable)

    # -- plumbing -----------------------------------------------------------------------------
    def _send(self, payload: dict) -> None:
        process = self._process
        if process is None or process.stdin is None:
            raise self.error(f"{self.name} is not running.")
        try:
            process.stdin.write(json.dumps(payload) + "\n")
            process.stdin.flush()
        except (OSError, ValueError) as err:
            raise self.error(self._died_message()) from err

    def _await(
        self,
        ident: int,
        timeout: float,
        cancel: CancelToken | None,
        events: Event | None,
        *,
        stoppable: bool,
    ) -> dict:
        """Wait for the reply to *ident*, passing its events on and watching *cancel*."""
        last = time.monotonic()
        stopped_at: float | None = None
        while True:
            if cancel is not None and cancel.cancelled and stopped_at is None:
                if not stoppable:
                    self.kill()
                    raise Cancelled
                self._send({"id": ident, "op": "stop"})
                stopped_at = time.monotonic()
            if stopped_at is not None and time.monotonic() - stopped_at > STOP_GRACE_S:
                self.kill()
                raise Cancelled
            try:
                line = self._lines.get(timeout=POLL_S)
            except queue.Empty:
                if time.monotonic() - last >= timeout:
                    silent = f"{self.name} did not answer within {timeout:.0f} s."
                    raise self.error(silent) from None
                continue
            if line is None:
                if stopped_at is not None:
                    raise Cancelled
                raise self.error(self._died_message())
            if line.get("id") != ident:
                log.debug("ignored a stale line %r", line.get("id"))
                continue
            last = time.monotonic()
            event = line.get("event")
            if isinstance(event, dict):
                if events is not None:
                    events(event)
                continue
            if stopped_at is not None:
                raise Cancelled(line)
            return line

    def _died_message(self) -> str:
        """What went wrong, with the exit code and the end of stderr once the process is gone
        (its stdout closing comes a moment before it can be reaped)."""
        process = self._process
        code = None
        if process is not None:
            try:
                code = process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                code = process.poll()
        time.sleep(POLL_S)  # the stderr reader's last lines
        tail = "".join(self._errors)[-TAIL_CHARS:].strip()
        head = f"{self.name} stopped"
        if code not in (None, 0):
            head += f" (exit code {code})"
        return f"{head}.\n{tail}" if tail else f"{head}."


def _read_lines(process: subprocess.Popen[str], into: queue.Queue[dict | None]) -> None:
    """Every JSON object *process* writes, into *into*; None when its stdout closes."""
    try:
        if process.stdout is None:
            return
        for raw in process.stdout:
            text = raw.strip()
            if not text:
                continue
            try:
                found = json.loads(text)
            except ValueError:
                log.debug("ignored a line that is not JSON: %.200s", text)
                continue
            if isinstance(found, dict):
                into.put(found)
    except (OSError, ValueError):
        pass
    finally:
        into.put(None)  # end of stream: the worker is gone


def _read_errors(process: subprocess.Popen[str], into: list[str]) -> None:
    """What *process* writes to stderr, the last few hundred lines kept."""
    if process.stderr is None:
        return
    try:
        for line in process.stderr:
            into.append(line)
            if len(into) > 400:
                del into[:200]
    except (OSError, ValueError):
        pass


W = TypeVar("W", bound=LineWorker)


class WarmWorkers(Generic[W]):
    """The warm workers of one kind, by key, started on first use.

    *single* keeps one at a time: a new key retires the others (a worker that holds gigabytes
    of models should not linger beside its replacement).
    """

    def __init__(self, make: Callable[[Any], W], *, single: bool = False) -> None:
        self._make = make
        self._single = single
        self._workers: dict[Any, W] = {}
        self._lock = threading.Lock()

    def get(self, key: Any) -> W:
        with self._lock:
            stale = []
            if self._single:
                stale = [w for k, w in self._workers.items() if k != key]
                self._workers = {k: w for k, w in self._workers.items() if k == key}
            found = self._workers.get(key)
            if found is None:
                found = self._make(key)
                self._workers[key] = found
        for other in stale:
            other.close()
        return found

    def close_all(self) -> None:
        """Stop every warm worker (the app closing, or the paths changing in Settings)."""
        with self._lock:
            workers = list(self._workers.values())
            self._workers.clear()
        for worker in workers:
            worker.close()


# -- the interpreter ----------------------------------------------------------------------------
#: chemsuite's own Python settings, which a worker's interpreter must not see
_OWN_SETTINGS = (
    "PYTHONPATH",
    "PYTHONHOME",
    "CONDA_PREFIX",
    "CONDA_DEFAULT_ENV",
    "CONDA_PYTHON_EXE",
)


def child_environment(python: str) -> dict[str, str]:
    """The environment for a worker run by *python*: chemsuite's own Python settings removed.

    A conda environment gets its own folders (``Library\\bin``, ``Scripts``, ``bin``) in front of
    ``PATH``, so the libraries beside its interpreter load without activating it: RDKit's and
    HDF5's, and on Windows OpenMM's CUDA plugin, which fails with error 126 otherwise.
    """
    env = os.environ.copy()
    for key in _OWN_SETTINGS:
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    root = conda_root(python)
    if root is None:
        return env
    extra = [root, root / "Library" / "bin", root / "Scripts", root / "bin"]
    prefix = os.pathsep.join(str(path) for path in extra if path.is_dir())
    if prefix:
        env["PATH"] = prefix + os.pathsep + env.get("PATH", "")
    env["CONDA_PREFIX"] = str(root)
    return env


def conda_root(python: str) -> Path | None:
    """The environment root of *python* when it looks like a conda or venv layout."""
    for parent in Path(python).parents:
        if (parent / "conda-meta").is_dir() or (parent / "pyvenv.cfg").is_file():
            return parent
    return None


def is_app_python(python: str) -> bool:
    """True when *python* is the interpreter running chemsuite itself."""
    try:
        return Path(python).resolve() == Path(sys.executable).resolve()
    except OSError:
        return False
