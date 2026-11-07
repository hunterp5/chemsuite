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

"""Driving the AiZynthFinder worker: start it, keep it warm, search target by target.

One worker per (interpreter, config) is kept alive in this process (:func:`worker_for`), so
the ~10–20 s model load is paid once per session instead of once per job. A cancelled search
is the one case where the worker is killed: a tree search cannot be interrupted from outside,
so the process goes and the targets already finished are kept (``Cancelled(partial=…)``).
"""

from __future__ import annotations

import json
import logging
import os
import queue
import subprocess
import threading
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from chemsuite.core.jobs import Cancelled, CancelToken, Progress, check

from . import worker as worker_script
from .config import SearchOptions
from .routes import compact_routes, summary_columns

log = logging.getLogger(__name__)

#: how long to wait for the worker to start, and for one search to come back
START_TIMEOUT_S = 180.0
SEARCH_TIMEOUT_MARGIN_S = 300.0
#: how often a blocking wait looks at the cancel token
POLL_S = 0.05
#: tail of the worker's stderr kept for an error message
_TAIL_CHARS = 4000


class RetroError(RuntimeError):
    """The worker could not be started, died, or answered nothing usable."""


@dataclass(frozen=True)
class SearchReply:
    """What one target came back with. *error* is set when that target failed."""

    target: str
    columns: dict[str, str] = field(default_factory=dict)
    routes: list[dict] = field(default_factory=list)
    error: str = ""

    @property
    def ok(self) -> bool:
        return not self.error


class Worker:
    """One AiZynthFinder process, talking JSON lines (see :mod:`.worker`)."""

    def __init__(self, python: str, config: Path | str) -> None:
        self.python = str(python)
        self.config = str(config)
        self._process: subprocess.Popen[str] | None = None
        self._replies: queue.Queue[dict | None] = queue.Queue()
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
        command = worker_script.command(self.python, self.config)
        try:
            self._process = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                env=child_environment(self.python),
                cwd=str(Path(self.config).parent if Path(self.config).parent.is_dir() else "."),
                **_no_window_kwargs(),
            )
        except OSError as err:
            raise RetroError(f"Could not start {self.python}: {err}") from err
        self._replies = queue.Queue()
        self._errors = []
        threading.Thread(target=self._read_replies, daemon=True).start()
        threading.Thread(target=self._read_errors, daemon=True).start()
        ready = self._await_reply(0, START_TIMEOUT_S, cancel)
        if not ready.get("ready"):
            raise RetroError(self._died_message())

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
            _kill(process)

    # -- requests -----------------------------------------------------------------------------
    def request(
        self, op: str, body: dict | None = None, *, timeout: float, cancel: CancelToken | None
    ) -> dict:
        """Send one request and return its reply. Raises :class:`RetroError` on a dead worker.

        A reply with ``ok: false`` is returned as-is: a target that fails is the caller's to
        report, not a reason to tear the worker down.
        """
        with self._lock:
            self.start(cancel)
            process = self._process
            if process is None or process.stdin is None:
                raise RetroError("The AiZynthFinder worker is not running.")
            ident = self._next_id
            self._next_id += 1
            payload = {"id": ident, "op": op, **(body or {})}
            try:
                process.stdin.write(json.dumps(payload) + "\n")
                process.stdin.flush()
            except (OSError, ValueError) as err:
                raise RetroError(self._died_message()) from err
            return self._await_reply(ident, timeout, cancel)

    def hello(self, cancel: CancelToken | None = None) -> dict:
        """The worker's AiZynthFinder version. Loads no models."""
        reply = self.request("hello", timeout=START_TIMEOUT_S, cancel=cancel)
        if not reply.get("ok"):
            raise RetroError(str(reply.get("error") or "AiZynthFinder could not be imported."))
        return reply

    def search_one(
        self, target: str, options: SearchOptions, cancel: CancelToken | None = None
    ) -> SearchReply:
        """Search one target. A failure for this target comes back in ``error``."""
        body = {"target": target, **options.request()}
        timeout = float(options.time_limit) + SEARCH_TIMEOUT_MARGIN_S
        reply = self.request("search", body, timeout=timeout, cancel=cancel)
        if not reply.get("ok"):
            log.warning("retrosynthesis failed for %s: %s", target, reply.get("traceback") or "")
            return SearchReply(target, error=_message(reply))
        stats = reply.get("stats") if isinstance(reply.get("stats"), dict) else {}
        return SearchReply(
            target=str(reply.get("target") or target),
            columns=summary_columns(stats),
            routes=compact_routes(reply.get("trees"), stats.get("stock_info")),
        )

    # -- plumbing -----------------------------------------------------------------------------
    def _read_replies(self) -> None:
        process = self._process
        if process is None or process.stdout is None:
            return
        try:
            for line in process.stdout:
                text = line.strip()
                if not text:
                    continue
                try:
                    self._replies.put(json.loads(text))
                except ValueError:
                    log.debug("ignored a line that is not a reply: %.200s", text)
        except (OSError, ValueError):
            pass
        finally:
            self._replies.put(None)  # end of stream: the worker is gone

    def _read_errors(self) -> None:
        process = self._process
        if process is None or process.stderr is None:
            return
        try:
            for line in process.stderr:
                self._errors.append(line)
                if len(self._errors) > 400:
                    del self._errors[:200]
        except (OSError, ValueError):
            pass

    def _await_reply(self, ident: int, timeout: float, cancel: CancelToken | None) -> dict:
        """Wait for the reply to *ident*, looking at *cancel* as it goes."""
        waited = 0.0
        while True:
            check(cancel)
            try:
                reply = self._replies.get(timeout=POLL_S)
            except queue.Empty:
                waited += POLL_S
                if waited >= timeout:
                    raise RetroError(
                        f"AiZynthFinder did not answer within {timeout:.0f} s."
                    ) from None
                continue
            if reply is None:
                raise RetroError(self._died_message())
            if reply.get("id") == ident:
                return reply
            log.debug("ignored a stale reply %r", reply.get("id"))

    def _died_message(self) -> str:
        tail = "".join(self._errors)[-_TAIL_CHARS:].strip()
        code = self._process.poll() if self._process is not None else None
        head = "The AiZynthFinder worker stopped"
        if code not in (None, 0):
            head += f" (exit code {code})"
        return f"{head}.\n{tail}" if tail else f"{head}."


# -- the warm worker ----------------------------------------------------------------------------
# One at a time: a loaded AiZynthFinder holds its models and the whole stock in memory (about a
# gigabyte with the ZINC stock), so changing either path in Settings retires the old process
# rather than leaving it loaded beside the new one.
_WORKERS: dict[tuple[str, str], Worker] = {}
_WORKERS_LOCK = threading.Lock()


def worker_for(python: str, config: Path | str) -> Worker:
    """The warm worker for this interpreter and config, started on first use."""
    key = (str(python), str(config))
    with _WORKERS_LOCK:
        stale = [other for name, other in _WORKERS.items() if name != key]
        for name in [name for name in _WORKERS if name != key]:
            del _WORKERS[name]
        found = _WORKERS.get(key)
        if found is None:
            found = Worker(*key)
            _WORKERS[key] = found
    for other in stale:
        other.close()
    return found


def close_workers() -> None:
    """Stop every warm worker (the app closing, or the paths changing in Settings)."""
    with _WORKERS_LOCK:
        workers = list(_WORKERS.values())
        _WORKERS.clear()
    for item in workers:
        item.close()


def search(
    python: str,
    config: Path | str,
    targets: Sequence[str],
    options: SearchOptions,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    reuse: bool = True,
) -> list[SearchReply]:
    """Search each of *targets* in turn, one reply per target, in order.

    Cancelling raises ``Cancelled`` carrying the replies already collected, so the rows that
    finished are written (the project rule for every per-row tool). The worker is killed on
    the way out, since a running tree search cannot be stopped any other way.
    """
    worker = worker_for(python, config) if reuse else Worker(python, config)
    done: list[SearchReply] = []
    total = len(targets)
    _report(progress, 0, total)
    try:
        for index, target in enumerate(targets, start=1):
            done.append(worker.search_one(target, options, cancel))
            _report(progress, index, total)
    except Cancelled:
        worker.close()
        raise Cancelled(done) from None
    except RetroError:
        worker.close()
        if not done:
            raise
        log.exception("the retrosynthesis worker stopped after %d of %d targets", len(done), total)
        return done
    finally:
        if not reuse:
            worker.close()
    return done


def probe(python: str, config: Path | str | None) -> tuple[bool, str]:
    """``(ok, message)`` for the Settings ▸ External Tools **Test** button.

    Checks the path, refuses chemsuite's own interpreter, then asks a throwaway worker for the
    AiZynthFinder version: no models are loaded and no warm worker is left behind.
    """
    from .config import config_problem, is_app_python

    if not python:
        return False, "Set the interpreter of the environment AiZynthFinder is installed in."
    if is_app_python(python):
        return False, "That is chemsuite's own interpreter; AiZynthFinder needs its own."
    if config is not None:
        path = Path(config)
        if not path.is_file():
            return False, "No such file."
        try:
            problem = config_problem(path.read_text(encoding="utf-8"))
        except OSError as err:
            return False, f"Could not read {path.name}: {err}"
        if problem:
            return False, problem
    worker = Worker(python, config if config is not None else "")
    try:
        reply = worker.hello()
    except RetroError as err:
        return False, str(err).splitlines()[0]
    finally:
        worker.close()
    version = str(reply.get("version") or "?")
    where = f", config {Path(config).name}" if config is not None else ""
    return True, f"AiZynthFinder {version} on Python {reply.get('python', '?')}{where}."


def unique_targets(values: Iterable[str]) -> list[str]:
    """The distinct, non-blank targets in order: the same SMILES is searched once."""
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = (value or "").strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def child_environment(python: str) -> dict[str, str]:
    """The worker's environment: chemsuite's own Python settings removed.

    A conda environment gets its own ``bin``/``Library\\bin`` in front of ``PATH`` so that the
    RDKit and HDF5 libraries beside its interpreter load without activating the environment.
    """
    env = os.environ.copy()
    for key in (
        "PYTHONPATH",
        "PYTHONHOME",
        "CONDA_PREFIX",
        "CONDA_DEFAULT_ENV",
        "CONDA_PYTHON_EXE",
    ):
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    root = _conda_root(python)
    if root is None:
        return env
    extra = [root, root / "Library" / "bin", root / "Scripts", root / "bin"]
    prefix = os.pathsep.join(str(path) for path in extra if path.is_dir())
    if prefix:
        env["PATH"] = prefix + os.pathsep + env.get("PATH", "")
    env["CONDA_PREFIX"] = str(root)
    return env


def _conda_root(python: str) -> Path | None:
    """The environment root of *python* when it looks like a conda or venv layout."""
    path = Path(python)
    for parent in path.parents:
        if (parent / "conda-meta").is_dir() or (parent / "pyvenv.cfg").is_file():
            return parent
    return None


def _no_window_kwargs() -> dict:
    """Keep Windows from flashing a console window for the worker."""
    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    return {"creationflags": flags} if flags else {}


def _kill(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    try:
        process.kill()
        process.wait(timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        log.warning("could not stop the AiZynthFinder worker (pid %s)", process.pid)


def _message(reply: dict) -> str:
    text = str(reply.get("error") or "").strip()
    return text or "AiZynthFinder reported no result for this structure."


def _report(progress: Progress | None, done: int, total: int) -> None:
    if progress is not None:
        progress(done, total)
