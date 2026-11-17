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

import logging
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from chemsuite.core.jobs import Cancelled, CancelToken, Progress
from chemsuite.core.jobs.lines import LineWorker, WarmWorkers, WorkerError, child_environment

from . import worker as worker_script
from .config import SearchOptions
from .routes import compact_routes, summary_columns

log = logging.getLogger(__name__)

#: how long to wait for the worker to start, and for one search to come back
START_TIMEOUT_S = 180.0
SEARCH_TIMEOUT_MARGIN_S = 300.0


class RetroError(WorkerError):
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


class Worker(LineWorker):
    """One AiZynthFinder process, talking JSON lines (see :mod:`.worker`)."""

    error = RetroError

    def __init__(self, python: str, config: Path | str) -> None:
        self.python = str(python)
        self.config = str(config)
        folder = Path(self.config).parent
        super().__init__(
            worker_script.command(self.python, self.config),
            name="The AiZynthFinder worker",
            env=child_environment(self.python),
            cwd=str(folder if folder.is_dir() else "."),
            start_timeout=START_TIMEOUT_S,
        )

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


# -- the warm worker ----------------------------------------------------------------------------
# One at a time: a loaded AiZynthFinder holds its models and the whole stock in memory (about a
# gigabyte with the ZINC stock), so changing either path in Settings retires the old process
# rather than leaving it loaded beside the new one.
_WORKERS: WarmWorkers[Worker] = WarmWorkers(lambda key: Worker(*key), single=True)


def worker_for(python: str, config: Path | str) -> Worker:
    """The warm worker for this interpreter and config, started on first use."""
    return _WORKERS.get((str(python), str(config)))


def close_workers() -> None:
    """Stop every warm worker (the app closing, or the paths changing in Settings)."""
    _WORKERS.close_all()


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
    """``(ok, message)`` for the Settings ▸ External **Test** button.

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


def _message(reply: dict) -> str:
    text = str(reply.get("error") or "").strip()
    return text or "AiZynthFinder reported no result for this structure."


def _report(progress: Progress | None, done: int, total: int) -> None:
    if progress is not None:
        progress(done, total)
