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

"""What a plot computes before it can draw (cores, networks, trees), found in the background.

A :class:`Finder` keys each request by its settings and rows. A request already answered is
answered again at once; a new one starts a job and :attr:`Finder.ready` says when it is done.
A newer request makes an older job's result go unheard.
"""

from __future__ import annotations

from collections.abc import Callable, Hashable
from concurrent.futures import Executor

import numpy as np
from PySide6.QtCore import QObject, Signal

from chemsuite.core.table import Table

from .spec import PlotSpec

#: a job runner: ``run(label, work(job), on_done(result))``
Runner = Callable[[str, Callable, Callable], object]
#: what a job returns: ``(result, problem, learn)``; *learn* runs on the GUI thread
Outcome = tuple[object, str, Callable[[], None] | None]


class Finder(QObject):
    """Computes a plot's data in the background when it can run jobs (without a runner, at
    once). Subclasses give :meth:`key` and :meth:`work`."""

    ready = Signal()
    waiting = "Working…"

    def __init__(
        self,
        table: Table,
        run: Runner | None = None,
        workers: Callable[[int], Executor | None] | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self.table = table
        self._run = run
        self._workers = workers
        self._key: Hashable | None = None
        self._result: object = None
        self.message = self.waiting

    def key(self, spec: PlotSpec, oids: np.ndarray) -> Hashable | None:
        """What the result depends on; None (with :attr:`message` saying why) when it cannot
        be found."""
        raise NotImplementedError

    def work(self, spec: PlotSpec, oids: np.ndarray) -> Callable[[object], Outcome]:
        """The job for the request :meth:`key` just keyed (read the table here, on the GUI
        thread; the job it returns runs in the background)."""
        raise NotImplementedError

    def job_title(self, spec: PlotSpec) -> str:
        """The job's name in the status bar and the Log: the plot being calculated."""
        return "Calculating plot…"

    def executor(self, rows: int) -> Executor | None:
        return self._workers(rows) if self._workers is not None else None

    def result(self, spec: PlotSpec, oids: np.ndarray) -> object:
        """The result for *spec* over *oids*; None while it is being found (then
        :attr:`ready` says when) or when it cannot be, :attr:`message` saying why."""
        key = self.key(spec, oids)
        if key is None:
            self._key = self._result = None
            return None
        if key == self._key:
            return self._result
        self._key, self._result, self.message = key, None, self.waiting
        job = self.work(spec, oids)

        def run(handle=None) -> tuple:
            return (key, *job(handle))

        if self._run is None:
            self._done(run())
        else:
            self._run(self.job_title(spec), run, self._done)
        return self._result

    def close(self) -> None:
        """The plot closed: a job still running finishes unheard."""
        self._key = None

    def _done(self, outcome: tuple) -> None:
        key, result, problem, learn = outcome
        if learn is not None:
            learn()
        if key != self._key:
            return  # a newer request is under way
        self._result = result
        self.message = f"⚠ {problem}" if problem else ""
        self.ready.emit()
