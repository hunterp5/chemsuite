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
"""Building the slow parts of shown proteins off the GUI thread.

Surfaces, interactions and the cartoons of large proteins are built as quick background jobs
(see :class:`~chemsuite.protein.look.display.Display` and its ``pending``/``accept``). Until a part
arrives the view keeps drawing the last one, so nothing flickers.

A result is used only if, when it arrives, the session still shows that display, the display
still holds the structure the work was handed (an edit replaces it), and the part has not
changed since (its version). Anything else is dropped: the work for the newer version is
already on its way. Asking again for the version being built does nothing; a newer version
cancels the older job, which never runs if it has not started.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import shiboken6
from PySide6.QtCore import QObject

from chemsuite.protein.look.display import Build, Display

if TYPE_CHECKING:
    from chemsuite.app.shell.jobs import Job

    from .session import ProteinSession

#: runs ``work(job)`` in the background and ``on_done(result)`` on the GUI thread
Run = Callable[[str, Callable[[Any], Any], Callable[[Any], None]], "Job | None"]
TITLES = {
    "surface": "Building a protein surface",
    "pocket": "Building a pocket surface",
    "interactions": "Finding interactions",
    "cartoon": "Building a cartoon",
}


class Builder(QObject):
    def __init__(self, session: ProteinSession, run: Run) -> None:
        super().__init__(session)
        self.session = session
        self._run = run
        self._running: dict[tuple[str, str], tuple[Build, object, Job]] = {}

    def schedule(self, shown: tuple[Display, ...] | list[Display]) -> None:
        """Start the work the shown displays are waiting for."""
        for display in shown:
            for build in display.pending():
                self._start(display, build)

    def busy(self) -> bool:
        return bool(self._running)

    def cancel_all(self) -> None:
        for _build, _structure, job in list(self._running.values()):
            job.cancel()
        self._running.clear()

    def _start(self, display: Display, build: Build) -> None:
        key = (display.id, build.part)
        structure = display.structure
        running = self._running.get(key)
        if running is not None:
            old, old_structure, job = running
            if old.version == build.version and old_structure is structure:
                return
            job.cancel()
        alive = shiboken6.isValid  # the viewer can close while the job runs

        def done(result: Any) -> None:
            if alive(self):
                self._finished(key, display, structure, build, result)

        job = self._run(TITLES[build.part], lambda _job: build(), done)
        if job is None:
            return
        self._running[key] = (build, structure, job)
        job.done.connect(lambda: alive(self) and self._forget(key, job))

    def _forget(self, key: tuple[str, str], job: Job) -> None:
        running = self._running.get(key)
        if running is not None and running[2] is job:
            del self._running[key]

    def _finished(
        self, key: tuple[str, str], display: Display, structure, build: Build, result: Any
    ) -> None:
        running = self._running.get(key)
        if running is not None and running[0] is build:
            del self._running[key]
        if display not in self.session.all_displays() or display.structure is not structure:
            return
        if display.accept(build.part, build.version, result):
            self.session.built()
