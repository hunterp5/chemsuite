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
"""The workspace session's events as Qt signals, for one open viewer window.

The session (:class:`~chemsuite.protein.workspace.session.WorkbenchSession`) is Qt-free and
lives with the document; a window's widgets connect to these signals instead, so Qt drops each
connection when its widget goes, and :meth:`SessionSignals.close` disconnects the window from
the session when the window closes.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QObject, Signal

from chemsuite.protein.workspace.session import EVENTS, WorkbenchSession


class SessionSignals(QObject):
    proteinsChanged = Signal()
    sceneChanged = Signal()
    lookChanged = Signal()
    selectionChanged = Signal()
    hovered = Signal(str, int)
    focusChanged = Signal()
    measurementsChanged = Signal()
    labelsChanged = Signal()
    edited = Signal()

    def __init__(self, session: WorkbenchSession, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.session = session
        relays = {
            "proteinsChanged": (session.proteinsChanged, self.proteinsChanged),
            "sceneChanged": (session.sceneChanged, self.sceneChanged),
            "lookChanged": (session.lookChanged, self.lookChanged),
            "selectionChanged": (session.selectionChanged, self.selectionChanged),
            "hovered": (session.hovered, self.hovered),
            "focusChanged": (session.focusChanged, self.focusChanged),
            "measurementsChanged": (session.measurementsChanged, self.measurementsChanged),
            "labelsChanged": (session.labelsChanged, self.labelsChanged),
            "edited": (session.edited, self.edited),
        }
        assert set(relays) == set(EVENTS), "every session event has a signal"
        self._disconnect: list[Callable[[], None]] = [
            event.connect(signal.emit) for event, signal in relays.values()
        ]

    def close(self) -> None:
        """Stop relaying: the window is closing (the session lives on)."""
        for disconnect in self._disconnect:
            disconnect()
        self._disconnect = []
