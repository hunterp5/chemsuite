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
"""Events without Qt: what the workspace session tells whoever draws it.

An :class:`Event` is a list of callbacks. The Qt workbench connects each one to a Qt signal of
the same name (``protein/ui/session.py``), so widgets connect to signals that Qt drops when the
widget goes, and the window disconnects from the session when it closes, while the session
lives on with the document.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any


class Event:
    def __init__(self) -> None:
        self._callbacks: list[Callable[..., Any]] = []

    def connect(self, callback: Callable[..., Any]) -> Callable[[], None]:
        """Call *callback* on every :meth:`emit`; returns a function that disconnects it."""
        self._callbacks.append(callback)

        def disconnect() -> None:
            if callback in self._callbacks:
                self._callbacks.remove(callback)

        return disconnect

    def emit(self, *args: Any) -> None:
        for callback in list(self._callbacks):
            callback(*args)

    def __len__(self) -> int:
        return len(self._callbacks)
