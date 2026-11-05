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

"""Background work without Qt: cancellation, progress reporting and process-pool helpers.

A job is a plain function that accepts optional ``progress`` and ``cancel`` arguments. The Qt
job runner (``chemsuite.app``) decides which thread or process runs it; the job itself only
reports and checks.
"""

from __future__ import annotations

import threading
from collections.abc import Callable

#: ``progress(done, total)``; *total* may be 0 when unknown.
Progress = Callable[[int, int], None]


class Cancelled(Exception):
    """Raised inside a job when its :class:`CancelToken` has been cancelled.

    A row-by-row computation may carry what it finished before stopping as *partial* (results
    for the first rows, in order), so the app can keep that work instead of discarding it.
    """

    def __init__(self, partial: object = None) -> None:
        super().__init__()
        self.partial = partial


class CancelToken:
    """Thread-safe cancellation flag shared between the GUI and a running job."""

    def __init__(self) -> None:
        self._event = threading.Event()

    def cancel(self) -> None:
        self._event.set()

    @property
    def cancelled(self) -> bool:
        return self._event.is_set()

    def raise_if_cancelled(self) -> None:
        if self._event.is_set():
            raise Cancelled


def check(cancel: CancelToken | None) -> None:
    """Raise :class:`Cancelled` if *cancel* is set (a no-op for None)."""
    if cancel is not None:
        cancel.raise_if_cancelled()
