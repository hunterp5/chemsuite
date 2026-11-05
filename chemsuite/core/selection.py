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

"""The selected rows of one table: the single source of truth for every view.

The table view, plots and the browser all read and write this store. Each change records its
*origin* (e.g. ``"table"``, ``"plot:3"``) so a view can ignore the echo of a change it made
itself instead of redrawing. Rows deleted from the table leave the selection automatically.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np

from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev

SelectionListener = Callable[[frozenset[int], str], None]


class Selection:
    def __init__(self, table: Table) -> None:
        self._table = table
        self._oids: frozenset[int] = frozenset()
        self._listeners: list[SelectionListener] = []
        self._unsubscribe = table.subscribe(self._on_table_event)

    @property
    def oids(self) -> frozenset[int]:
        return self._oids

    def __len__(self) -> int:
        return len(self._oids)

    def __contains__(self, oid: object) -> bool:
        return oid in self._oids

    def as_array(self) -> np.ndarray:
        """Selected oids in table storage order (ascending)."""
        return np.fromiter(sorted(self._oids), dtype=np.int64, count=len(self._oids))

    def set(self, oids: Iterable[int], origin: str) -> None:
        new = frozenset(int(o) for o in oids)
        if new != self._oids:
            self._oids = new
            self._notify(origin)

    def add(self, oids: Iterable[int], origin: str) -> None:
        self.set(self._oids | {int(o) for o in oids}, origin)

    def remove(self, oids: Iterable[int], origin: str) -> None:
        self.set(self._oids - {int(o) for o in oids}, origin)

    def clear(self, origin: str) -> None:
        self.set((), origin)

    def invert(self, origin: str) -> None:
        """Select every row (including filtered-out rows) that is not selected now."""
        self.set(set(self._table.oids.tolist()) - self._oids, origin)

    def subscribe(self, listener: SelectionListener) -> Callable[[], None]:
        self._listeners.append(listener)

        def unsubscribe() -> None:
            if listener in self._listeners:
                self._listeners.remove(listener)

        return unsubscribe

    def close(self) -> None:
        self._unsubscribe()
        self._listeners.clear()

    def _notify(self, origin: str) -> None:
        for listener in list(self._listeners):
            listener(self._oids, origin)

    def _on_table_event(self, event: ev.TableEvent) -> None:
        if not isinstance(event, Changed) or not self._oids:
            return
        removed = [c.oids for c in event.changes if isinstance(c, ev.RowsRemoved)]
        if removed:
            gone = set(np.concatenate(removed).tolist())
            if gone & self._oids:
                self.set(self._oids - gone, "table-edit")
