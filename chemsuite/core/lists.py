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

"""Named row lists ("actives", "to make", "sent for assay") saved with the document.

A list is a set of row oids. Deleted rows simply stop counting: a list keeps their oids so
undoing the delete brings them back, and saving drops them. List edits are not table edits
(they are not undoable) but they do make the document need saving; :attr:`RowLists.revision`
counts every change so the owner can compare it with the revision it last saved.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from typing import Any

import numpy as np

from chemsuite.core.table import Table

EXTRAS_KEY = "lists"
ListsListener = Callable[[], None]


class ListError(ValueError):
    """An invalid list name (blank, or taken when renaming)."""


def _as_oids(oids: Iterable[int] | np.ndarray) -> np.ndarray:
    if isinstance(oids, np.ndarray):
        return np.unique(oids.astype(np.int64, copy=False))
    return np.unique(np.fromiter(oids, dtype=np.int64))


class RowLists:
    def __init__(self, table: Table) -> None:
        self._table = table
        self._lists: dict[str, np.ndarray] = {}  # name -> sorted unique oids
        self._listeners: list[ListsListener] = []
        self.revision = 0

    # -- reading ------------------------------------------------------------------------------
    @property
    def names(self) -> list[str]:
        return list(self._lists)

    def __contains__(self, name: object) -> bool:
        return name in self._lists

    def __len__(self) -> int:
        return len(self._lists)

    def oids(self, name: str) -> np.ndarray:
        """The list's rows that exist in the table now (sorted oids)."""
        stored = self._lists[name]
        return stored[np.isin(stored, self._table.oids, assume_unique=True)]

    def count(self, name: str) -> int:
        return len(self.oids(name))

    # -- editing ------------------------------------------------------------------------------
    def add(self, name: str, oids: Iterable[int] | np.ndarray) -> int:
        """Add rows to list *name*, creating it if needed; returns how many were new."""
        name = self._check(name)
        new = _as_oids(oids)
        old = self._lists.get(name)
        if old is None:
            self._lists[name] = new
            added = len(new)
        else:
            merged = np.union1d(old, new)
            added = len(merged) - len(old)
            self._lists[name] = merged
        if old is None or added:
            self._changed()
        return added

    def remove(self, name: str, oids: Iterable[int] | np.ndarray) -> int:
        """Remove rows from list *name* (which stays, even if empty); returns how many."""
        old = self._lists[name]
        kept = np.setdiff1d(old, _as_oids(oids), assume_unique=True)
        removed = len(old) - len(kept)
        if removed:
            self._lists[name] = kept
            self._changed()
        return removed

    def rename(self, old: str, new: str) -> None:
        new = self._check(new)
        if new == old:
            return
        if new in self._lists:
            raise ListError(f"A list named “{new}” already exists.")
        self._lists = {(new if k == old else k): v for k, v in self._lists.items()}
        self._changed()

    def delete(self, name: str) -> None:
        del self._lists[name]
        self._changed()

    @staticmethod
    def _check(name: str) -> str:
        name = name.strip()
        if not name:
            raise ListError("A list needs a name.")
        return name

    # -- persistence --------------------------------------------------------------------------
    def to_json(self) -> list[dict[str, Any]]:
        """For the document's extras; deleted rows are dropped."""
        return [{"name": name, "oids": self.oids(name).tolist()} for name in self._lists]

    def load_json(self, data: Any) -> None:
        """Replace all lists from :meth:`to_json` output (bad entries are skipped)."""
        lists: dict[str, np.ndarray] = {}
        for entry in data if isinstance(data, list) else []:
            if not isinstance(entry, Mapping):
                continue
            name = str(entry.get("name", "")).strip()
            try:
                oids = _as_oids(int(o) for o in entry.get("oids", []))
            except (TypeError, ValueError):
                continue
            if name and name not in lists:
                lists[name] = oids
        self._lists = lists
        self._changed()

    # -- listeners ----------------------------------------------------------------------------
    def subscribe(self, listener: ListsListener) -> Callable[[], None]:
        self._listeners.append(listener)

        def unsubscribe() -> None:
            if listener in self._listeners:
                self._listeners.remove(listener)

        return unsubscribe

    def close(self) -> None:
        self._listeners.clear()

    def _changed(self) -> None:
        self.revision += 1
        for listener in list(self._listeners):
            listener()
