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
"""Everything that belongs to one open table: the table, its selection, model and row lists,
and the document parts tools keep with it.

A **document part** is state a tool keeps with the document rather than in a window: it is
saved with the document under its own key whether or not the tool's window is open, and its
``revision`` tells the document it needs saving. The protein workspace is one
(:mod:`chemsuite.tools.protein`). The app never knows what a part is: what the file holds
under a key no part has claimed yet is saved back as it was.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any, Protocol, TypeVar, cast

from chemsuite.core.lists import EXTRAS_KEY as LISTS_KEY
from chemsuite.core.lists import RowLists
from chemsuite.core.selection import Selection
from chemsuite.core.table import Table

from ..table.table_model import TableModel
from . import session

T = TypeVar("T")


class DocumentPart(Protocol):
    """State a tool keeps with a document (see the module docstring)."""

    @property
    def revision(self) -> int:
        """Changes so far that are worth saving."""
        ...

    def capture(self) -> Any:
        """What the document saves (JSON-friendly; empty: nothing)."""
        ...

    def close(self) -> None:
        """The document is closing."""
        ...


class Workspace:
    """Everything that belongs to one open table."""

    def __init__(
        self,
        table: Table,
        path: Path | None = None,
        extras: dict[str, Any] | None = None,
        *,
        unsaved: bool = False,
    ) -> None:
        self.table = table
        #: new content that was never saved (e.g. rows copied from another document)
        self._unsaved = unsaved
        self.path = path
        self.selection = Selection(table)
        self.model = TableModel(table)
        self.lists = RowLists(table)
        if extras and LISTS_KEY in extras:
            self.lists.load_json(extras[LISTS_KEY])
        #: the window session saved with the document (see ``session.PARTS``), put back when
        #: the window installs this workspace
        self.saved_session: dict[str, Any] = {
            key: value for key, value in (extras or {}).items() if key in session.KEYS
        }
        self._saved_lists = self.lists.revision
        #: what the file holds that neither the window nor the lists read: the saved state of
        #: document parts (:meth:`part`) not made yet, saved back as it was until they are
        self._kept: dict[str, Any] = {
            key: value
            for key, value in (extras or {}).items()
            if key not in session.KEYS and key != LISTS_KEY
        }
        self._parts: dict[str, DocumentPart] = {}
        self._saved_parts: dict[str, int] = {}

    @property
    def is_empty(self) -> bool:
        return not (len(self.table) or self.table.column_names or self.table.assets)

    @property
    def is_modified(self) -> bool:
        """Unsaved table edits or list changes, or content that was never saved."""
        return (
            self._unsaved
            or not self.table.history.is_clean
            or self.lists.revision != self._saved_lists
            or self.part_revisions() != self._saved_parts
        )

    def part(self, key: str, make: Callable[[Any], T]) -> T:
        """The document's part saved under *key* (a tool's state that lives with the document
        rather than a window: the protein workspace), made by ``make(saved)`` when first asked
        for, *saved* being what the file held under *key* (None if nothing)."""
        found = self._parts.get(key)
        if found is None:
            found = self._parts[key] = make(self._kept.pop(key, None))
            self._saved_parts[key] = found.revision  # a new part has nothing unsaved
        return cast(T, found)

    def part_revisions(self) -> dict[str, int]:
        return {key: part.revision for key, part in self._parts.items()}

    def extras(self) -> dict[str, Any]:
        """Document state saved beside the table."""
        out = dict(self._kept)
        for key, part in self._parts.items():
            data = part.capture()
            if data:
                out[key] = data
        if len(self.lists):
            out[LISTS_KEY] = self.lists.to_json()
        return out

    def mark_saved(
        self, checkpoint: Any, lists_revision: int, parts: dict[str, int] | None = None
    ) -> None:
        self.table.history.mark_clean(checkpoint)
        self._saved_lists = lists_revision
        self._saved_parts = dict(parts if parts is not None else self.part_revisions())
        self._unsaved = False

    def close(self) -> None:
        for part in self._parts.values():
            part.close()
        self._parts = {}
        self.model.close()
        self.selection.close()
        self.lists.close()
