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
"""Undo and redo, for any document part that changes through invertible operations.

Every edit the user can undo is a :class:`Command`: a label plus a list of primitive operations
(:class:`Op`). Each op captures exactly the data it needs to run in both directions and nothing
else. Undo applies the inverses in reverse order; redo replays the operations. The table's ops
are in :mod:`chemsuite.core.table.history`, the assets' in :mod:`chemsuite.core.assets`; a table
and the asset store it holds share one :class:`History`, so a document has one undo stack.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any


class Op(ABC):
    """A primitive, invertible change to a *target* (a table, an asset store): it calls the
    target's private primitives."""

    #: True when the op can change which rows are visible, their order, or the columns.
    structural: bool = True

    @abstractmethod
    def apply(self, target: Any) -> None: ...

    @abstractmethod
    def inverse(self) -> Op: ...

    @property
    def nbytes(self) -> int:
        return 0


@dataclass(eq=False)
class Command:
    label: str
    ops: list[Op]

    @property
    def structural(self) -> bool:
        return any(op.structural for op in self.ops)

    @property
    def nbytes(self) -> int:
        return sum(op.nbytes for op in self.ops)

    def inverted(self) -> Command:
        return Command(self.label, [op.inverse() for op in reversed(self.ops)])


class History:
    """Undo/redo stacks with a count and memory cap (oldest commands are dropped first)."""

    def __init__(self, *, limit: int = 100, max_bytes: int = 512 * 1024**2) -> None:
        self.limit = limit
        self.max_bytes = max_bytes
        self._undo: list[Command] = []
        self._redo: list[Command] = []
        self._clean: Command | None = None
        self._clean_valid = True

    # -- queries -----------------------------------------------------------------------------
    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    @property
    def undo_label(self) -> str | None:
        return self._undo[-1].label if self._undo else None

    @property
    def redo_label(self) -> str | None:
        return self._redo[-1].label if self._redo else None

    @property
    def nbytes(self) -> int:
        return sum(c.nbytes for c in self._undo) + sum(c.nbytes for c in self._redo)

    # -- clean state (document saved) --------------------------------------------------------
    def checkpoint(self) -> object:
        """A token for the current state; pass it to :meth:`mark_clean` once that state is saved."""
        return self._undo[-1] if self._undo else None

    def mark_clean(self, checkpoint: object = ...) -> None:
        """Mark the current state (or an earlier *checkpoint*) as the saved one."""
        self._clean = self.checkpoint() if checkpoint is ... else checkpoint  # type: ignore[assignment]
        self._clean_valid = True

    @property
    def is_clean(self) -> bool:
        top = self._undo[-1] if self._undo else None
        return self._clean_valid and top is self._clean

    # -- mutation (called by the target) -----------------------------------------------------------
    def push(self, command: Command) -> None:
        if self._clean is not None and self._clean in self._redo:
            self._clean_valid = False  # the saved state can no longer be reached by undo/redo
        self._undo.append(command)
        self._redo.clear()
        self._trim()

    def take_undo(self) -> Command:
        command = self._undo.pop()
        self._redo.append(command)
        return command

    def take_redo(self) -> Command:
        command = self._redo.pop()
        self._undo.append(command)
        return command

    def clear(self) -> None:
        self._undo.clear()
        self._redo.clear()
        self._clean = None
        self._clean_valid = True

    def _trim(self) -> None:
        total = sum(c.nbytes for c in self._undo)
        while self._undo and (len(self._undo) > self.limit or total > self.max_bytes):
            dropped = self._undo.pop(0)
            total -= dropped.nbytes
            if self._clean is None or dropped is self._clean:
                self._clean_valid = False  # undo can no longer get back to the saved state
