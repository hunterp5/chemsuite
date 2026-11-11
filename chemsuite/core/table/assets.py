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
"""Named objects stored with a document beside its table: proteins, later other large data.

An :class:`Asset` is immutable. Changing one replaces it through a :class:`PutAsset` op, so
asset changes undo and redo with the table's own history (one undo stack per document) and
mark the document modified. The core never interprets a payload: it is either the stored bytes
or an object that can turn itself into them (:class:`Packable`), which lets a package decode
an asset once and keep working with the decoded object.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol, runtime_checkable

from .history import Op

if TYPE_CHECKING:
    from .table import Table


@runtime_checkable
class Packable(Protocol):
    """A decoded payload: it can be stored again, and says how much memory it holds."""

    def to_bytes(self) -> bytes: ...

    @property
    def nbytes(self) -> int: ...


@dataclass(frozen=True, eq=False)
class Asset:
    id: str  # stable and unique in its document; table cells refer to it
    kind: str  # what reads it, e.g. "protein"
    name: str  # what the user sees (renaming never changes the id)
    payload: bytes | Packable
    meta: Mapping[str, str] = field(default_factory=dict)  # small facts: title, source

    def data(self) -> bytes:
        payload = self.payload
        return payload if isinstance(payload, bytes) else payload.to_bytes()

    @property
    def nbytes(self) -> int:
        payload = self.payload
        return len(payload) if isinstance(payload, bytes) else payload.nbytes


@dataclass(eq=False)
class PutAsset(Op):
    """Add (old is None), replace, or remove (new is None) one asset."""

    id: str
    old: Asset | None
    new: Asset | None
    structural: bool = field(default=False)

    def apply(self, table: Table) -> None:
        table._put_asset(self.id, self.new)

    def inverse(self) -> Op:
        return PutAsset(self.id, self.new, self.old)

    @property
    def nbytes(self) -> int:
        return self.new.nbytes if self.new is not None else 0
