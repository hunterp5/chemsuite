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
"""Named objects stored with a document: proteins, later other large data.

An :class:`Asset` is immutable. Changing one replaces it through a :class:`PutAsset` op, so
asset changes undo and redo like any other edit. The core never interprets a payload: it is
either the stored bytes or an object that can turn itself into them (:class:`Packable`), which
lets a package decode an asset once and keep working with the decoded object.

An :class:`AssetStore` holds a document's assets. Inside a table document the
:class:`~chemsuite.core.table.Table` holds the store and runs its ops, so adding a protein and
its rows is one undo step and its changes come out with the table's events. Standing alone (a
protein workspace with no table) the store keeps its own history and tells its own listeners.
Either way, what reads assets sees an :class:`AssetHost`.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Protocol, runtime_checkable

from .history import Command, History, Op

#: told which assets were added, replaced or removed
AssetListener = Callable[[tuple[str, ...]], None]


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


@dataclass(frozen=True)
class AssetsChanged:
    """Assets were added, replaced or removed."""

    ids: tuple[str, ...]


@dataclass(eq=False)
class PutAsset(Op):
    """Add (old is None), replace, or remove (new is None) one asset."""

    id: str
    old: Asset | None
    new: Asset | None
    structural: bool = field(default=False)

    def apply(self, target: AssetTarget) -> None:
        target._put_asset(self.id, self.new)

    def inverse(self) -> Op:
        return PutAsset(self.id, self.new, self.old)

    @property
    def nbytes(self) -> int:
        return self.new.nbytes if self.new is not None else 0


class AssetTarget(Protocol):
    """What :class:`PutAsset` runs on: an asset store, or the table holding one."""

    def _put_asset(self, asset_id: str, asset: Asset | None) -> None: ...


class AssetHost(Protocol):
    """What reads and edits a document's assets (the protein workbench) asks of it: a
    :class:`~chemsuite.core.table.Table` in a table document, an :class:`AssetStore` alone."""

    @property
    def assets(self) -> Mapping[str, Asset]: ...

    def unique_asset_id(self, wanted: str) -> str: ...

    def edit(self, label: str, ops: Sequence[Op]) -> None: ...

    def subscribe_assets(self, listener: AssetListener) -> Callable[[], None]: ...


class AssetStore:
    """A document's assets by id, in the order they were added."""

    def __init__(self, *, history: History | None = None) -> None:
        self._assets: dict[str, Asset] = {}
        #: undo and redo when the store stands alone (a table runs its store's ops itself)
        self.history = history or History()
        self._listeners: list[AssetListener] = []
        self._changed: list[str] = []

    # -- reading ---------------------------------------------------------------------------------
    @property
    def assets(self) -> Mapping[str, Asset]:
        """The assets by id (read only)."""
        return MappingProxyType(self._assets)

    def __len__(self) -> int:
        return len(self._assets)

    def unique_asset_id(self, wanted: str) -> str:
        """*wanted*, or ``wanted-2``, ``wanted-3``… if an asset already has that id."""
        candidate, n = wanted, 1
        while candidate in self._assets:
            n += 1
            candidate = f"{wanted}-{n}"
        return candidate

    def copy(self, keep: Callable[[Asset], bool] | None = None) -> AssetStore:
        """A new store (no history or listeners) sharing the assets *keep* accepts (all by
        default): assets are immutable, so sharing them is safe."""
        copy = AssetStore()
        copy._assets = {i: a for i, a in self._assets.items() if keep is None or keep(a)}
        return copy

    # -- changing --------------------------------------------------------------------------------
    def set_assets(self, assets: Iterable[Asset]) -> None:
        """Install assets read from a file. Not undoable, and no listener is told."""
        self._assets = {asset.id: asset for asset in assets}

    def store(self, asset_id: str, asset: Asset | None) -> None:
        """The primitive under :class:`PutAsset`: add, replace or (None) remove one asset. A
        table holding the store calls it from its own ``_put_asset``."""
        if asset is None:
            self._assets.pop(asset_id, None)
        else:
            self._assets[asset_id] = asset

    def _put_asset(self, asset_id: str, asset: Asset | None) -> None:
        self.store(asset_id, asset)
        self._changed.append(asset_id)

    # -- standing alone: edits, undo and listeners -----------------------------------------------
    def edit(self, label: str, ops: Sequence[Op]) -> None:
        """Run *ops* as one undoable step."""
        if ops:
            self._execute(Command(label, list(ops)))

    def undo(self) -> str | None:
        if not self.history.can_undo:
            return None
        command = self.history.take_undo()
        self._execute(command.inverted(), record=False)
        return command.label

    def redo(self) -> str | None:
        if not self.history.can_redo:
            return None
        command = self.history.take_redo()
        self._execute(command, record=False)
        return command.label

    def subscribe_assets(self, listener: AssetListener) -> Callable[[], None]:
        """Register *listener*; returns a function that unregisters it."""
        self._listeners.append(listener)

        def unsubscribe() -> None:
            if listener in self._listeners:
                self._listeners.remove(listener)

        return unsubscribe

    def _execute(self, command: Command, *, record: bool = True) -> None:
        self._changed = []
        applied: list[Op] = []
        try:
            for op in command.ops:
                op.apply(self)
                applied.append(op)
        finally:
            if len(applied) < len(command.ops):  # an op failed: leave the store as it was
                for op in reversed(applied):
                    op.inverse().apply(self)
                self._changed = []
        if record:
            self.history.push(command)
        changed, self._changed = tuple(dict.fromkeys(self._changed)), []
        if changed:
            for listener in list(self._listeners):
                listener(changed)
