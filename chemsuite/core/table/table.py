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

"""The compound table: one Qt-free source of truth for rows, columns and the visible view.

Rows
    Each row has a stable integer ``oid`` that never changes when the table is sorted, filtered
    or edited. Rows are stored in ascending-oid order (new rows get larger oids; undo puts deleted
    rows back where they were), so ``oid -> slot`` is a binary search, not a dictionary.

Columns
    Typed, columnar storage (:mod:`.columns`). Column order is display order; no position is
    reserved for any particular column.

View
    The visible rows are an index array computed as *filter → sort* and cached until something
    it depends on changes. Sorting and filtering never reorder the stored data.

Edits
    Every public mutator builds a :class:`~.history.Command` and runs it, so it is undoable, and
    emits one :class:`~.events.AboutToChange` / :class:`~.events.Changed` pair to listeners.
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Iterator, Mapping, Sequence
from contextlib import contextmanager
from types import MappingProxyType

import numpy as np

from . import events as ev
from .assets import Asset
from .columns import (
    NUMBER,
    REACTION,
    STRUCTURE,
    TEXT,
    Column,
    NumberColumn,
    ParsedColumn,
    as_slots,
    blank_column,
    coerce,
    column_from_values,
    to_number,
    to_reaction,
    to_structure,
    to_text,
)
from .filters import Filter
from .history import (
    Command,
    DeleteRows,
    History,
    InsertColumn,
    InsertRows,
    MoveColumn,
    Op,
    PutCells,
    RemoveColumn,
    RenameColumn,
    ReplaceColumn,
    SetStyle,
)
from .style import ColumnStyle

Listener = Callable[[ev.TableEvent], None]
SortKey = tuple[str, bool]  # (column, ascending)

#: Kinds a column is given on purpose, which incoming values of another kind do not replace.
_CHOSEN = (STRUCTURE, REACTION)


class Table:
    def __init__(self, *, history: History | None = None) -> None:
        self._oids = np.empty(0, dtype=np.int64)
        self._columns: dict[str, Column] = {}
        self._next_oid = 1
        self._filters: dict[str, Filter] = {}
        self._sort: tuple[SortKey, ...] = ()
        self._view: np.ndarray | None = None
        self._ordered: tuple[np.ndarray, np.ndarray] | None = None  # (view it goes with, oids)
        #: each filter's mask (see filter_mask); any change but a re-sort or a new filter clears
        #: them, since a filter on every column (the quick search) can read any cell
        self._masks: dict[str, np.ndarray] = {}
        self._listeners: list[Listener] = []
        self._depth = 0
        self._pending: list[ev.Change] = []
        self._pending_structural = False
        self._assets: dict[str, Asset] = {}
        self.history = history or History()

    @classmethod
    def from_columns(
        cls,
        columns: Mapping[str, Column],
        *,
        oids: np.ndarray | None = None,
        next_oid: int | None = None,
    ) -> Table:
        """Build a table directly from finished columns (file loads, documents). Not undoable."""
        lengths = {len(c) for c in columns.values()}
        n = len(oids) if oids is not None else (next(iter(lengths)) if lengths else 0)
        if lengths and lengths != {n}:
            raise ValueError("all columns must have the same length as oids")
        if any(not name for name in columns):
            raise ValueError("column names must not be empty")
        table = cls()
        if oids is None:
            table._oids = np.arange(1, n + 1, dtype=np.int64)
        else:
            table._oids = np.asarray(oids, dtype=np.int64).copy()
            if n > 1 and not (np.diff(table._oids) > 0).all():
                raise ValueError("oids must be strictly increasing")
        table._columns = dict(columns)
        last = int(table._oids[-1]) if n else 0
        table._next_oid = max(last + 1, next_oid or 1)
        return table

    def snapshot(self) -> Table:
        """An independent copy of the data (no history, filters or listeners).

        Cheap: numeric arrays are copied, strings and molecule blobs are immutable and shared.
        Use it to hand a consistent table to a background save or export while the user keeps
        editing (and background parsing keeps filling structures) on the GUI thread.
        """
        everything = np.arange(len(self._oids))
        copy = Table.from_columns(
            {name: column.styled_take(everything) for name, column in self._columns.items()},
            oids=self._oids,
            next_oid=self._next_oid,
        )
        copy._sort = self._sort
        copy._assets = dict(self._assets)  # assets are immutable: sharing them is safe
        return copy

    def rows_copy(self, slots: np.ndarray) -> Table:
        """A new, independent table of the rows at *slots*, in that order, with new oids.

        Columns keep their order and kinds; history, filters, sort and listeners are not
        copied. The assets the rows link to (:data:`~.assets.LINK_COLUMNS`) come too. Used for
        "Copy Shown Rows to New Document".
        """
        from .assets import LINK_COLUMNS

        slots = np.asarray(slots, dtype=np.int64)
        copy = Table.from_columns(
            {name: column.styled_take(slots) for name, column in self._columns.items()}
        )
        named = {
            (kind, str(text))
            for kind, column in LINK_COLUMNS.items()
            if column in self._columns
            for text in self._columns[column].displays(slots)
        }
        copy._assets = {i: a for i, a in self._assets.items() if (a.kind, i) in named}
        return copy

    @property
    def next_oid(self) -> int:
        """The oid the next new row will get (oids are never reused)."""
        return self._next_oid

    # =========================================================================================
    # Reading
    # =========================================================================================
    def __len__(self) -> int:
        return len(self._oids)

    @property
    def oids(self) -> np.ndarray:
        """All oids in storage order (read-only view)."""
        view = self._oids.view()
        view.flags.writeable = False
        return view

    @property
    def column_names(self) -> list[str]:
        return list(self._columns)

    def __contains__(self, name: object) -> bool:
        return name in self._columns

    def column(self, name: str) -> Column:
        """The column object. Treat it as read-only: mutate only through Table methods."""
        try:
            return self._columns[name]
        except KeyError:
            raise KeyError(f"no column named {name!r}") from None

    def kind(self, name: str) -> str:
        return self.column(name).kind

    def slots_of(self, oids: Sequence[int] | np.ndarray) -> np.ndarray:
        """Storage positions of *oids*; raises KeyError if any oid is not in the table."""
        oids = as_slots(oids)
        slots = np.searchsorted(self._oids, oids)
        if len(oids):
            if not len(self._oids):
                raise KeyError(f"oids not in table: {oids[:10].tolist()}")
            clipped = np.minimum(slots, len(self._oids) - 1)
            bad = (slots >= len(self._oids)) | (self._oids[clipped] != oids)
            if bad.any():
                missing = oids[bad][:10].tolist()
                raise KeyError(f"oids not in table: {missing}")
        return slots

    def has_oid(self, oid: int) -> bool:
        slot = int(np.searchsorted(self._oids, oid))
        return slot < len(self._oids) and int(self._oids[slot]) == oid

    def cell(self, name: str, oid: int) -> str:
        """Display text of one cell."""
        return self.column(name).display(int(self.slots_of([oid])[0]))

    def read(
        self, names: Sequence[str], oids: Sequence[int] | np.ndarray | None = None
    ) -> dict[str, np.ndarray]:
        """Typed values for computations (see :meth:`Column.values`), in *oids* order.

        With *oids* omitted, all rows in storage order (aligned with :attr:`oids`).
        """
        slots = None if oids is None else self.slots_of(oids)
        return {name: self.column(name).values(slots) for name in names}

    def unique_name(self, name: str) -> str:
        """*name*, or ``name (2)``, ``name (3)``… if it is taken."""
        if name not in self._columns:
            return name
        n = 2
        while f"{name} ({n})" in self._columns:
            n += 1
        return f"{name} ({n})"

    def unique_names(self, names: Sequence[str]) -> list[str]:
        """*names* with one shared suffix, `` (2)``, `` (3)``…, so that none is taken.

        A tool run that writes several columns gets the same suffix on all of them, so the
        columns of one run read as a set even when only some names were taken before.
        """
        n = 1
        while True:
            candidates = [name if n == 1 else f"{name} ({n})" for name in names]
            if not any(name in self._columns for name in candidates):
                return candidates
            n += 1

    @property
    def nbytes(self) -> int:
        return self._oids.nbytes + sum(c.nbytes for c in self._columns.values())

    # =========================================================================================
    # View: filters and sort
    # =========================================================================================
    @property
    def filters(self) -> dict[str, Filter]:
        return dict(self._filters)

    @property
    def sort_keys(self) -> tuple[SortKey, ...]:
        return self._sort

    def view(self) -> np.ndarray:
        """Storage slots of the visible rows, in display order (read-only)."""
        if self._view is None:
            mask: np.ndarray | None = None
            for key in self._filters:
                keep = self._mask(key)
                mask = keep if mask is None else (mask & keep)
            if mask is None:
                slots = np.arange(len(self._oids), dtype=np.int64)
            else:
                slots = np.flatnonzero(mask).astype(np.int64)
            view = self._sorted(slots)
            view.flags.writeable = False
            self._view = view
        return self._view

    def filter_mask(self, key: str) -> np.ndarray | None:
        """The rows filter *key* keeps (physical order, read-only), or None without it.

        Computed with the view, so asking for it (a match count) costs nothing extra.
        """
        if key not in self._filters:
            return None
        mask = self._mask(key).view()  # read-only, without freezing the filter's own array
        mask.flags.writeable = False
        return mask

    def _mask(self, key: str) -> np.ndarray:
        mask = self._masks.get(key)
        if mask is None:
            mask = self._masks[key] = self._filters[key].matches(self)
        return mask

    def visible_oids(self) -> np.ndarray:
        return self._oids[self.view()]

    def ordered_oids(self) -> np.ndarray:
        """Every row's oid in the sort order, the rows filters hide included (read-only).

        Kept until the view is next recomputed: every change that can reorder rows (sort,
        edits to a sorted column, rows added or removed) recomputes the view.
        """
        view = self.view()
        if self._ordered is None or self._ordered[0] is not view:
            ordered = self._oids[self._sorted(np.arange(len(self._oids), dtype=np.int64))]
            ordered.flags.writeable = False
            self._ordered = (view, ordered)
        return self._ordered[1]

    def set_sort(self, keys: Sequence[SortKey]) -> None:
        """Sort by columns, most significant first. Blanks always sort last. Stable."""
        keys = tuple((name, bool(asc)) for name, asc in keys)
        for name, _ in keys:
            if not self.column(name).sortable:
                raise ValueError(f"column {name!r} cannot be sorted")
        with self._changing(structural=True):
            self._sort = keys
            self._invalidate_view()

    def set_filter(self, key: str, flt: Filter | None) -> None:
        """Add, replace (same *key*) or, with None, remove a filter."""
        if flt is not None:
            for name in flt.columns:
                self.column(name)
        with self._changing(structural=True):
            if flt is None:
                self._filters.pop(key, None)
            else:
                self._filters[key] = flt
            self._masks.pop(key, None)
            self._invalidate_view()

    def clear_filters(self) -> None:
        with self._changing(structural=True):
            self._filters.clear()
            self._masks = {}
            self._invalidate_view()

    def _sorted(self, slots: np.ndarray) -> np.ndarray:
        if not self._sort or len(slots) < 2:
            return slots
        keys: list[np.ndarray] = []  # np.lexsort: last key is the most significant
        for name, ascending in reversed(self._sort):
            column = self._columns[name]
            blank = column.blank_mask()[slots]
            value = column.sort_key()[slots]
            value = np.where(blank, 0.0, value if ascending else -value)
            keys.append(value)
            keys.append(blank)  # blanks last, whatever the direction
        return slots[np.lexsort(keys)]

    def _invalidate_view(self) -> None:
        self._view = None
        self._pending.append(ev.ViewChanged())

    # =========================================================================================
    # Editing (each call is one undoable command)
    # =========================================================================================
    def add_rows(
        self,
        data: Mapping[str, Sequence[object] | Column] | None = None,
        *,
        count: int | None = None,
        label: str = "Add rows",
        undoable: bool = True,
        extra_ops: Sequence[Op] = (),
    ) -> np.ndarray:
        """Append rows. Unknown columns in *data* are created. Returns the new oids.

        Values are inferred as numbers or text; pass a :class:`Column` (e.g. a
        :class:`StructureColumn`) to give the kind explicitly.

        Pass ``undoable=False`` when loading a file, which should not fill the undo history.
        *extra_ops* run first in the same undo step (e.g. adding the protein the rows link to).
        """
        data = dict(data or {})
        lengths = {len(v) for v in data.values()}
        if count is not None:
            lengths.add(count)
        if len(lengths) > 1:
            raise ValueError(f"columns have different lengths: {sorted(lengths)}")
        n = lengths.pop() if lengths else 0
        oids = np.arange(self._next_oid, self._next_oid + n, dtype=np.int64)

        ops: list[Op] = list(extra_ops)
        incoming: dict[str, Column] = {}
        position = len(self._columns)
        for name, values in data.items():
            column = values if isinstance(values, Column) else column_from_values(values)
            if name not in self._columns:
                ops.append(InsertColumn(name, position, blank_column(column.kind, len(self))))
                position += 1
                incoming[name] = column
            else:
                incoming[name], _ = self._fit(name, column, ops)
        kinds = {op.name: op.column.kind for op in ops if isinstance(op, InsertColumn)}
        cells = {
            name: incoming[name]
            if name in incoming
            else blank_column(kinds.get(name) or self._columns[name].kind, n)
            for name in [*self._columns, *kinds]
        }
        ops.append(InsertRows(oids, cells))
        self._execute(Command(label, ops), record=undoable)
        return oids

    def remove_rows(
        self,
        oids: Sequence[int] | np.ndarray,
        *,
        label: str = "Delete rows",
        extra_ops: Sequence[Op] = (),
    ) -> None:
        """Delete rows. *extra_ops* run first in the same undo step (e.g. removing the protein
        the rows link to)."""
        oids = np.unique(as_slots(oids))
        ops: list[Op] = list(extra_ops)
        if len(oids):
            slots = self.slots_of(oids)
            snapshot = {name: col.take(slots) for name, col in self._columns.items()}
            ops.append(DeleteRows(oids, snapshot))
        if ops:
            self._execute(Command(label, ops))

    def duplicate_rows(
        self, oids: Sequence[int] | np.ndarray, *, label: str = "Duplicate rows"
    ) -> np.ndarray:
        """Add a copy of each row in *oids* (in that order) as new rows; returns their oids."""
        slots = self.slots_of(np.asarray(oids, dtype=np.int64))
        if not len(slots):
            return np.empty(0, dtype=np.int64)
        copies = {name: column.take(slots) for name, column in self._columns.items()}
        return self.add_rows(copies, label=label)

    def set_cells(
        self,
        name: str,
        oids: Sequence[int] | np.ndarray,
        values: Sequence[object] | Column,
        *,
        label: str = "Edit cells",
    ) -> None:
        """Write *values* into column *name* at *oids*.

        A number column that receives text becomes a text column (in the same undoable step).
        """
        oids = as_slots(oids)
        if len(oids) != len(values):
            raise ValueError("oids and values must have the same length")
        if not len(oids):
            return
        ops: list[Op] = []
        incoming = values if isinstance(values, Column) else column_from_values(values)
        self._put_ops(name, oids, incoming, ops)
        self._execute(Command(label, ops))

    def set_block(
        self,
        updates: Mapping[str, tuple[Sequence[int] | np.ndarray, Sequence[object] | Column]],
        *,
        label: str = "Edit cells",
    ) -> None:
        """Write cells in several columns as **one** undoable step (paste, find and replace).

        *updates* maps a column name to ``(oids, values)``; each column follows the same rules as
        :meth:`set_cells`.
        """
        ops: list[Op] = []
        for name, (oids, values) in updates.items():
            oids = as_slots(oids)
            if len(oids) != len(values):
                raise ValueError(f"{name!r}: oids and values differ in length")
            if len(oids):
                incoming = values if isinstance(values, Column) else column_from_values(values)
                self._put_ops(name, oids, incoming, ops)
        if ops:
            self._execute(Command(label, ops))

    def fill_structures(self, name: str, oids: np.ndarray, blobs: Sequence[bytes]) -> None:
        """Store background parse results for a structure or reaction column.

        Not an edit (no visible text changes, nothing to undo), but listeners are told the
        cells changed so a view can redraw the molecules and flag invalid ones.
        """
        column = self.column(name)
        if not isinstance(column, ParsedColumn):
            raise TypeError(f"column {name!r} is not a structure or reaction column")
        oids = as_slots(oids)
        with self._changing(structural=False):
            column.fill(self.slots_of(oids), blobs)
            self._pending.append(ev.CellsChanged(name, oids))

    def add_column(
        self,
        name: str,
        values: Sequence[object] | Column | None = None,
        *,
        kind: str | None = None,
        position: int | None = None,
        label: str | None = None,
    ) -> None:
        """Add a column; *values* are in storage order (one per row) or omitted for blanks."""
        if name in self._columns:
            raise ValueError(f"column {name!r} already exists")
        if not name:
            raise ValueError("column name must not be empty")
        if values is None:
            column = blank_column(kind or TEXT, len(self))
        else:
            if len(values) != len(self):
                raise ValueError(f"expected {len(self)} values, got {len(values)}")
            column = values if isinstance(values, Column) else column_from_values(values)
            if kind is not None:
                coerced = coerce(column, kind)
                if coerced is None:
                    raise ValueError(f"values cannot be stored losslessly as {kind}")
                column = coerced
        position = (
            len(self._columns) if position is None else max(0, min(position, len(self._columns)))
        )
        self._execute(
            Command(label or f"Add column {name}", [InsertColumn(name, position, column)])
        )

    def remove_column(self, name: str, *, label: str | None = None) -> None:
        column = self.column(name)
        position = self.column_names.index(name)
        self._execute(
            Command(label or f"Delete column {name}", [RemoveColumn(name, position, column)])
        )

    def remove_columns(self, names: Sequence[str], *, label: str | None = None) -> None:
        """Remove several columns as one undoable step."""
        for name in names:
            self.column(name)  # raises KeyError for unknown names
        order = self.column_names
        positions = sorted({order.index(n) for n in names})
        if not positions:
            return
        # Right to left, so each recorded position is still valid when undone left to right.
        ops: list[Op] = [
            RemoveColumn(order[p], p, self._columns[order[p]]) for p in reversed(positions)
        ]
        default = (
            f"Delete column {order[positions[0]]}"
            if len(positions) == 1
            else f"Delete {len(positions)} columns"
        )
        self._execute(Command(label or default, ops))

    def move_column(self, name: str, position: int) -> None:
        """Move column *name* to *position* (its index after the move)."""
        self.column(name)  # raises KeyError for an unknown name
        old = self.column_names.index(name)
        position = max(0, min(int(position), len(self._columns) - 1))
        if position != old:
            self._execute(Command(f"Move column {name}", [MoveColumn(name, old, position)]))

    def rename_column(self, old: str, new: str) -> None:
        self.column(old)
        if new == old:
            return
        if not new:
            raise ValueError("column name must not be empty")
        if new in self._columns:
            raise ValueError(f"column {new!r} already exists")
        self._execute(Command(f"Rename column {old}", [RenameColumn(old, new)]))

    def convert_column(self, name: str, kind: str, *, force: bool = False) -> int:
        """Change a column's kind. Returns how many cells were blanked (only with *force*)."""
        column = self.column(name)
        if column.kind == kind:
            return 0
        if kind == TEXT:
            new, blanked = to_text(column), 0
        elif kind == NUMBER:
            converted = to_number(column, force=force)
            if converted is None:
                raise ValueError(f"column {name!r} has values that are not plain numbers")
            new, blanked = converted
        elif kind == STRUCTURE:
            new, blanked = to_structure(column), 0
        elif kind == REACTION:
            new, blanked = to_reaction(column), 0
        else:
            raise ValueError(f"unknown column kind {kind!r}")
        new.style = column.style.for_kind(kind)
        self._execute(Command(f"Convert {name} to {kind}", [ReplaceColumn(name, column, new)]))
        return blanked

    # -- how columns are shown ----------------------------------------------------------------
    def set_style(self, name: str, style: ColumnStyle, *, label: str | None = None) -> None:
        """Change column *name*'s precision or colours as one undoable step."""
        old = self.column(name).style
        if style != old:
            self._execute(Command(label or f"Format {name}", [SetStyle(name, old, style)]))

    def set_log10(self, name: str, on: bool) -> int:
        """Convert a number column's values to ``log10(x)`` (*on*) or back to ``10**x``.

        One undoable step that also sets the column's log10 marker. Blank cells stay blank.
        Raises ValueError if *on* and a value is zero or negative. Returns how many cells
        changed.
        """
        column = self.column(name)
        if column.kind != NUMBER:
            raise ValueError(f"column {name!r} is not a number column")
        if column.style.log10 == on:
            return 0
        values = column.values()
        present = ~np.isnan(values)
        if on:
            if (values[present] <= 0).any():
                raise ValueError(f"column {name!r} has values that are zero or negative")
            converted = np.log10(values, where=present, out=np.full(len(values), np.nan))
        else:
            with np.errstate(over="ignore"):
                raw = np.power(10.0, values, where=present, out=np.full(len(values), np.nan))
            # 10**log10(x) is x to about 1e-15; 12 significant digits drops that noise.
            converted = np.array([float(f"{v:.12g}") for v in raw.tolist()])
        new = NumberColumn.from_floats(converted)
        new.style = column.style.with_(log10=on)
        title = "Logarithmic" if on else "Undo logarithmic"
        self._execute(Command(f"{title} {name}", [ReplaceColumn(name, column, new)]))
        return int(present.sum())

    def duplicate_column(self, name: str, new_name: str | None = None) -> str:
        """Copy column *name* (cells and style) to its right as one undoable step."""
        column = self.column(name)
        new_name = new_name or self.unique_name(f"{name} copy")
        if new_name in self._columns:
            raise ValueError(f"column {new_name!r} already exists")
        position = self.column_names.index(name) + 1
        copy = column.styled_take(np.arange(len(column)))
        self._execute(Command(f"Duplicate column {name}", [InsertColumn(new_name, position, copy)]))
        return new_name

    def apply_results(
        self,
        results: Mapping[str, tuple[Sequence[int] | np.ndarray, Sequence[object] | Column]],
        *,
        label: str = "Add results",
        position: int | None = None,
        overwrite: bool | Collection[str] = False,
    ) -> dict[str, str]:
        """Write a tool's output columns as one undoable step.

        *results* maps a requested column name to ``(oids, values)``. An existing column is
        reused when every cell in it is blank, or always with *overwrite* (rerunning a tool
        updates its columns; rows outside *oids* keep their values). *overwrite* may also name
        the only columns to overwrite. A structure column is only
        ever overwritten by structures (regenerated conformers). Otherwise the result goes
        to a new column named ``name (2)``, ``name (3)``…, inserted at *position* (default: the
        end). Returns ``{requested name: actual name}``.
        """
        ops: list[Op] = []
        names: dict[str, str] = {}
        taken: set[str] = set()

        def replace(name: str) -> bool:
            return overwrite if isinstance(overwrite, bool) else name in overwrite

        end = len(self._columns)
        position = end if position is None else max(0, min(position, end))
        for requested, (oids, values) in results.items():
            oids = as_slots(oids)
            if len(oids) != len(values):
                raise ValueError(f"{requested!r}: oids and values differ in length")
            slots = self.slots_of(oids)
            incoming = values if isinstance(values, Column) else column_from_values(values)
            existing = self._columns.get(requested)
            reusable = existing is not None and (
                existing.is_all_blank()
                or (
                    replace(requested)
                    and (existing.kind not in _CHOSEN or incoming.kind == existing.kind)
                )
            )
            if reusable and requested not in taken:
                self._put_ops(requested, oids, incoming, ops)
                names[requested] = requested
            else:
                name = requested
                n = 2
                while name in self._columns or name in taken:
                    name = f"{requested} ({n})"
                    n += 1
                full = blank_column(incoming.kind, len(self))
                full.put(slots, incoming)
                ops.append(InsertColumn(name, position, full))
                position += 1
                names[requested] = name
            taken.add(names[requested])
        if ops:
            self._execute(Command(label, ops))
        return names

    def edit(self, label: str, ops: Sequence[Op]) -> None:
        """Run *ops* (built by a package that extends the table, e.g. asset changes) as one
        undoable step."""
        if ops:
            self._execute(Command(label, list(ops)))

    # -- assets ----------------------------------------------------------------------------------
    @property
    def assets(self) -> Mapping[str, Asset]:
        """The document's assets by id, in the order they were added (read only)."""
        return MappingProxyType(self._assets)

    def set_assets(self, assets: Sequence[Asset]) -> None:
        """Install assets read from a file. Not undoable, like :meth:`from_columns`."""
        self._assets = {asset.id: asset for asset in assets}

    def unique_asset_id(self, wanted: str) -> str:
        """*wanted*, or ``wanted-2``, ``wanted-3``… if an asset already has that id."""
        candidate, n = wanted, 1
        while candidate in self._assets:
            n += 1
            candidate = f"{wanted}-{n}"
        return candidate

    # -- undo / redo -----------------------------------------------------------------------------
    def undo(self) -> str | None:
        """Undo the last command; returns its label (None if there was nothing to undo)."""
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

    # =========================================================================================
    # Listeners
    # =========================================================================================
    def subscribe(self, listener: Listener) -> Callable[[], None]:
        """Register *listener*; returns a function that unregisters it."""
        self._listeners.append(listener)

        def unsubscribe() -> None:
            if listener in self._listeners:
                self._listeners.remove(listener)

        return unsubscribe

    @contextmanager
    def batch(self, *, structural: bool = True) -> Iterator[None]:
        """Group several operations into one AboutToChange/Changed pair."""
        with self._changing(structural):
            yield

    @contextmanager
    def _changing(self, structural: bool) -> Iterator[None]:
        outer = self._depth == 0
        if outer:
            self._pending = []
            self._pending_structural = structural
            self._emit(ev.AboutToChange(structural))
        elif structural and not self._pending_structural:
            raise RuntimeError("structural change inside a non-structural batch")
        self._depth += 1
        try:
            yield
        finally:
            self._depth -= 1
            if outer:
                changes, self._pending = tuple(self._pending), []
                if any(not isinstance(c, ev.ViewChanged) for c in changes):
                    self._masks = {}  # cells, rows, columns or styles changed
                self._emit(ev.Changed(self._pending_structural, changes))

    def _emit(self, event: ev.TableEvent) -> None:
        for listener in list(self._listeners):
            listener(event)

    # =========================================================================================
    # Command execution and primitives (called by history.Op subclasses only)
    # =========================================================================================
    def _execute(self, command: Command, *, record: bool = True) -> None:
        structural = command.structural or any(
            isinstance(op, PutCells) and self._affects_view(op.column) for op in command.ops
        )
        applied: list[Op] = []
        with self._changing(structural):
            try:
                for op in command.ops:
                    op.apply(self)
                    applied.append(op)
            except BaseException:
                for op in reversed(applied):  # leave the table as it was
                    op.inverse().apply(self)
                raise
            # Record before `Changed` goes out, so listeners see the history (undo label,
            # clean state) as of this change, not one step behind.
            if record:
                self.history.push(command)

    def _affects_view(self, name: str) -> bool:
        return any(name == key for key, _ in self._sort) or any(
            name in flt.columns for flt in self._filters.values()
        )

    def _fit(self, name: str, incoming: Column, ops: list[Op]) -> tuple[Column, Column]:
        """Make *incoming* storable in column *name*.

        Returns ``(incoming as the column's kind, the column as it will be)``. When the kinds
        differ, a structure column reads incoming text as structures, an all-blank column takes
        the incoming kind, and a number column receiving text is demoted to text (the last two
        append a ``ReplaceColumn`` op to *ops*).
        """
        target = self._columns[name]
        if target.kind == incoming.kind:
            return incoming, target
        if target.kind == STRUCTURE:  # a structure column is chosen on purpose: keep its kind
            return to_structure(incoming), target
        if target.kind == REACTION:  # so is a reaction column
            return to_reaction(incoming), target
        if target.is_all_blank() and not incoming.is_all_blank():
            replacement = blank_column(incoming.kind, len(target))
            replacement.style = target.style.for_kind(incoming.kind)
            ops.append(ReplaceColumn(name, target, replacement))
            return incoming, replacement
        fitted = coerce(incoming, target.kind)
        if fitted is not None:
            return fitted, target
        replacement = to_text(target)  # a number column receiving text
        replacement.style = target.style.for_kind(TEXT)
        ops.append(ReplaceColumn(name, target, replacement))
        return to_text(incoming), replacement

    def _put_ops(self, name: str, oids: np.ndarray, incoming: Column, ops: list[Op]) -> None:
        slots = self.slots_of(oids)
        fitted, current = self._fit(name, incoming, ops)
        ops.append(PutCells(name, oids, fitted, current.take(slots)))

    def _put_asset(self, asset_id: str, asset: Asset | None) -> None:
        if asset is None:
            self._assets.pop(asset_id, None)
        else:
            self._assets[asset_id] = asset
        self._pending.append(ev.AssetsChanged((asset_id,)))

    def _insert_rows(self, oids: np.ndarray, cells: Mapping[str, Column]) -> None:
        if set(cells) != set(self._columns):
            raise ValueError("inserted rows must provide every column")
        positions = np.searchsorted(self._oids, oids)
        self._oids = np.insert(self._oids, positions, oids)
        for name, column in self._columns.items():
            column.insert(positions, cells[name])
        if len(oids):
            self._next_oid = max(self._next_oid, int(oids.max()) + 1)
        self._view = None
        self._pending.append(ev.RowsInserted(oids))

    def _delete_rows(self, oids: np.ndarray) -> None:
        slots = self.slots_of(oids)
        self._oids = np.delete(self._oids, slots)
        for column in self._columns.values():
            column.delete(slots)
        self._view = None
        self._pending.append(ev.RowsRemoved(oids))

    def _put_cells(self, name: str, oids: np.ndarray, source: Column) -> None:
        self._columns[name].put(self.slots_of(oids), source)
        if self._affects_view(name):
            self._view = None
        self._pending.append(ev.CellsChanged(name, oids))

    def _insert_column(self, name: str, position: int, column: Column) -> None:
        if len(column) != len(self._oids):
            raise ValueError(f"column {name!r} has {len(column)} rows, table has {len(self._oids)}")
        items = list(self._columns.items())
        items.insert(position, (name, column))
        self._columns = dict(items)
        self._pending.append(ev.ColumnInserted(name, position))

    def _remove_column(self, name: str) -> None:
        position = self.column_names.index(name)
        del self._columns[name]
        self._pending.append(ev.ColumnRemoved(name, position))
        dropped_filters = [k for k, f in self._filters.items() if name in f.columns]
        for key in dropped_filters:
            del self._filters[key]
        sort = tuple(k for k in self._sort if k[0] != name)
        if dropped_filters or sort != self._sort:
            self._sort = sort
            self._invalidate_view()

    def _move_column(self, name: str, position: int) -> None:
        items = list(self._columns.items())
        old = [k for k, _ in items].index(name)
        items.insert(position, items.pop(old))
        self._columns = dict(items)
        self._pending.append(ev.ColumnMoved(name, old, position))

    def _rename_column(self, old: str, new: str) -> None:
        self._columns = {(new if k == old else k): v for k, v in self._columns.items()}
        self._filters = {k: f.renamed(old, new) for k, f in self._filters.items()}
        self._sort = tuple(((new if n == old else n), asc) for n, asc in self._sort)
        self._pending.append(ev.ColumnRenamed(old, new))

    def _replace_column(self, name: str, column: Column) -> None:
        if len(column) != len(self._oids):
            raise ValueError("replacement column has the wrong length")
        previous = self._columns[name]
        self._columns[name] = column
        if column.kind != previous.kind:
            # Sort keys and filters were built for the old kind; drop those that no longer fit.
            sort = tuple(k for k in self._sort if k[0] != name or column.sortable)
            stale = [k for k, f in self._filters.items() if name in f.columns]
            for key in stale:
                del self._filters[key]
            if stale or sort != self._sort:
                self._sort = sort
                self._invalidate_view()
        if self._affects_view(name):
            self._view = None
        self._pending.append(ev.ColumnKindChanged(name, column.kind))

    def _set_style(self, name: str, style: ColumnStyle) -> None:
        self._columns[name].style = style
        self._pending.append(ev.ColumnStyleChanged(name))
