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

"""Open and save compound files: format dispatch, table assembly and structure parsing."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from concurrent.futures import Executor
from pathlib import Path

import numpy as np

from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import parse_structures
from chemsuite.core.table import Table
from chemsuite.core.table.columns import (
    Column,
    StructureColumn,
    blank_column,
    column_from_values,
    concat_columns,
)

from . import delimited, sdf, xlsx
from .formats import CSV, SDF, SMILES, TSV, XLSX, Batch, split_suffix

PARSE_CHUNK = 2_000


class TableBuilder:
    """Assemble a table from batches whose columns may differ.

    Each batch is converted to compact typed columns as it arrives (so a large file never sits
    in memory as Python strings), then all parts of a column are joined once in :meth:`build`.
    Columns that first appear in a later batch are blank in earlier rows.
    """

    def __init__(self) -> None:
        self._parts: dict[str, list[Column]] = {}
        self._rows = 0

    @property
    def rows(self) -> int:
        return self._rows

    def add(self, batch: Batch) -> None:
        lengths = {len(values) for values in batch.values()}
        if len(lengths) > 1:
            raise ValueError("batch columns have different lengths")
        n = lengths.pop() if lengths else 0
        if not n:
            return
        for name, values in batch.items():
            column = values if isinstance(values, Column) else column_from_values(values)
            parts = self._parts.get(name)
            if parts is None:
                parts = self._parts[name] = [blank_column(column.kind, self._rows)]
            parts.append(column)
        for name, parts in self._parts.items():
            if name not in batch:
                parts.append(blank_column(parts[-1].kind, n))
        self._rows += n

    def build(self) -> Table:
        columns = {name: concat_columns(parts) for name, parts in self._parts.items()}
        return Table.from_columns(columns) if columns else Table()


def read_batches(
    path: str | Path,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    structure_column: str | None = delimited.AUTO,
) -> Iterable[Batch]:
    fmt, _ = split_suffix(path)
    if fmt == SDF:
        return sdf.read_sdf(path, executor=executor, progress=progress, cancel=cancel)
    if fmt in (CSV, TSV):
        return delimited.read_delimited(
            path, structure_column=structure_column, progress=progress, cancel=cancel
        )
    if fmt == SMILES:
        return delimited.read_smiles(path, progress=progress, cancel=cancel)
    if fmt == XLSX:
        return xlsx.read_xlsx(path, structure_column=structure_column, cancel=cancel)
    raise AssertionError(fmt)  # split_suffix only returns known formats


def load_table(
    path: str | Path,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    structure_column: str | None = delimited.AUTO,
) -> Table:
    """Read a compound file into a new table (not undoable; the history starts empty).

    Structures from text formats are left unparsed; call :func:`parse_pending_structures`.
    """
    builder = TableBuilder()
    for batch in read_batches(
        path,
        executor=executor,
        progress=progress,
        cancel=cancel,
        structure_column=structure_column,
    ):
        builder.add(batch)
    table = builder.build()
    table.history.mark_clean()
    return table


def parse_pending_structures(
    table: Table,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    columns: Sequence[str] | None = None,
) -> int:
    """Parse every not-yet-parsed structure cell; returns how many cells were invalid.

    Texts are snapshotted per chunk, parsed (in worker processes if *executor* is given) and
    written back with :meth:`Table.fill_structures`, which must run on the thread that owns the
    table. A row deleted or edited meanwhile is skipped: its text no longer matches.
    """
    names = columns or [
        n for n in table.column_names if isinstance(table.column(n), StructureColumn)
    ]
    invalid = 0
    for name in names:
        column = table.column(name)
        assert isinstance(column, StructureColumn)
        pending = np.flatnonzero(column.pending_mask())
        oid_chunks = [
            table.oids[pending[i : i + PARSE_CHUNK]].copy()
            for i in range(0, len(pending), PARSE_CHUNK)
        ]
        text_chunks = [
            column.displays(pending[i : i + PARSE_CHUNK])
            for i in range(0, len(pending), PARSE_CHUNK)
        ]
        done = 0
        for oids, texts, blobs in zip(
            oid_chunks,
            text_chunks,
            ordered_map(parse_structures, text_chunks, executor),
            strict=True,
        ):
            check(cancel)
            keep = [
                i for i, oid in enumerate(oids.tolist()) if _unchanged(table, name, oid, texts[i])
            ]
            if keep:
                table.fill_structures(name, oids[keep], [blobs[i] for i in keep])
                invalid += sum(1 for i in keep if not blobs[i])
            done += len(oids)
            if progress is not None:
                progress(done, len(pending))
    return invalid


def _unchanged(table: Table, name: str, oid: int, text: str) -> bool:
    if name not in table or not table.has_oid(oid):
        return False
    column = table.column(name)
    return isinstance(column, StructureColumn) and table.cell(name, oid) == text


def save_table(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | np.ndarray | None = None,
    columns: Sequence[str] | None = None,
    structure: str | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Export to the format given by *path*'s suffix (written atomically)."""
    fmt, gzipped = split_suffix(path)
    if gzipped:
        raise ValueError("saving compressed files is not supported")
    check(cancel)
    if fmt == SDF:
        sdf.write_sdf(
            table,
            path,
            oids=oids,
            structure=structure,
            columns=columns,
            progress=progress,
            cancel=cancel,
        )
    elif fmt in (CSV, TSV):
        delimited.write_delimited(
            table,
            path,
            delimiter="," if fmt == CSV else "\t",
            oids=oids,
            columns=columns,
            progress=progress,
            cancel=cancel,
        )
    elif fmt == SMILES:
        structure = structure or sdf.first_structure_column(table)
        name_column = next((c for c in ("Name", "Title", "ID") if c in table), None)
        delimited.write_smiles(table, path, structure=structure, name_column=name_column, oids=oids)
    elif fmt == XLSX:
        xlsx.write_xlsx(table, path, oids=oids, columns=columns, cancel=cancel)
