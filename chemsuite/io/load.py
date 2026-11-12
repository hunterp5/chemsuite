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
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import parse_structures
from chemsuite.core.reactions import parse_reactions
from chemsuite.core.table import Table
from chemsuite.core.table.columns import (
    Column,
    ParsedColumn,
    ReactionColumn,
    blank_column,
    column_from_values,
    concat_columns,
)

from . import delimited, rdf, sdf, xlsx
from .formats import CSV, RDF, RXN, SDF, SMILES, TSV, XLSX, Batch, split_suffix

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
    if fmt in (RDF, RXN):
        return rdf.read_rdf(path, executor=executor, progress=progress, cancel=cancel)
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


@dataclass(frozen=True, eq=False)
class ParseChunk:
    """Structure or reaction texts snapshotted from the table, with their parse results once
    computed."""

    column: str
    oids: np.ndarray
    texts: list[str]
    blobs: list[bytes] | None = None
    reactions: bool = False


def pending_structure_chunks(
    table: Table, columns: Sequence[str] | None = None, *, chunk: int = PARSE_CHUNK
) -> list[ParseChunk]:
    """Snapshot every not-yet-parsed structure and reaction cell (call on the thread that owns
    the table)."""
    names = columns or [n for n in table.column_names if isinstance(table.column(n), ParsedColumn)]
    out: list[ParseChunk] = []
    for name in names:
        column = table.column(name)
        if not isinstance(column, ParsedColumn):
            raise TypeError(f"column {name!r} is not a structure or reaction column")
        reactions = isinstance(column, ReactionColumn)
        pending = np.flatnonzero(column.pending_mask())
        for i in range(0, len(pending), chunk):
            slots = pending[i : i + chunk]
            out.append(
                ParseChunk(name, table.oids[slots].copy(), column.displays(slots), None, reactions)
            )
    return out


def parse_chunks(
    chunks: Iterable[ParseChunk],
    *,
    executor: Executor | None = None,
    cancel: CancelToken | None = None,
) -> Iterable[ParseChunk]:
    """Parse snapshotted chunks (safe on any thread; in worker processes with *executor*)."""
    chunks = list(chunks)
    jobs = ((c.reactions, c.texts) for c in chunks)
    for parsed, blobs in zip(
        chunks, ordered_map(_parse_job, jobs, executor, cancel=cancel), strict=True
    ):
        check(cancel)
        yield ParseChunk(parsed.column, parsed.oids, parsed.texts, blobs, parsed.reactions)


def _parse_job(job: tuple[bool, list[str]]) -> list[bytes]:
    reactions, texts = job
    return parse_reactions(texts) if reactions else parse_structures(texts)


def apply_parsed(table: Table, chunk: ParseChunk) -> int:
    """Store one parsed chunk; returns how many of its cells were invalid.

    Rows deleted or edited since the snapshot are skipped: their text no longer matches, so the
    result would describe a structure that is not there any more.
    """
    if chunk.blobs is None:
        raise ValueError("chunk has not been parsed")
    column = table.column(chunk.column) if chunk.column in table else None
    if not isinstance(column, ParsedColumn):
        return 0
    if isinstance(column, ReactionColumn) != chunk.reactions:
        return 0  # converted to another kind since the snapshot
    keep = [
        i
        for i, oid in enumerate(chunk.oids.tolist())
        if table.has_oid(oid) and table.cell(chunk.column, oid) == chunk.texts[i]
    ]
    if not keep:
        return 0
    table.fill_structures(chunk.column, chunk.oids[keep], [chunk.blobs[i] for i in keep])
    return sum(1 for i in keep if not chunk.blobs[i])


def parse_pending_structures(
    table: Table,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    columns: Sequence[str] | None = None,
) -> int:
    """Snapshot, parse and apply in one call (for scripts and tests); returns invalid count.

    The GUI runs :func:`parse_chunks` on a worker thread and :func:`apply_parsed` on the GUI
    thread instead, so the window stays responsive while a large file is parsed.
    """
    chunks = pending_structure_chunks(table, columns)
    total = sum(len(c.oids) for c in chunks)
    invalid = done = 0
    for parsed in parse_chunks(chunks, executor=executor, cancel=cancel):
        invalid += apply_parsed(table, parsed)
        done += len(parsed.oids)
        if progress is not None:
            progress(done, total)
    return invalid


def save_table(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | np.ndarray | None = None,
    columns: Sequence[str] | None = None,
    structure: str | None = None,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Export to the format given by *path*'s suffix (written atomically). *executor* lays out
    reactions for RD files in worker processes."""
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
    elif fmt == RDF:
        rdf.write_rdf(
            table,
            path,
            oids=oids,
            columns=columns,
            executor=executor,
            progress=progress,
            cancel=cancel,
        )
    elif fmt == RXN:
        rows = table.oids if oids is None else np.asarray(oids)
        if len(rows) != 1:
            raise ValueError("an RXN file holds one reaction: save several as an RD file")
        rdf.write_rxn(table, path, int(rows[0]))
