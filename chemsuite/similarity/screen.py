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
"""Screening a library file (a vendor catalog, say) for compounds like the references.

The file is streamed in batches and never held in memory: each batch is scored with the cheap
ranking fingerprint, and only the best rows are kept, at most ``top_k`` of them, in a heap. The
rest are dropped as they go by. When the file is done, the full search (the expensive methods
among them) runs on those survivors alone, and the hits come back as a table, best first.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.table import Table
from chemsuite.core.table.columns import Column, ParsedColumn
from chemsuite.io.load import TableBuilder, read_batches

from . import fingerprint
from .model import (
    CONSENSUS_COLUMN,
    FINGERPRINT,
    NEAREST_COLUMN,
    Reference,
    Scores,
    Search,
    check_references,
)
from .run import RANKING, score_cells

#: rows kept by default; a vendor catalog has millions, and a chemist looks at a few hundred
TOP_K = 1000


@dataclass
class Hits:
    """What a screen found: the hits as a table (best first, with the score columns), how many
    rows the file held, how many of them could be read, and whether the scan was stopped."""

    table: Table
    scanned: int = 0
    readable: int = 0
    cancelled: bool = False  # the scan was stopped, so the file was not read to its end

    def __len__(self) -> int:
        return len(self.table)


def structure_name(batch) -> str | None:
    """The name of the batch's structure column (None: it has none)."""
    for name, values in batch.items():
        if isinstance(values, ParsedColumn):
            return name
    return None


def cells(column: ParsedColumn) -> list[bytes | str | None]:
    """Every row of *column* as the scorers take it (a blob, unparsed text, or None)."""
    return column.worker_inputs(np.arange(len(column)))


def slice_batch(batch, rows: np.ndarray) -> dict:
    """*batch* cut down to *rows* (a batch's values are columns or plain lists)."""
    out: dict[str, object] = {}
    for name, values in batch.items():
        if isinstance(values, Column):
            out[name] = values.take(rows)
        else:
            out[name] = [values[i] for i in rows.tolist()]
    return out


def _lengths(batch) -> int:
    return max((len(values) for values in batch.values()), default=0)


class _Best:
    """The best *limit* rows seen so far, by score (0: keep every row that qualifies).

    Each kept row is held as its own one-row batch, so the rows that lose are dropped with
    their data and a file of any size costs only what is kept.
    """

    def __init__(self, limit: int) -> None:
        self.limit = limit
        self._rows: list[tuple[float, int, dict]] = []  # (score, arrival, batch)
        self._seen = 0

    def __len__(self) -> int:
        return len(self._rows)

    def offer(self, batch, scores: np.ndarray, threshold: float) -> None:
        """Keep the rows of *batch* whose score reaches *threshold* and beats what is held."""
        import heapq

        qualifies = np.flatnonzero(np.nan_to_num(scores, nan=-np.inf) >= threshold)
        for row in qualifies.tolist():
            self._seen += 1
            entry = (float(scores[row]), self._seen, slice_batch(batch, np.array([row])))
            if not self.limit or len(self._rows) < self.limit:
                heapq.heappush(self._rows, entry)
            elif entry[0] > self._rows[0][0]:  # the heap's root is the weakest kept row
                heapq.heapreplace(self._rows, entry)

    def batches(self) -> list[dict]:
        """The kept rows, best first."""
        return [row[2] for row in sorted(self._rows, key=lambda r: (-r[0], r[1]))]


def _scan(
    path: str | Path,
    references: Sequence[Reference],
    search: Search,
    limit: int,
    threshold: float,
    executor: Executor | None,
    progress: Callable[[int, int], None] | None,
    cancel: CancelToken | None,
) -> tuple[_Best, int, int, bool]:
    """Stage one: every row of the file ranked by the cheap fingerprint, the best ones kept.

    Returns ``(the keeper, rows scanned, rows that could be read, whether it was cancelled)``.
    A cancel stops the scan and keeps what it has (the file is not read again).
    """
    method = next((m for m in search.methods if m.kind == FINGERPRINT), RANKING)
    keeper = _Best(limit)
    scanned = readable = 0
    stopped = False
    try:
        for batch in read_batches(path, executor=executor, progress=progress, cancel=cancel):
            name = structure_name(batch)
            if name is None:
                raise ValueError("the file has no structure column")
            scanned += _lengths(batch)
            fps = compute_fingerprints(cells(batch[name]), method.fingerprint)
            readable += int(fps.valid.sum())
            scores, _nearest = fingerprint.score(fps, references, method, search.fuse)
            keeper.offer(batch, scores, threshold)
            if cancel is not None and cancel.cancelled:
                stopped = True
                break
    except Cancelled:  # the reader stops on its own; the rows already kept are still good
        stopped = True
    return keeper, scanned, readable, stopped


def _kept_table(batches: Sequence[dict]) -> Table:
    builder = TableBuilder()
    for batch in batches:
        builder.add(batch)
    return builder.build()


def screen_file(
    path: str | Path,
    references: Sequence[Reference],
    search: Search,
    *,
    top_k: int = TOP_K,
    threshold: float = 0.0,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> Hits:
    """The compounds in *path* that are most like *references*, as a table, best first.

    At most *top_k* rows are kept (0: every row at or above *threshold*), ranked while the file
    streams by the search's own first fingerprint method (or a plain ECFP4 Tanimoto). Every
    method then scores those rows, so the score columns, the consensus and the nearest
    reference are filled as they are for a table search.

    Cancelling stops the scan and scores what was kept by then, so the hits so far come back
    rather than being lost (``Hits.cancelled`` says so); the file is never read twice.
    """
    check_references(references, search)
    keeper, scanned, readable, stopped = _scan(path, references, search, top_k, threshold,
                                               executor, progress, cancel)  # fmt: skip
    table = _kept_table(keeper.batches())
    if not len(table):
        return Hits(table, scanned, readable, stopped)
    name = next(n for n in table.column_names if table.kind(n) == "structure")
    column = table.column(name)
    found = score_cells(cells(column), references, search, executor=executor, progress=progress)
    return Hits(_with_scores(table, found, search, references), scanned, readable, stopped)


def _with_scores(
    table: Table, found: Scores, search: Search, references: Sequence[Reference]
) -> Table:
    """*table* with the search's score columns added, sorted by the first of them, best first.

    The nearest reference is named only when there is more than one to tell apart.
    """
    from chemsuite.core.table.columns import NumberColumn, column_from_values

    everything = np.arange(len(table))
    columns: dict[str, Column] = {
        name: table.column(name).take(everything) for name in table.column_names
    }
    for column in search.columns:
        values = found.columns.get(column)
        if values is not None:
            columns[column] = NumberColumn.from_floats(values, 3)
    if found.consensus is not None:
        columns[CONSENSUS_COLUMN] = NumberColumn.from_floats(found.consensus, 3)
    if len(references) > 1:
        names = [references[i].name if i >= 0 else "" for i in found.nearest.tolist()]
        columns[NEAREST_COLUMN] = column_from_values(names)
    out = Table.from_columns(columns)
    first = next(c for c in search.columns if c in found.columns)
    order = np.argsort(-np.nan_to_num(found.columns[first], nan=-np.inf), kind="stable")
    return out.rows_copy(order)
