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
"""Running a similarity search over a list of structure cells (a table's rows)."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import Executor

import numpy as np

from chemsuite.analysis.fingerprints import Cell, Fingerprints, compute_fingerprints
from chemsuite.core.jobs import Cancelled, CancelToken

from . import fingerprint
from .fusion import consensus
from .model import FINGERPRINT, Method, Reference, Scores, Search

Progress = Callable[[int, int], None]
#: the progress bar's resolution: each method gets an equal share of it
STEPS = 1000


def _stage(progress: Progress | None, index: int, stages: int) -> Progress | None:
    """A progress callback for stage *index* of *stages*, reporting on the whole search."""
    if progress is None:
        return None

    def report(done: int, total: int) -> None:
        share = done / total if total else 1.0
        progress(int(STEPS * (index + share) / stages), STEPS)

    return report


def score_cells(
    cells: Sequence[Cell],
    references: Sequence[Reference],
    search: Search,
    *,
    fingerprints: Mapping[str, Fingerprints] | None = None,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Scores:
    """Every cell scored against *references* by each of *search*'s methods.

    *fingerprints* holds fingerprints of *cells* already at hand (by kind, row for row); the
    others are computed (with *executor* when given) and returned in ``Scores.fingerprints``.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` is the :class:`Scores` of
    the first rows (those every method finished), or None if there are none.
    """
    if not references:
        raise ValueError("add a reference to compare with")
    n = len(cells)
    known = dict(fingerprints or {})
    columns: dict[str, np.ndarray] = {}
    nearest = np.full(n, -1, dtype=np.int64)
    stages = len(search.methods)
    for index, method in enumerate(search.methods):
        report = _stage(progress, index, stages)
        try:
            values, which = _score(method, cells, references, search, known, executor,
                                   report, cancel)  # fmt: skip
        except Cancelled as stop:
            partial = stop.partial
            rows = len(partial[0]) if partial else 0
            if index:  # the methods before this one finished every row
                columns[method.column] = partial[0] if partial else np.empty(0)
            elif partial:
                columns[method.column], nearest[:rows] = partial
            raise Cancelled(_finish(search, columns, nearest, known).head(rows) if rows
                            else None) from None  # fmt: skip
        columns[method.column] = values
        if index == 0:
            nearest = which
    return _finish(search, columns, nearest, known)


def _finish(
    search: Search,
    columns: dict[str, np.ndarray],
    nearest: np.ndarray,
    known: dict[str, Fingerprints],
) -> Scores:
    rows = min(len(v) for v in columns.values()) if columns else 0
    columns = {name: values[:rows] for name, values in columns.items()}
    combined = None
    if len(search.methods) > 1 and len(columns) == len(search.methods):
        combined = consensus(search.methods, [columns[m.column] for m in search.methods],
                             search.consensus)  # fmt: skip
    return Scores(columns, nearest[:rows], combined, known)


def _score(
    method: Method,
    cells: Sequence[Cell],
    references: Sequence[Reference],
    search: Search,
    known: dict[str, Fingerprints],
    executor: Executor | None,
    progress: Progress | None,
    cancel: CancelToken | None,
) -> tuple[np.ndarray, np.ndarray]:
    if method.kind != FINGERPRINT:
        raise ValueError(f"{method.column} is not available yet")
    fps = known.get(method.fingerprint)
    if fps is None:
        fps = compute_fingerprints(cells, method.fingerprint, executor=executor, cancel=cancel)
        known[method.fingerprint] = fps
    return fingerprint.score(fps, references, method, search.fuse, progress=progress,
                             cancel=cancel)  # fmt: skip
