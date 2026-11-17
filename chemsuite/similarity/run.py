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
"""Running a similarity search over a list of structure cells (a table's rows).

The cheap methods (fingerprints, scaffolds) run first and score every row, whatever order the
methods were chosen in. An expensive one (the common substructure) then scores only the rows a
fingerprint ranks best, as many as ``Search.prefilter`` (0: every row), best first; the rest of
its column stays blank. The ranking comes from the search's own first fingerprint method, or,
when it has none, from a plain ECFP4 Tanimoto computed for the purpose.

Cancelling keeps what is finished, on the rows it belongs to: a cheap method that finished part
way gives its first rows, and an expensive one gives the prefiltered rows it reached, each
score left where its row is.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import Executor

import numpy as np

from chemsuite.analysis.fingerprints import Cell, Fingerprints, compute_fingerprints
from chemsuite.core.jobs import Cancelled, CancelToken

from . import fingerprint, shape3d, substructure
from .fusion import consensus
from .model import (
    FINGERPRINT,
    SUBSTRUCTURE,
    THREE_D,
    Method,
    Reference,
    Scores,
    Search,
    check_references,
)

Progress = Callable[[int, int], None]
#: the progress bar's resolution: each method gets an equal share of it
STEPS = 1000
#: what ranks the rows for an expensive method when no fingerprint method is chosen
RANKING = Method()


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

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` is the :class:`Scores` so
    far, or None if there are none: every row when only an expensive method was cut short (its
    column blank where it did not reach), else the rows the cheap methods all finished.
    """
    check_references(references, search)
    known = dict(fingerprints or {})
    columns: dict[str, np.ndarray] = {}
    nearest: dict[str, np.ndarray] = {}
    order = sorted(range(len(search.methods)), key=lambda i: search.methods[i].expensive)
    for stage, index in enumerate(order):
        method = search.methods[index]
        report = _stage(progress, stage, len(order))
        try:
            values, which = _score(method, cells, references, search, known, executor,
                                   report, cancel)  # fmt: skip
        except Cancelled as stop:
            if stop.partial is not None:
                columns[method.column], nearest[method.column] = stop.partial
            raise Cancelled(_finish(search, columns, nearest, known)
                            if columns else None) from None  # fmt: skip
        columns[method.column] = values
        nearest[method.column] = which
    return _finish(search, columns, nearest, known)


def _finish(
    search: Search,
    columns: dict[str, np.ndarray],
    nearest: dict[str, np.ndarray],
    known: dict[str, Fingerprints],
) -> Scores:
    """The scores so far, in the order the methods were chosen.

    Every column is cut to the rows each of them covers (an expensive method's blank rows
    still count: it was asked to skip them). Each row's nearest reference comes from the
    primary method, or from the next method that found one, so a row the expensive primary
    method skipped still names one. The consensus needs every method's score for the row.
    """
    rows = min((len(v) for v in columns.values()), default=0)
    found = [m for m in search.methods if m.column in columns]
    out = {m.column: columns[m.column][:rows] for m in found}
    which = np.full(rows, -1, dtype=np.int64)
    for method in reversed(found):  # the primary method's wins, else the next one's
        mine = nearest[method.column][:rows]
        which = np.where(mine >= 0, mine, which)
    combined = None
    if len(search.methods) > 1 and len(found) == len(search.methods):
        combined = consensus(search.methods, [out[m.column] for m in search.methods],
                             search.consensus)  # fmt: skip
        scored = np.column_stack([np.isfinite(out[m.column]) for m in search.methods]).all(axis=1)
        combined = np.where(scored, combined, np.nan)  # every method must have scored the row
    return Scores(out, which, combined, known)


def _prefilter(
    cells: Sequence[Cell],
    references: Sequence[Reference],
    search: Search,
    known: dict[str, Fingerprints],
    executor: Executor | None,
    cancel: CancelToken | None,
) -> np.ndarray | None:
    """The rows an expensive method should score, best first (None: every row).

    Ranked by the search's first fingerprint method, else by :data:`RANKING`; rows that could
    not be fingerprinted come last, so they are dropped first.
    """
    if not search.prefilter or search.prefilter >= len(cells):
        return None
    method = next((m for m in search.methods if m.kind == FINGERPRINT), RANKING)
    fps = known.get(method.fingerprint)
    if fps is None:
        fps = compute_fingerprints(cells, method.fingerprint, executor=executor, cancel=cancel)
        known[method.fingerprint] = fps
    ranking, _which = fingerprint.score(fps, references, method, search.fuse, cancel=cancel)
    order = np.argsort(-np.nan_to_num(ranking, nan=-np.inf), kind="stable")
    return order[: search.prefilter]


def _spread(
    found: tuple[np.ndarray, np.ndarray], rows: np.ndarray, total: int
) -> tuple[np.ndarray, np.ndarray]:
    """Scores for the prefiltered *rows* put back where those rows are; the rest stay blank."""
    score, nearest = found
    out = np.full(total, np.nan)
    which = np.full(total, -1, dtype=np.int64)
    out[rows[: len(score)]] = score
    which[rows[: len(nearest)]] = nearest
    return out, which


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
    if method.kind in SUBSTRUCTURE or method.kind in THREE_D:
        rows = None
        if method.expensive:
            try:
                rows = _prefilter(cells, references, search, known, executor, cancel)
            except Cancelled:
                raise Cancelled from None  # nothing of this method's own was scored
        chosen = cells if rows is None else [cells[i] for i in rows.tolist()]
        own = fingerprint.self_rows(references, len(cells))
        mine = own if rows is None else own[rows]
        try:
            if method.kind in THREE_D:
                found = shape3d.score(chosen, references, method, search.fuse, mine,
                                      executor=executor, progress=progress,
                                      cancel=cancel)  # fmt: skip
            else:
                found = substructure.score(chosen, references, method, search.fuse, mine,
                                           executor=executor, progress=progress,
                                           cancel=cancel)  # fmt: skip
        except Cancelled as stop:
            if stop.partial is None or not method.expensive:
                raise  # a cheap method's partial is its first rows, as the cheap methods' are
            # an expensive method's column is blank where it did not reach, whole rows kept
            where = rows if rows is not None else np.arange(len(cells))
            raise Cancelled(_spread(stop.partial, where, len(cells))) from None
        return found if rows is None else _spread(found, rows, len(cells))
    if method.kind != FINGERPRINT:
        raise ValueError(f"{method.column} is not available yet")
    fps = known.get(method.fingerprint)
    if fps is None:
        fps = compute_fingerprints(cells, method.fingerprint, executor=executor, cancel=cancel)
        known[method.fingerprint] = fps
    return fingerprint.score(fps, references, method, search.fuse, progress=progress,
                             cancel=cancel)  # fmt: skip
