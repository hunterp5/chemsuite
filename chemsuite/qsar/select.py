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
"""Choosing features out of the pool — on a training fold, never on the whole set.

The pool is of order a few thousand columns on a few hundred rows, so selection is not a
nicety: most columns are noise, and some are the same information twice. The cascade is

1. **flat** columns go — no spread means nothing to learn from;
2. a **univariate filter** cuts the pool to a shortlist, cheaply, by each column's own
   relation to y (F-test, or mutual information when the relation may not be monotonic);
3. **correlated** columns go, keeping whichever of a pair relates to y more strongly — on the
   shortlist, because a full correlation matrix of the raw pool is tens of millions of cells;
4. an **embedded** pass fits one forest and keeps the columns it actually used.

Every function here takes a training fold and returns indices into the pool's columns. Calling
any of them on the whole set is the mistake this module exists to make hard: a score measured
after selecting on all the rows is not a score, which is what the shuffled-y test in
``tests/qsar/test_search.py`` checks.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: The univariate filters. The F-test asks about a straight-line relation, mutual information
#: about any relation at all, at a few times the cost.
FILTERS: dict[str, str] = {"f_test": "F-test", "mutual_info": "Mutual information"}
DEFAULT_FILTER = "f_test"
#: The shortlist the filter cuts to, as a multiple of the features finally kept.
SHORTLIST = 4
#: Trees in the embedded pass: enough for stable importances, few enough to be cheap.
EMBEDDED_TREES = 100


@dataclass(frozen=True)
class Selection:
    """Which columns of the pool a fold kept, and what each stage dropped."""

    kept: np.ndarray  # int64 indices into the pool's columns, in order
    dropped: dict[str, int] = field(default_factory=dict)
    #: the univariate score of every column of the pool (NaN where it was never scored)
    scores: np.ndarray | None = None

    def __len__(self) -> int:
        return len(self.kept)


def flat_columns(x: np.ndarray) -> np.ndarray:
    """The columns with nothing in them: no finite value, or every finite value the same."""
    x = np.asarray(x, dtype=np.float64)
    finite = np.isfinite(x)
    nothing = ~finite.any(axis=0)
    same = np.zeros(x.shape[1], dtype=np.bool_)
    columns = np.flatnonzero(~nothing)  # nanmin warns on a column with no value at all
    if columns.size:
        usable = np.where(finite[:, columns], x[:, columns], np.nan)
        same[columns] = np.nanmin(usable, axis=0) == np.nanmax(usable, axis=0)
    return np.flatnonzero(nothing | same)


def univariate_scores(x: np.ndarray, y: np.ndarray, task: str, kind: str) -> np.ndarray:
    """How strongly each column relates to y on its own. Higher is more interesting; a column
    the filter cannot score at all gets 0, so it ranks last rather than raising."""
    if kind not in FILTERS:
        raise ValueError(f"unknown filter {kind!r}")
    from sklearn.feature_selection import (
        f_classif,
        f_regression,
        mutual_info_classif,
        mutual_info_regression,
    )

    classification = task == "classification"
    if kind == "mutual_info":
        fn = mutual_info_classif if classification else mutual_info_regression
        scores = fn(x, y, random_state=0)
    else:
        fn = f_classif if classification else f_regression
        scores = fn(x, y)[0]
    return np.nan_to_num(np.asarray(scores, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)


def _ranked(columns: np.ndarray, scores: np.ndarray) -> np.ndarray:
    """*columns* best score first, ties by column order so a fit is reproducible."""
    order = np.lexsort((columns, -scores[columns]))
    return columns[order]


def decorrelated(x: np.ndarray, columns: np.ndarray, threshold: float) -> np.ndarray:
    """Walk *columns* in the order given and keep a column only when it correlates with every
    column already kept below *threshold*. The order decides which of a pair survives, so the
    caller passes them best-first."""
    taken = np.asarray(x, dtype=np.float64)[:, columns]
    spread = taken.std(axis=0)
    normal = (taken - taken.mean(axis=0)) / np.where(spread > 0, spread, 1.0)
    rows = max(len(taken), 1)
    kept: list[int] = []
    for i in range(len(columns)):
        if kept:
            against = np.abs(normal[:, kept].T @ normal[:, i]) / rows
            if against.max() >= threshold:
                continue
        kept.append(i)
    return columns[np.asarray(kept, dtype=np.int64)]


def embedded(
    x: np.ndarray, y: np.ndarray, columns: np.ndarray, task: str, keep: int, seed: int
) -> np.ndarray:
    """The columns one forest actually used, best first, cut to *keep*.

    A forest is the cheapest honest judge of a column *in company*: the univariate filter
    cannot see that two columns only matter together, or that a third is redundant.
    """
    if len(columns) <= keep:
        return columns
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

    forest = (RandomForestClassifier if task == "classification" else RandomForestRegressor)(
        n_estimators=EMBEDDED_TREES, random_state=seed, n_jobs=1
    )
    forest.fit(np.asarray(x, dtype=np.float64)[:, columns], y)
    importance = np.asarray(forest.feature_importances_, dtype=np.float64)
    order = np.lexsort((columns, -importance))[:keep]
    return np.sort(columns[order])


def select(
    x: np.ndarray,
    y: np.ndarray,
    *,
    task: str,
    max_features: int = 200,
    correlation: float = 0.95,
    kind: str = DEFAULT_FILTER,
    seed: int = 0,
) -> Selection:
    """The cascade, on a training fold. Returns the pool columns to keep, in order.

    *x* is already filled (:meth:`~.prepare.Prepare.filled`): the filters cannot read a hole,
    and the fill values have to come from this same fold.

    Never called with rows a model will be scored on: *x* and *y* are a training fold's.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y)
    if len(x) != len(y):
        raise ValueError(f"{len(x)} rows of features against {len(y)} values")
    if x.size and not np.isfinite(x).all():
        raise ValueError("selection needs the holes filled first (Prepare.filled)")
    dropped: dict[str, int] = {}
    columns = np.arange(x.shape[1], dtype=np.int64)
    columns = np.setdiff1d(columns, flat_columns(x))
    dropped["flat"] = x.shape[1] - len(columns)
    if not len(columns):
        return Selection(columns, dropped)

    scores = np.full(x.shape[1], np.nan)
    scores[columns] = univariate_scores(x[:, columns], y, task, kind)
    shortlist = _ranked(columns, scores)[: max(max_features * SHORTLIST, max_features)]
    dropped["filter"] = len(columns) - len(shortlist)

    kept = decorrelated(x, shortlist, correlation)
    dropped["correlated"] = len(shortlist) - len(kept)

    final = embedded(x, y, np.sort(kept), task, max_features, seed)
    dropped["unused"] = len(kept) - len(final)
    return Selection(np.sort(final), dropped, scores)
