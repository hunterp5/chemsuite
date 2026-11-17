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
"""Several methods' scores combined into one consensus score per row.

Each method's score is first put on a 0–1 scale where 1 is identical (shape + features adds
two Tanimotos, so it is halved; the field similarity can be negative, so it is clipped at 0).
Then, row by row, over the methods that scored the row:

- **mean**: the mean of those scores;
- **rank**: the mean of each method's percentile rank (0 the least similar row, 1 the most;
  ties share their mean rank), which ignores how differently the methods spread their scores;
- **normalized**: the mean after stretching each method's scores to fill 0–1 (a method that
  gives every row the same score says nothing and is left out).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .model import COLOR, FIELD, Method


def unit(method: Method, values: np.ndarray) -> np.ndarray:
    """*values* (that method's scores) on a 0–1 scale; NaN stays NaN."""
    values = np.asarray(values, dtype=np.float64)
    if method.kind == COLOR:
        values = values / 2.0
    elif method.kind == FIELD:
        values = np.clip(values, 0.0, 1.0)
    return values


def percentile_ranks(values: np.ndarray) -> np.ndarray:
    """Each finite value's rank among the finite values, 0 (lowest) to 1 (highest), ties
    averaged; NaN stays NaN. A single finite value ranks 1."""
    out = np.full(len(values), np.nan)
    finite = np.flatnonzero(np.isfinite(values))
    if not len(finite):
        return out
    unique, inverse, counts = np.unique(values[finite], return_inverse=True, return_counts=True)
    mean_rank = np.cumsum(counts) - (counts - 1) / 2.0  # 1-based, ties averaged
    ranks = mean_rank[inverse]
    out[finite] = (ranks - 1) / (len(finite) - 1) if len(finite) > 1 else 1.0
    return out


def _stretched(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    if not finite.any():
        return values
    low, high = values[finite].min(), values[finite].max()
    if high <= low:
        return np.full(len(values), np.nan)
    return (values - low) / (high - low)


def consensus(methods: Sequence[Method], scores: Sequence[np.ndarray], how: str) -> np.ndarray:
    """One score per row from each method's *scores* (row for row), combined as *how*
    (``mean``, ``rank`` or ``normalized``); NaN for a row no method scored."""
    columns = [unit(m, s) for m, s in zip(methods, scores, strict=True)]
    if how == "rank":
        columns = [percentile_ranks(c) for c in columns]
    elif how == "normalized":
        columns = [_stretched(c) for c in columns]
    elif how != "mean":
        raise ValueError(f"unknown consensus: {how!r}")
    stack = np.column_stack(columns)
    counted = np.isfinite(stack).sum(axis=1)
    total = np.where(np.isfinite(stack), stack, 0.0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counted > 0, total / np.maximum(counted, 1), np.nan)
