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
"""Iterative refinement of a multiple alignment (tree-dependent restricted partitioning).

Each edge of the guide tree splits the sequences in two; both halves are taken out of the
alignment, their all-gap columns dropped, and realigned as two profiles. The new alignment is
kept when its weighted sum-of-pairs score is higher (Edgar, Nucleic Acids Res. 2004, 32, 1792;
MAFFT's FFT-NS-i and MUSCLE both refine this way). Passes repeat until one changes nothing.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from . import profile as profiles
from . import tree

#: passes over every edge at most
PASSES = 2


def refine(
    rows: np.ndarray,
    codes: Sequence[np.ndarray],
    guide: tree.GuideTree,
    score: Callable[[np.ndarray], float],
    realign: Callable[[list[int], np.ndarray, list[int], np.ndarray], np.ndarray],
    *,
    passes: int = PASSES,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """``rows`` (sequences x columns of positions) improved edge by edge.

    ``score`` is the objective of a full rows matrix; ``realign(members_a, rows_a, members_b,
    rows_b)`` returns the two groups' rows stacked on their new common columns.
    """
    n = rows.shape[0]
    if n < 3:
        return rows
    best = score(rows)
    everyone = np.arange(n)
    for _ in range(passes):
        changed = False
        # every node but the root, deepest first: each is one side of an edge
        for node in range(2 * n - 3, -1, -1):
            check(cancel)
            inside = tree.leaves(guide.merges, n, node)
            outside = list(np.setdiff1d(everyone, inside))
            if not outside:
                continue
            rows_a = _compact(rows[inside])
            rows_b = _compact(rows[outside])
            stacked = realign(inside, rows_a, outside, rows_b)
            trial = np.empty_like(stacked)
            trial[inside + outside] = stacked
            value = score(trial)
            if value > best + 1e-9:
                rows, best, changed = trial, value, True
        if not changed:
            break
    return rows


def _compact(rows: np.ndarray) -> np.ndarray:
    """``rows`` without the columns where every one of them has a gap."""
    return rows[:, (rows >= 0).any(axis=0)]


def letters_of(rows: np.ndarray, codes: Sequence[np.ndarray]) -> np.ndarray:
    """(sequences x columns) of matrix indices, -1 for a gap."""
    out = np.full(rows.shape, -1, dtype=np.int64)
    for s, row in enumerate(rows):
        present = row >= 0
        out[s, present] = codes[s][row[present]]
    return out


def sum_of_pairs(
    rows: np.ndarray, codes: Sequence[np.ndarray], weights: np.ndarray, substitution: np.ndarray,
    gap_open: float, gap_extend: float,
) -> float:  # fmt: skip
    from . import kernels  # numba

    return float(kernels.sum_of_pairs(
        letters_of(rows, codes), np.asarray(weights, dtype=np.float64), substitution,
        gap_open, gap_extend, profiles.TERMINAL, profiles.TERMINAL_EXTEND,
    ))  # fmt: skip
