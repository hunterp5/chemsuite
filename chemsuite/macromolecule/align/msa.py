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
"""Multiple sequence alignment: progressive, on a guide tree, profile against profile.

1. Distances between every pair of sequences (:mod:`.tree`) give a UPGMA guide tree and
   ClustalW's sequence weights.
2. Up the tree, each pair of clusters is aligned as two profiles (:mod:`.profile`) with the
   shared DP (:mod:`.kernels`), and merged.

The result is a matrix of positions, one row per sequence and one column per alignment column,
holding the position in that sequence or -1 for a gap.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from . import matrices, pairwise, refine, tree
from . import profile as profiles

#: what a multiple alignment uses by default, tuned on BAliBASE 3's training half (A0 review)
MATRIX, GAP_OPEN, GAP_EXTEND = "BLOSUM45", 15.0, 1.0


def align_strings(sequences: Sequence[str], **options) -> list[str]:
    """``sequences`` gapped with ``-`` so every one has the same length."""
    rows = align(sequences, **options)
    return [
        "".join(seq[k] if k >= 0 else "-" for k in row)
        for seq, row in zip(sequences, rows, strict=True)
    ]


def align(
    sequences: Sequence[str],
    *,
    matrix: str = MATRIX,
    gap_open: float = GAP_OPEN,
    gap_extend: float = GAP_EXTEND,
    offset: float = 0.0,
    passes: int = refine.PASSES,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """The (sequences x columns) matrix of positions of a progressive alignment."""
    n = len(sequences)
    table = matrices.matrix(matrix)
    codes = [table.encode(s) for s in sequences]
    if n == 0:
        return np.zeros((0, 0), np.int32)
    if n == 1:
        return np.arange(len(sequences[0]), dtype=np.int32)[None, :]

    def identities(first: np.ndarray, second: np.ndarray) -> np.ndarray:
        return _identities(codes, table.scores, matrix, first, second, cancel)

    guide = tree.guide_tree(sequences, identities)
    letters = len(table.letters)
    clusters: dict[int, tuple[list[int], np.ndarray]] = {
        s: ([s], np.arange(len(codes[s]), dtype=np.int32)[None, :]) for s in range(n)
    }
    for k, (left, right, _, _) in enumerate(guide.merges):
        check(cancel)
        members_a, rows_a = clusters.pop(int(left))
        members_b, rows_b = clusters.pop(int(right))
        clusters[n + k] = _merge(members_a, rows_a, members_b, rows_b, codes, guide.weights,
                                 table.scores + offset, letters, gap_open, gap_extend)  # fmt: skip
    ((members, rows),) = clusters.values()
    out = np.empty_like(rows)
    out[members] = rows
    if passes <= 0:
        return out
    substitution = table.scores + offset

    def realign(members_a, rows_a, members_b, rows_b) -> np.ndarray:
        return _merge(members_a, rows_a, members_b, rows_b, codes, guide.weights, substitution,
                      letters, gap_open, gap_extend)[1]  # fmt: skip

    def objective(rows: np.ndarray) -> float:
        return refine.sum_of_pairs(rows, codes, guide.weights, substitution, gap_open, gap_extend)

    return refine.refine(out, codes, guide, objective, realign, passes=passes, cancel=cancel)


#: pairs handed to the compiled loop at a time, so a cancel is seen within a fraction of a second
_CHUNK = 256


def _identities(
    codes: Sequence[np.ndarray], substitution: np.ndarray, matrix: str, first: np.ndarray,
    second: np.ndarray, cancel: CancelToken | None,
) -> np.ndarray:  # fmt: skip
    """Each pair's identity in its free-end-gap alignment, with the matrix's own gap costs."""
    from . import kernels  # numba

    gap_open, gap_extend = matrices.GAPS[matrix]
    joined = np.concatenate(codes) if codes else np.zeros(0, np.int64)
    starts = np.concatenate(([0], np.cumsum([len(c) for c in codes]))).astype(np.int64)
    out = np.empty(len(first))
    for k in range(0, len(first), _CHUNK):
        check(cancel)
        out[k : k + _CHUNK] = kernels.pair_identities(
            joined, starts, substitution, gap_open, gap_extend,
            first[k : k + _CHUNK].astype(np.int64), second[k : k + _CHUNK].astype(np.int64),
        )  # fmt: skip
    return out


def _merge(
    members_a: list[int], rows_a: np.ndarray, members_b: list[int], rows_b: np.ndarray,
    codes: Sequence[np.ndarray], weights: np.ndarray, substitution: np.ndarray, letters: int,
    gap_open: float, gap_extend: float,
) -> tuple[list[int], np.ndarray]:  # fmt: skip
    """Align two clusters' profiles and lay both clusters' rows out on the result."""
    pa = profiles.build(codes, rows_a, members_a, weights, letters)
    pb = profiles.build(codes, rows_b, members_b, weights, letters)
    scores = profiles.column_scores(pa, pb, substitution)
    row_open, row_extend = profiles.gap_costs(pa, gap_open, gap_extend)
    col_open, col_extend = profiles.gap_costs(pb, gap_open, gap_extend)
    _, cols_a, cols_b = pairwise.gotoh(scores, row_open, row_extend, col_open, col_extend,
                                       row_occ=1.0 - pa.gaps, col_occ=1.0 - pb.gaps)  # fmt: skip
    return members_a + members_b, np.vstack((_spread(rows_a, cols_a), _spread(rows_b, cols_b)))


def _spread(rows: np.ndarray, columns: np.ndarray) -> np.ndarray:
    """``rows`` laid out on new columns (``columns[k]`` is the old column at k, -1 a gap)."""
    out = np.full((rows.shape[0], len(columns)), -1, dtype=np.int32)
    present = columns >= 0
    out[:, present] = rows[:, columns[present]]
    return out
