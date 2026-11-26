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
2. Up to :data:`CONSISTENCY` sequences, a library of the residue pairs the pairwise alignments
   agree on (:mod:`.library`) is built and extended through every third sequence.
3. Up the tree, each pair of clusters is aligned as two profiles (:mod:`.profile`) with the
   shared DP (:mod:`.kernels`), scored by substitution plus the library, and merged.
4. The alignment is refined edge by edge on the same score (:mod:`.refine`).

The result is a matrix of positions, one row per sequence and one column per alignment column,
holding the position in that sequence or -1 for a gap.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from . import matrices, pairwise, refine, tree
from . import profile as profiles

#: what a multiple alignment uses by default, tuned on BAliBASE 3's training half (A0 review)
MATRIX, GAP_OPEN, GAP_EXTEND = "BLOSUM62", 11.0, 1.0
#: the library's weight against the substitution score (0: no consistency)
LIBRARY = 4.0
#: the most sequences consistency is used for; above it, progressive alignment and refinement
CONSISTENCY = 60


def align_strings(sequences: Sequence[str], **options) -> list[str]:
    """``sequences`` gapped with ``-`` so every one has the same length."""
    rows = align(sequences, **options)
    return [
        "".join(seq[k] if k >= 0 else "-" for k in row)
        for seq, row in zip(sequences, rows, strict=True)
    ]


@dataclass(frozen=True, eq=False)
class _Scoring:
    """How columns are scored, shared by the progressive stage, refinement and its objective."""

    codes: Sequence[np.ndarray]
    weights: np.ndarray
    substitution: np.ndarray
    letters: int
    gap_open: float
    gap_extend: float
    library: float = 0.0  # its weight; 0 for none
    starts: np.ndarray | None = None
    partner: np.ndarray | None = None
    pair_weight: np.ndarray | None = None


def align(
    sequences: Sequence[str],
    *,
    matrix: str = MATRIX,
    gap_open: float = GAP_OPEN,
    gap_extend: float = GAP_EXTEND,
    offset: float = 0.0,
    library: float = LIBRARY,
    passes: int = refine.PASSES,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """The (sequences x columns) matrix of positions of a multiple alignment."""
    n = len(sequences)
    table = matrices.matrix(matrix)
    codes = [table.encode(s) for s in sequences]
    if n == 0:
        return np.zeros((0, 0), np.int32)
    if n == 1:
        return np.arange(len(sequences[0]), dtype=np.int32)[None, :]
    joined = np.concatenate(codes)
    starts = np.concatenate(([0], np.cumsum([len(c) for c in codes]))).astype(np.int64)
    use_library = library > 0 and n <= CONSISTENCY
    partner = pair_weight = identity = None
    if use_library:
        from . import library as lib  # numba

        pair_open, pair_extend = matrices.GAPS[matrix]
        check(cancel)
        paths, identity = lib.direct(joined, starts, table.scores, pair_open, pair_extend)
        check(cancel)
        partner, pair_weight = lib.extend(paths, identity, starts, lib.KEEP)

    def identities(first: np.ndarray, second: np.ndarray) -> np.ndarray:
        if identity is not None:
            return identity[first, second, 0]
        return _identities(joined, starts, table.scores, matrix, first, second, cancel)

    guide = tree.guide_tree(sequences, identities)
    scoring = _Scoring(
        codes, guide.weights, table.scores + offset, len(table.letters), gap_open, gap_extend,
        library if use_library else 0.0, starts, partner, pair_weight,
    )  # fmt: skip
    clusters: dict[int, tuple[list[int], np.ndarray]] = {
        s: ([s], np.arange(len(codes[s]), dtype=np.int32)[None, :]) for s in range(n)
    }
    for k, (left, right, _, _) in enumerate(guide.merges):
        check(cancel)
        members_a, rows_a = clusters.pop(int(left))
        members_b, rows_b = clusters.pop(int(right))
        merged = _merge(members_a, rows_a, members_b, rows_b, scoring)
        clusters[n + k] = (members_a + members_b, merged)
    ((members, rows),) = clusters.values()
    out = np.empty_like(rows)
    out[members] = rows
    if passes <= 0:
        return out

    def realign(members_a, rows_a, members_b, rows_b) -> np.ndarray:
        return _merge(members_a, rows_a, members_b, rows_b, scoring)

    def objective(rows: np.ndarray) -> float:
        return _objective(rows, scoring)

    return refine.refine(out, codes, guide, objective, realign, passes=passes, cancel=cancel)


#: pairs handed to the compiled loop at a time, so a cancel is seen within a fraction of a second
_CHUNK = 256


def _identities(
    joined: np.ndarray, starts: np.ndarray, substitution: np.ndarray, matrix: str,
    first: np.ndarray, second: np.ndarray, cancel: CancelToken | None,
) -> np.ndarray:  # fmt: skip
    """Each pair's identity in its free-end-gap alignment, with the matrix's own gap costs."""
    from . import kernels  # numba

    gap_open, gap_extend = matrices.GAPS[matrix]
    out = np.empty(len(first))
    for k in range(0, len(first), _CHUNK):
        check(cancel)
        out[k : k + _CHUNK] = kernels.pair_identities(
            joined, starts, substitution, gap_open, gap_extend,
            first[k : k + _CHUNK].astype(np.int64), second[k : k + _CHUNK].astype(np.int64),
        )  # fmt: skip
    return out


def _group_weights(members: Sequence[int], weights: np.ndarray) -> np.ndarray:
    """A cluster's sequence weights, summing to 1 (as its profile's are)."""
    w = np.asarray([weights[s] for s in members], dtype=np.float64)
    return w / w.sum() if w.sum() > 0 else np.full(len(members), 1.0 / max(len(members), 1))


def _merge(
    members_a: list[int], rows_a: np.ndarray, members_b: list[int], rows_b: np.ndarray,
    scoring: _Scoring,
) -> np.ndarray:  # fmt: skip
    """Align two clusters' profiles; both clusters' rows laid out on the result, stacked."""
    pa = profiles.build(scoring.codes, rows_a, members_a, scoring.weights, scoring.letters)
    pb = profiles.build(scoring.codes, rows_b, members_b, scoring.weights, scoring.letters)
    scores = profiles.column_scores(pa, pb, scoring.substitution)
    if scoring.library > 0:
        from . import library as lib  # numba

        scores += scoring.library * lib.profile_scores(
            scoring.partner, scoring.pair_weight, scoring.starts,
            np.asarray(members_a, np.int64), rows_a, _group_weights(members_a, scoring.weights),
            np.asarray(members_b, np.int64), rows_b, _group_weights(members_b, scoring.weights),
        )  # fmt: skip
    row_open, row_extend = profiles.gap_costs(pa, scoring.gap_open, scoring.gap_extend)
    col_open, col_extend = profiles.gap_costs(pb, scoring.gap_open, scoring.gap_extend)
    _, cols_a, cols_b = pairwise.gotoh(
        scores, row_open, row_extend, col_open, col_extend,
        row_occ=1.0 - pa.gaps, col_occ=1.0 - pb.gaps,
    )  # fmt: skip
    return np.vstack((_spread(rows_a, cols_a), _spread(rows_b, cols_b)))


def _objective(rows: np.ndarray, scoring: _Scoring) -> float:
    """The weighted sum-of-pairs score: substitutions and gaps, plus the library's pairs."""
    value = refine.sum_of_pairs(
        rows, scoring.codes, scoring.weights, scoring.substitution, scoring.gap_open,
        scoring.gap_extend,
    )  # fmt: skip
    if scoring.library > 0:
        from . import library as lib  # numba

        value += scoring.library * lib.sum_of_pairs(
            scoring.partner,
            scoring.pair_weight,
            scoring.starts,
            rows.astype(np.int64),
            np.asarray(scoring.weights, dtype=np.float64),
        )
    return value


def _spread(rows: np.ndarray, columns: np.ndarray) -> np.ndarray:
    """``rows`` laid out on new columns (``columns[k]`` is the old column at k, -1 a gap)."""
    out = np.full((rows.shape[0], len(columns)), -1, dtype=np.int32)
    present = columns >= 0
    out[:, present] = rows[:, columns[present]]
    return out
