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
"""Pairwise sequence alignment: global, local or global with free end gaps.

Sequences are one-letter strings, scored by a substitution matrix (:mod:`.matrices`) and affine
gaps (a gap of k residues costs ``open + k * extend``). The result keeps, for each alignment
column, the residue of either sequence there (-1 for a gap), so it maps straight onto chains.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import matrices

GLOBAL, LOCAL, ENDS_FREE = "global", "local", "ends-free"
MODES = (GLOBAL, LOCAL, ENDS_FREE)


@dataclass(frozen=True, eq=False)
class PairAlignment:
    score: float
    a: np.ndarray  # int64: position in the first sequence at each column, -1 for a gap
    b: np.ndarray  # int64: the same for the second
    identical: int  # aligned pairs with the same letter
    similar: int  # aligned pairs scoring above zero
    aligned: int  # columns with a residue of both

    @property
    def identity(self) -> float:
        """Identical pairs over aligned pairs (0 when nothing is aligned)."""
        return self.identical / self.aligned if self.aligned else 0.0

    @property
    def similarity(self) -> float:
        return self.similar / self.aligned if self.aligned else 0.0

    def strings(self, a: str, b: str) -> tuple[str, str]:
        """The two sequences gapped with ``-``, for display and tests."""
        return _gapped(a, self.a), _gapped(b, self.b)


def align(
    a: str,
    b: str,
    *,
    mode: str = ENDS_FREE,
    matrix: str = matrices.DEFAULT,
    gap_open: float | None = None,
    gap_extend: float | None = None,
) -> PairAlignment:
    """Align two one-letter sequences; gap costs default to BLAST's for the matrix."""
    table = matrices.matrix(matrix)
    default_open, default_extend = matrices.GAPS[matrix]
    gap_open = default_open if gap_open is None else gap_open
    gap_extend = default_extend if gap_extend is None else gap_extend
    pair = table.pair_scores(a, b)
    best, cols_a, cols_b = run_kernel(pair, gap_open, gap_extend, mode)
    both = (cols_a >= 0) & (cols_b >= 0)
    ia, ib = cols_a[both], cols_b[both]
    letters_a = np.frombuffer(a.upper().encode("ascii", "replace"), dtype=np.uint8)
    letters_b = np.frombuffer(b.upper().encode("ascii", "replace"), dtype=np.uint8)
    return PairAlignment(
        score=float(best),
        a=cols_a,
        b=cols_b,
        identical=int((letters_a[ia] == letters_b[ib]).sum()),
        similar=int((pair[ia, ib] > 0).sum()),
        aligned=int(both.sum()),
    )


def run_kernel(
    pair: np.ndarray, gap_open: float, gap_extend: float, mode: str
) -> tuple[float, np.ndarray, np.ndarray]:
    """:func:`.kernels.gotoh` on a matrix of pair scores with the same gap costs everywhere,
    except that free end gaps cost nothing at the ends."""
    n, m = pair.shape
    row_open, row_extend = np.full(n + 1, float(gap_open)), np.full(n + 1, float(gap_extend))
    col_open, col_extend = np.full(m + 1, float(gap_open)), np.full(m + 1, float(gap_extend))
    if mode == ENDS_FREE:
        for costs in (row_open, row_extend, col_open, col_extend):
            costs[0] = costs[-1] = 0.0
    return gotoh(pair, row_open, row_extend, col_open, col_extend, local=mode == LOCAL)


def gotoh(
    pair: np.ndarray,
    row_open: np.ndarray,
    row_extend: np.ndarray,
    col_open: np.ndarray,
    col_extend: np.ndarray,
    *,
    row_occ: np.ndarray | None = None,
    col_occ: np.ndarray | None = None,
    local: bool = False,
) -> tuple[float, np.ndarray, np.ndarray]:
    """The compiled DP (:func:`.kernels.gotoh`), imported only when something is aligned;
    occupancies default to 1 (sequences rather than profiles)."""
    from . import kernels  # numba

    n, m = pair.shape
    row_occ = np.ones(n) if row_occ is None else row_occ
    col_occ = np.ones(m) if col_occ is None else col_occ

    best, cols_a, cols_b = kernels.gotoh(
        np.ascontiguousarray(pair, dtype=np.float64),
        np.ascontiguousarray(row_open, dtype=np.float64),
        np.ascontiguousarray(row_extend, dtype=np.float64),
        np.ascontiguousarray(col_open, dtype=np.float64),
        np.ascontiguousarray(col_extend, dtype=np.float64),
        np.ascontiguousarray(row_occ, dtype=np.float64),
        np.ascontiguousarray(col_occ, dtype=np.float64),
        local,
    )
    return float(best), cols_a, cols_b


def _gapped(sequence: str, columns: np.ndarray) -> str:
    return "".join(sequence[k] if k >= 0 else "-" for k in columns)
