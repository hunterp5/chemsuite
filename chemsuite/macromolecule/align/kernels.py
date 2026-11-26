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
"""The dynamic programming behind every alignment, compiled with numba.

:func:`gotoh` aligns two sequences (or profiles) from their (n, m) matrix of pair scores with
affine gaps (Gotoh, J. Mol. Biol. 1982, 162, 705) whose costs may change along each side: a
gap of k positions at a place costs that place's ``open + k * extend``. Global alignment
(Needleman & Wunsch, J. Mol. Biol. 1970, 48, 443) with zero costs at the ends is global with
free end gaps; local alignment is Smith & Waterman's (J. Mol. Biol. 1981, 147, 195).

A gap "at row boundary i" sits between rows i - 1 and i (0 before the first, n after the last);
it holds columns against nothing. Column boundaries likewise hold rows against nothing. Each
step of a gap is weighted by the occupancy of the row or column it holds (1 for a sequence; the
share of a profile's sequences with a residue there), so under sum-of-pairs a column most
sequences skip costs little to set against a gap, just as it scores little against a residue.

Kernels release the GIL (``nogil``) and are cached on disk, so a job thread runs them without
stalling the window and only the first use on a machine pays the compile.
"""

from __future__ import annotations

import numba
import numpy as np

GLOBAL, LOCAL = 0, 1

_NEG = -1.0e18
# the state a cell was reached from, kept per state for the traceback
_FROM_M, _FROM_X, _FROM_Y, _START = 0, 1, 2, 3


@numba.njit(cache=True, nogil=True)
def gotoh(score, row_open, row_extend, col_open, col_extend, row_occ, col_occ, local):  # noqa: C901 - one compiled DP
    """Best alignment of the rows against the columns of ``score``.

    ``row_open``/``row_extend`` (n + 1) cost a gap at each row boundary, ``col_open``/
    ``col_extend`` (m + 1) one at each column boundary; ``row_occ`` (n) and ``col_occ`` (m)
    weight each gap step by what it holds. Returns ``(best, a, b)``: ``a[k]`` and
    ``b[k]`` are the row and column at alignment column ``k``, -1 for a gap. Every row and
    column appears once: in local mode the unaligned ends are laid out as gap columns before and
    after the aligned part (the rows' first).
    """
    n, m = score.shape
    pm = np.zeros((n + 1, m + 1), np.int8)
    px = np.zeros((n + 1, m + 1), np.int8)  # a row against a gap (moving down)
    py = np.zeros((n + 1, m + 1), np.int8)  # a column against a gap (moving right)
    m_prev = np.full(m + 1, _NEG)
    x_prev = np.full(m + 1, _NEG)
    y_prev = np.full(m + 1, _NEG)
    m_cur = np.full(m + 1, _NEG)
    x_cur = np.full(m + 1, _NEG)
    y_cur = np.full(m + 1, _NEG)
    m_prev[0] = 0.0
    if not local:
        held = 0.0
        for j in range(1, m + 1):
            held += col_occ[j - 1]
            y_prev[j] = -(row_open[0] * col_occ[0] + held * row_extend[0])
            py[0, j] = _FROM_M if j == 1 else _FROM_Y
    best = 0.0
    end_i, end_j = 0, 0
    held_rows = 0.0
    for i in range(1, n + 1):
        held_rows += row_occ[i - 1]
        occ_i = row_occ[i - 1]
        m_cur[0] = _NEG
        y_cur[0] = _NEG
        x_cur[0] = _NEG
        if not local:
            x_cur[0] = -(col_open[0] * row_occ[0] + held_rows * col_extend[0])
            px[i, 0] = _FROM_M if i == 1 else _FROM_X
        r_open = row_open[i] + row_extend[i]
        r_ext = row_extend[i]
        for j in range(1, m + 1):
            # match: from the diagonal
            v, p = m_prev[j - 1], _FROM_M
            if x_prev[j - 1] > v:
                v, p = x_prev[j - 1], _FROM_X
            if y_prev[j - 1] > v:
                v, p = y_prev[j - 1], _FROM_Y
            if local and v < 0.0:
                v, p = 0.0, _START
            m_cur[j] = v + score[i - 1, j - 1]
            pm[i, j] = p
            # row i - 1 against a gap at column boundary j: from above
            c_open = (col_open[j] + col_extend[j]) * occ_i
            c_ext = col_extend[j] * occ_i
            v, p = m_prev[j] - c_open, _FROM_M
            if x_prev[j] - c_ext > v:
                v, p = x_prev[j] - c_ext, _FROM_X
            if y_prev[j] - c_open > v:
                v, p = y_prev[j] - c_open, _FROM_Y
            x_cur[j] = v
            px[i, j] = p
            # column j - 1 against a gap at row boundary i: from the left
            occ_j = col_occ[j - 1]
            v, p = m_cur[j - 1] - r_open * occ_j, _FROM_M
            if y_cur[j - 1] - r_ext * occ_j > v:
                v, p = y_cur[j - 1] - r_ext * occ_j, _FROM_Y
            if x_cur[j - 1] - r_open * occ_j > v:
                v, p = x_cur[j - 1] - r_open * occ_j, _FROM_X
            y_cur[j] = v
            py[i, j] = p
            if local and m_cur[j] > best:
                best, end_i, end_j = m_cur[j], i, j
        m_prev, m_cur = m_cur, m_prev
        x_prev, x_cur = x_cur, x_prev
        y_prev, y_cur = y_cur, y_prev
    state = _FROM_M
    if not local:
        best, end_i, end_j = m_prev[m], n, m
        if x_prev[m] > best:
            best, state = x_prev[m], _FROM_X
        if y_prev[m] > best:
            best, state = y_prev[m], _FROM_Y
    return best, *_traceback(pm, px, py, n, m, end_i, end_j, state)


@numba.njit(cache=True, nogil=True)
def _traceback(pm, px, py, n, m, end_i, end_j, state):  # pragma: no cover - compiled
    """The path back from ``(end_i, end_j)``, with whatever it leaves out as gap columns."""
    a_out = np.empty(n + m, np.int64)
    b_out = np.empty(n + m, np.int64)
    k = 0
    for jj in range(m - 1, end_j - 1, -1):
        a_out[k], b_out[k] = -1, jj
        k += 1
    for ii in range(n - 1, end_i - 1, -1):
        a_out[k], b_out[k] = ii, -1
        k += 1
    i, j = end_i, end_j
    while i > 0 or j > 0:
        if state == _FROM_M:
            if i == 0 or j == 0:
                break
            p = pm[i, j]
            a_out[k], b_out[k] = i - 1, j - 1
            i -= 1
            j -= 1
            k += 1
            if p == _START:
                break
            state = p
        elif state == _FROM_X:
            if i == 0:
                break
            p = px[i, j]
            a_out[k], b_out[k] = i - 1, -1
            i -= 1
            k += 1
            state = p
        else:
            if j == 0:
                break
            p = py[i, j]
            a_out[k], b_out[k] = -1, j - 1
            j -= 1
            k += 1
            state = p
    for jj in range(j - 1, -1, -1):
        a_out[k], b_out[k] = -1, jj
        k += 1
    for ii in range(i - 1, -1, -1):
        a_out[k], b_out[k] = ii, -1
        k += 1
    return a_out[:k][::-1].copy(), b_out[:k][::-1].copy()
