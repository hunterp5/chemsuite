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
affine gaps (Gotoh, J. Mol. Biol. 1982, 162, 705): a gap of k positions costs
``gap_open + k * gap_extend``. Three modes share it: global (Needleman & Wunsch, J. Mol. Biol.
1970, 48, 443), local (Smith & Waterman, J. Mol. Biol. 1981, 147, 195) and global with free end
gaps, for chains whose termini were not modelled or were cut differently.

Kernels release the GIL (``nogil``) and are cached on disk, so a job thread runs them without
stalling the window and only the first use on a machine pays the compile.
"""

from __future__ import annotations

import numba
import numpy as np

GLOBAL, LOCAL, ENDS_FREE = 0, 1, 2

_NEG = -1.0e18
# the state a cell was reached from, kept per state for the traceback
_FROM_M, _FROM_X, _FROM_Y, _START = 0, 1, 2, 3


@numba.njit(cache=True, nogil=True)
def gotoh(score, gap_open, gap_extend, mode):  # noqa: C901 - one compiled DP: fill, then traceback
    """Best alignment of rows against columns of ``score``.

    Returns ``(best, a, b)``: ``a[k]`` and ``b[k]`` are the row and column at alignment column
    ``k``, -1 for a gap. Every row and column appears once: in local mode the unaligned ends
    are laid out as gap columns before and after the aligned part (the rows' first).
    """
    n, m = score.shape
    o = gap_open + gap_extend  # the first position of a gap
    e = gap_extend
    free_ends = mode == ENDS_FREE
    local = mode == LOCAL
    pm = np.zeros((n + 1, m + 1), np.int8)
    px = np.zeros((n + 1, m + 1), np.int8)
    py = np.zeros((n + 1, m + 1), np.int8)
    m_prev = np.full(m + 1, _NEG)
    x_prev = np.full(m + 1, _NEG)
    y_prev = np.full(m + 1, _NEG)
    m_cur = np.full(m + 1, _NEG)
    x_cur = np.full(m + 1, _NEG)
    y_cur = np.full(m + 1, _NEG)
    m_prev[0] = 0.0
    if not local:
        for j in range(1, m + 1):
            y_prev[j] = 0.0 if free_ends else -(gap_open + j * e)
            py[0, j] = _FROM_M if j == 1 else _FROM_Y
    best = _NEG if not local else 0.0
    end_i, end_j, end_state = (0, 0, 0) if local else (n, m, 0)
    last_col_m = np.full(n + 1, _NEG)  # free ends: the best of each row's last cell
    last_col_x = np.full(n + 1, _NEG)
    last_col_y = np.full(n + 1, _NEG)
    for i in range(1, n + 1):
        m_cur[0] = _NEG
        y_cur[0] = _NEG
        if local:
            x_cur[0] = _NEG
        else:
            x_cur[0] = 0.0 if free_ends else -(gap_open + i * e)
            px[i, 0] = _FROM_M if i == 1 else _FROM_X
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
            # a row's residue against a gap: from above
            v, p = m_prev[j] - o, _FROM_M
            if x_prev[j] - e > v:
                v, p = x_prev[j] - e, _FROM_X
            if y_prev[j] - o > v:
                v, p = y_prev[j] - o, _FROM_Y
            x_cur[j] = v
            px[i, j] = p
            # a column's residue against a gap: from the left
            v, p = m_cur[j - 1] - o, _FROM_M
            if y_cur[j - 1] - e > v:
                v, p = y_cur[j - 1] - e, _FROM_Y
            if x_cur[j - 1] - o > v:
                v, p = x_cur[j - 1] - o, _FROM_X
            y_cur[j] = v
            py[i, j] = p
            if local and m_cur[j] > best:
                best, end_i, end_j, end_state = m_cur[j], i, j, _FROM_M
        last_col_m[i], last_col_x[i], last_col_y[i] = m_cur[m], x_cur[m], y_cur[m]
        m_prev, m_cur = m_cur, m_prev
        x_prev, x_cur = x_cur, x_prev
        y_prev, y_cur = y_cur, y_prev
    if mode == GLOBAL or (free_ends and (n == 0 or m == 0)):
        best, end_state = m_prev[m], _FROM_M
        if x_prev[m] > best:
            best, end_state = x_prev[m], _FROM_X
        if y_prev[m] > best:
            best, end_state = y_prev[m], _FROM_Y
        end_i, end_j = n, m
    elif free_ends:
        best = _NEG
        for j in range(0, m + 1):  # the last row (m_prev and friends after the swap)
            for state, value in ((_FROM_M, m_prev[j]), (_FROM_X, x_prev[j]), (_FROM_Y, y_prev[j])):
                if value > best:
                    best, end_i, end_j, end_state = value, n, j, state
        for i in range(0, n + 1):
            for state, value in ((_FROM_M, last_col_m[i]), (_FROM_X, last_col_x[i]),
                                 (_FROM_Y, last_col_y[i])):  # fmt: skip
                if value > best:
                    best, end_i, end_j, end_state = value, i, m, state
    # traceback, then the ends left out of it
    a_out = np.empty(n + m, np.int64)
    b_out = np.empty(n + m, np.int64)
    k = 0
    for jj in range(m - 1, end_j - 1, -1):
        a_out[k], b_out[k] = -1, jj
        k += 1
    for ii in range(n - 1, end_i - 1, -1):
        a_out[k], b_out[k] = ii, -1
        k += 1
    i, j, state = end_i, end_j, end_state
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
    return best, a_out[:k][::-1].copy(), b_out[:k][::-1].copy()
