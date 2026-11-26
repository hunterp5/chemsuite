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
"""Consistency: a library of the residue pairs the pairwise alignments agree on, extended
through every third sequence (Notredame, Higgins & Heringa, J. Mol. Biol. 2000, 302, 205;
MAFFT's L-INS-i, Katoh et al., Nucleic Acids Res. 2005, 33, 511).

For each pair of sequences, a global alignment with free end gaps and a local one each put
residue pairs in the library, weighted by that alignment's identity. Extension adds, for every
third sequence u, the smaller of the two weights that pair s_i with u_k and u_k with t_j: a
pair both sequences reach through u is more likely right. Each residue keeps its :data:`KEEP`
strongest partners in every other sequence.

Layout: sequence s's residues sit end to end at ``starts[s]:starts[s + 1]``, and
``partner[t, starts[s] + i, k]`` is the k-th strongest residue of t for residue i of s (-1 for
none), with its ``weight``.
"""

from __future__ import annotations

import numba
import numpy as np

from .kernels import gotoh

#: the strongest partners kept for each residue in each other sequence
KEEP = 6


@numba.njit(cache=True, nogil=True)
def _pair_paths(a, b, substitution, gap_open, gap_extend):  # pragma: no cover - compiled
    """The free-end-gap and the local alignment of two sequences: who pairs with whom, and
    each alignment's identity."""
    n, m = len(a), len(b)
    score = np.empty((n, m))
    for i in range(n):
        for j in range(m):
            score[i, j] = substitution[a[i], b[j]]
    row_open = np.full(n + 1, gap_open)
    row_extend = np.full(n + 1, gap_extend)
    col_open = np.full(m + 1, gap_open)
    col_extend = np.full(m + 1, gap_extend)
    row_open[0] = row_open[n] = row_extend[0] = row_extend[n] = 0.0
    col_open[0] = col_open[m] = col_extend[0] = col_extend[m] = 0.0
    ones_n, ones_m = np.ones(n), np.ones(m)
    paths = np.full((2, n), -1, np.int32)
    identity = np.zeros(2)
    for source in range(2):
        _, cols_a, cols_b = gotoh(
            score, row_open, row_extend, col_open, col_extend, ones_n, ones_m, source == 1
        )
        same = 0
        aligned = 0
        for c in range(len(cols_a)):
            if cols_a[c] >= 0 and cols_b[c] >= 0:
                paths[source, cols_a[c]] = cols_b[c]
                aligned += 1
                if a[cols_a[c]] == b[cols_b[c]]:
                    same += 1
        identity[source] = same / aligned if aligned else 0.0
    return paths, identity


@numba.njit(cache=True, nogil=True, parallel=True)
def direct(codes, starts, substitution, gap_open, gap_extend):  # pragma: no cover - compiled
    """Every pair's two alignments, both directions.

    ``paths[t, starts[s] + i, source]`` is the residue of t that residue i of s pairs with in
    that source's alignment (0 free-end-gap, 1 local), -1 for none; ``identity[s, t, source]``
    is that alignment's identity.
    """
    count = len(starts) - 1
    total = starts[count]
    paths = np.full((count, total, 2), -1, np.int32)
    identity = np.zeros((count, count, 2))
    for p in numba.prange(count * count):
        s = p // count
        t = p % count
        if t <= s:
            continue
        a = codes[starts[s] : starts[s + 1]]
        b = codes[starts[t] : starts[t + 1]]
        if len(a) == 0 or len(b) == 0:
            continue
        forward, ident = _pair_paths(a, b, substitution, gap_open, gap_extend)
        for source in range(2):
            identity[s, t, source] = ident[source]
            identity[t, s, source] = ident[source]
            for i in range(len(a)):
                j = forward[source, i]
                if j >= 0:
                    paths[t, starts[s] + i, source] = j
                    paths[s, starts[t] + j, source] = i
    return paths, identity


@numba.njit(cache=True, nogil=True, parallel=True)
def extend(paths, identity, starts, keep):  # noqa: C901 - one compiled pass  # pragma: no cover
    """The extended library, ``(partner, weight)``, each (count, total, keep).

    A weight is the direct weight plus every third sequence's path, divided by the number of
    sequences so its scale does not grow with the family.
    """
    count = len(starts) - 1
    total = starts[count]
    partner = np.full((count, total, keep), -1, np.int32)
    weight = np.zeros((count, total, keep), np.float32)
    for p in numba.prange(count * count):
        s = p // count
        t = p % count
        la = starts[s + 1] - starts[s]
        lb = starts[t + 1] - starts[t]
        if s == t or la == 0 or lb == 0:
            continue
        # each residue of s reaches at most 2 + 4 (count - 2) residues of t: note which
        acc = np.zeros((la, lb))
        cap = 2 + 4 * count
        touched = np.empty((la, cap), np.int64)
        touches = np.zeros(la, np.int64)
        for i in range(la):
            for source in range(2):
                j = paths[t, starts[s] + i, source]
                if j >= 0:
                    if acc[i, j] == 0.0:
                        touched[i, touches[i]] = j
                        touches[i] += 1
                    acc[i, j] += identity[s, t, source] + 1e-12
        for u in range(count):
            if u in (s, t):
                continue
            for i in range(la):
                for src1 in range(2):
                    k = paths[u, starts[s] + i, src1]
                    if k < 0:
                        continue
                    w1 = identity[s, u, src1]
                    for src2 in range(2):
                        j = paths[t, starts[u] + k, src2]
                        if j >= 0:
                            if acc[i, j] == 0.0:
                                touched[i, touches[i]] = j
                                touches[i] += 1
                            acc[i, j] += min(w1, identity[u, t, src2]) + 1e-12
        for i in range(la):
            for slot in range(keep):
                best_j = -1
                best_w = 0.0
                for c in range(touches[i]):
                    j = touched[i, c]
                    if acc[i, j] > best_w:
                        best_w = acc[i, j]
                        best_j = j
                if best_j < 0:
                    break
                partner[t, starts[s] + i, slot] = best_j
                weight[t, starts[s] + i, slot] = best_w / count
                acc[i, best_j] = -1.0  # taken (and never counted as untouched again)
    return partner, weight


@numba.njit(cache=True, nogil=True)
def profile_scores(partner, weight, starts, members_a, rows_a, w_a, members_b, rows_b, w_b):
    """(La, Lb) library scores between two clusters' columns.

    For every s in A and t in B, each library pair (s_i, t_j) adds its weight times
    ``w_a[x] * w_b[y]`` (x, y: their places in the clusters) at the columns holding i and j.
    """
    la = rows_a.shape[1]
    lb = rows_b.shape[1]
    out = np.zeros((la, lb))
    for y in range(len(members_b)):
        t = members_b[y]
        column_of = np.full(starts[t + 1] - starts[t], -1, np.int64)
        for c in range(lb):
            if rows_b[y, c] >= 0:
                column_of[rows_b[y, c]] = c
        for x in range(len(members_a)):
            s = members_a[x]
            scale = w_a[x] * w_b[y]
            for ca in range(la):
                i = rows_a[x, ca]
                if i < 0:
                    continue
                base = starts[s] + i
                for slot in range(partner.shape[2]):
                    j = partner[t, base, slot]
                    if j < 0:
                        break
                    out[ca, column_of[j]] += scale * weight[t, base, slot]
    return out


@numba.njit(cache=True, nogil=True)
def sum_of_pairs(partner, weight, starts, rows, weights):  # pragma: no cover - compiled
    """The library's weighted sum-of-pairs score of a full alignment (``rows``: sequences x
    columns of positions): every library pair the alignment puts in one column counts."""
    count, width = rows.shape
    column = np.full(starts[count], -1, np.int64)
    for s in range(count):
        for c in range(width):
            if rows[s, c] >= 0:
                column[starts[s] + rows[s, c]] = c
    total = 0.0
    for s in range(count):
        for t in range(s + 1, count):
            pair = 0.0
            for i in range(starts[s + 1] - starts[s]):
                c = column[starts[s] + i]
                for slot in range(partner.shape[2]):
                    j = partner[t, starts[s] + i, slot]
                    if j < 0:
                        break
                    if column[starts[t] + j] == c:
                        pair += weight[t, starts[s] + i, slot]
            total += weights[s] * weights[t] * pair
    return total
