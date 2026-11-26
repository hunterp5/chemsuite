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
"""Alignment A0: how fast is an affine-gap DP in numpy, batched numpy and numba?

Times a Gotoh global alignment (score and traceback pointers, the work a real aligner does) on
random 400-residue sequences, one pair at a time and all pairs of 50 sequences, so the A0 review
can say whether the own MSA (docs/ALIGN_PLAN.md) fits its speed gate in numpy alone.

    python scripts/review/align_kernel_bench.py [--length 400] [--count 50]
"""

from __future__ import annotations

import argparse
import itertools
import time

import numpy as np

NEG = -1.0e9


def numpy_pair(a: np.ndarray, b: np.ndarray, sub: np.ndarray, gap_open: float, gap_extend: float):
    """Row-by-row Gotoh; the horizontal gap state as a running maximum along the row."""
    n, m = len(a), len(b)
    j = np.arange(m + 1, dtype=np.float64)
    h = -gap_open - (j - 1) * gap_extend
    h[0] = 0.0
    vert = np.full(m + 1, NEG)
    pointers = np.empty((n + 1, m + 1), dtype=np.int8)
    for i in range(1, n + 1):
        diag = h[:-1] + sub[a[i - 1], b]
        vert = np.maximum(vert - gap_extend, h - gap_open)
        best = np.empty(m + 1)
        best[0] = -gap_open - (i - 1) * gap_extend
        best[1:] = np.maximum(diag, vert[1:])
        # horizontal: max over k < j of best[k] - open - (j - k - 1) * extend
        run = np.maximum.accumulate(best + j * gap_extend)
        horiz = np.full(m + 1, NEG)
        horiz[1:] = run[:-1] - gap_open - (j[1:] - 1) * gap_extend
        h = np.maximum(best, horiz)
        pointers[i] = (h == horiz).astype(np.int8) * 2 + (h == vert)
    return h[-1], pointers


def numpy_batch(seqs_a: np.ndarray, seqs_b: np.ndarray, sub, gap_open, gap_extend):
    """The same recurrence over many equal-length pairs at once (batch on the first axis)."""
    batch, n = seqs_a.shape
    m = seqs_b.shape[1]
    j = np.arange(m + 1, dtype=np.float64)
    h = np.broadcast_to(-gap_open - (j - 1) * gap_extend, (batch, m + 1)).copy()
    h[:, 0] = 0.0
    vert = np.full((batch, m + 1), NEG)
    pointers = np.empty((batch, n + 1, m + 1), dtype=np.int8)
    for i in range(1, n + 1):
        diag = h[:, :-1] + sub[seqs_a[:, i - 1, None], seqs_b]
        vert = np.maximum(vert - gap_extend, h - gap_open)
        best = np.empty((batch, m + 1))
        best[:, 0] = -gap_open - (i - 1) * gap_extend
        best[:, 1:] = np.maximum(diag, vert[:, 1:])
        run = np.maximum.accumulate(best + j * gap_extend, axis=1)
        horiz = np.full((batch, m + 1), NEG)
        horiz[:, 1:] = run[:, :-1] - gap_open - (j[1:] - 1) * gap_extend
        h = np.maximum(best, horiz)
        pointers[:, i] = (h == horiz).astype(np.int8) * 2 + (h == vert)
    return h[:, -1], pointers


def make_numba():
    import numba

    @numba.njit(cache=False, nogil=True)
    def pair(a, b, sub, gap_open, gap_extend):
        n, m = len(a), len(b)
        match = np.empty(m + 1)
        vert = np.empty(m + 1)
        pointers = np.empty((n + 1, m + 1), dtype=np.int8)
        match[0] = 0.0
        for jj in range(1, m + 1):
            match[jj] = -gap_open - (jj - 1) * gap_extend
            vert[jj] = NEG
        for i in range(1, n + 1):
            diag_prev = match[0]
            match[0] = -gap_open - (i - 1) * gap_extend
            horiz = NEG
            ai = a[i - 1]
            for jj in range(1, m + 1):
                v = max(vert[jj] - gap_extend, match[jj] - gap_open)
                horiz = max(horiz - gap_extend, match[jj - 1] - gap_open)
                d = diag_prev + sub[ai, b[jj - 1]]
                diag_prev = match[jj]
                best = d
                p = 0
                if v > best:
                    best = v
                    p = 1
                if horiz > best:
                    best = horiz
                    p = 2
                vert[jj] = v
                match[jj] = best
                pointers[i, jj] = p
        return match[m], pointers

    return pair


def timed(fn, *args, repeat: int = 3) -> float:
    best = float("inf")
    for _ in range(repeat):
        start = time.perf_counter()
        fn(*args)
        best = min(best, time.perf_counter() - start)
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--length", type=int, default=400)
    parser.add_argument("--count", type=int, default=50)
    args = parser.parse_args()
    rng = np.random.default_rng(7)
    sub = rng.integers(-4, 12, size=(21, 21)).astype(np.float64)
    sub = (sub + sub.T) / 2
    seqs = rng.integers(0, 20, size=(args.count, args.length))
    pairs = list(itertools.combinations(range(args.count), 2))
    a, b = seqs[0], seqs[1]

    one = timed(numpy_pair, a, b, sub, 10.0, 1.0)
    print(f"numpy, one pair {args.length}x{args.length}: {one * 1000:.1f} ms")
    print(f"numpy, all {len(pairs)} pairs (estimated): {one * len(pairs):.1f} s")

    ia = np.array([p[0] for p in pairs])
    ib = np.array([p[1] for p in pairs])
    batch = timed(numpy_batch, seqs[ia], seqs[ib], sub, 10.0, 1.0, repeat=1)
    print(f"numpy batched, all {len(pairs)} pairs: {batch:.1f} s")

    pair = make_numba()
    start = time.perf_counter()
    pair(a, b, sub, 10.0, 1.0)
    print(f"numba compile: {time.perf_counter() - start:.2f} s")
    one_nb = timed(pair, a, b, sub, 10.0, 1.0)
    print(f"numba, one pair: {one_nb * 1000:.2f} ms")

    def all_numba():
        for x, y in pairs:
            pair(seqs[x], seqs[y], sub, 10.0, 1.0)

    print(f"numba, all {len(pairs)} pairs: {timed(all_numba, repeat=1):.2f} s")

    s_np, p_np = numpy_pair(a, b, sub, 10.0, 1.0)
    s_nb, _ = pair(a, b, sub, 10.0, 1.0)
    print(f"scores agree: {s_np == s_nb} ({s_np} vs {s_nb})")


if __name__ == "__main__":
    main()
