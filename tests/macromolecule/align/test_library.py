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
"""The consistency library: pair-HMM posteriors against every path enumerated, and the matrix's
own background and scale."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.macromolecule.align import library, matrices

DELTA, EPSILON = 0.05, 0.6


def _brute_force(a: np.ndarray, b: np.ndarray, odds: np.ndarray) -> np.ndarray:
    """Match posteriors by summing the probability of every path through the pair HMM."""
    n, m = len(a), len(b)
    stay, leave = 1 - 2 * DELTA, 1 - EPSILON
    moves = {"M": {"M": stay, "X": DELTA, "Y": DELTA},
             "X": {"M": leave, "X": EPSILON}, "Y": {"M": leave, "Y": EPSILON}}  # fmt: skip
    post = np.zeros((n, m))
    total = 0.0

    def walk(i: int, j: int, state: str, p: float, matched: list[tuple[int, int]]) -> None:
        nonlocal total
        if i == n and j == m:
            total += p
            for x, y in matched:
                post[x, y] += p
            return
        for nxt, t in moves[state].items():
            if nxt == "M" and i < n and j < m:
                walk(i + 1, j + 1, "M", p * t * odds[a[i], b[j]], [*matched, (i, j)])
            elif nxt == "X" and i < n:
                walk(i + 1, j, "X", p * t, matched)
            elif nxt == "Y" and j < m:
                walk(i, j + 1, "Y", p * t, matched)

    walk(0, 0, "M", 1.0, [])
    return post / total


@pytest.mark.parametrize(("a", "b"), [("WAC", "WC"), ("MKT", "KTA"), ("HEAG", "PAW")])
def test_posteriors_match_every_path_summed(a: str, b: str) -> None:
    table = matrices.matrix("BLOSUM62")
    _, odds = matrices.odds("BLOSUM62")
    codes_a, codes_b = table.encode(a), table.encode(b)
    ours = library._posterior(codes_a, codes_b, np.ascontiguousarray(odds), DELTA, EPSILON)
    assert ours == pytest.approx(_brute_force(codes_a, codes_b, odds), abs=1e-9)


@pytest.mark.parametrize("name", ["BLOSUM45", "BLOSUM62", "BLOSUM80"])
def test_a_blosum_matrix_implies_a_background_that_sums_to_one(name: str) -> None:
    lam, odds = matrices.odds(name)
    core = matrices.matrix(name).scores[:20, :20]
    background = np.linalg.solve(np.exp(lam * core), np.ones(20))
    assert background.sum() == pytest.approx(1.0)
    assert (background > 0).all()
    assert np.allclose(odds[:20, :20], np.exp(lam * core))


def test_blosum62s_scale_is_near_the_published_one() -> None:
    lam, _ = matrices.odds("BLOSUM62")
    assert lam == pytest.approx(0.318, abs=0.01)  # Altschul's ungapped lambda for BLOSUM62
