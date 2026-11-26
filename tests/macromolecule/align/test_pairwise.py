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
"""Pairwise alignment: hand-checked cases, every residue placed once, and gemmi's global score."""

from __future__ import annotations

import gemmi
import numpy as np
import pytest

from chemsuite.macromolecule.align import pairwise as pw
from chemsuite.macromolecule.structure.sequence import one_letter

AMINO = "ARNDCQEGHILKMFPSTWYV"
NAMES = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
         "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")  # fmt: skip
THREE = {one_letter(n): n for n in NAMES}


def test_a_gap_of_k_costs_open_plus_k_extends() -> None:
    result = pw.align("AAAAWWWWWAAAA", "AAAAAAAA", mode=pw.GLOBAL, gap_open=10, gap_extend=1)
    assert result.strings("AAAAWWWWWAAAA", "AAAAAAAA") == ("AAAAWWWWWAAAA", "AAAA-----AAAA")
    assert result.score == 8 * 4 - (10 + 5 * 1)  # BLOSUM62 A/A is 4


def test_free_end_gaps_cost_nothing_and_global_ones_do() -> None:
    a, b = "MKTAYIAKQRQISFVKSHFSRQ", "GGGKTAYIAKQRQLSFVKSHFS"
    ends = pw.align(a, b, mode=pw.ENDS_FREE)
    full = pw.align(a, b, mode=pw.GLOBAL)
    assert ends.strings(a, b) == ("--MKTAYIAKQRQISFVKSHFSRQ", "GGGKTAYIAKQRQLSFVKSHFS--")
    assert full.score == ends.score - 2 * (11 + 2 * 1)


def test_local_keeps_only_the_shared_part_and_lays_the_rest_out_as_gaps() -> None:
    a, b = "MKTAYIAKQRQISFVKSHFSRQ", "GGGKTAYIAKQRQLSFVKSHFS"
    result = pw.align(a, b, mode=pw.LOCAL)
    assert result.strings(a, b) == ("M---KTAYIAKQRQISFVKSHFSRQ", "-GGGKTAYIAKQRQLSFVKSHFS--")
    assert result.aligned == 19 and result.identical == 18


def test_local_with_nothing_in_common_aligns_nothing() -> None:
    result = pw.align("WWWW", "PPPP", mode=pw.LOCAL)
    assert result.score == 0 and result.aligned == 0
    assert result.strings("WWWW", "PPPP") == ("WWWW----", "----PPPP")


@pytest.mark.parametrize("mode", pw.MODES)
def test_every_residue_is_placed_once_in_order(mode: str) -> None:
    rng = np.random.default_rng(3)
    for _ in range(50):
        a = "".join(rng.choice(list(AMINO), rng.integers(0, 40)))
        b = "".join(rng.choice(list(AMINO), rng.integers(0, 40)))
        result = pw.align(a, b, mode=mode)
        assert list(result.a[result.a >= 0]) == list(range(len(a)))
        assert list(result.b[result.b >= 0]) == list(range(len(b)))
        assert not ((result.a < 0) & (result.b < 0)).any()


def test_swapping_the_sequences_keeps_the_score() -> None:
    a, b = "HEAGAWGHEE", "PAWHEAE"
    for mode in pw.MODES:
        assert pw.align(a, b, mode=mode).score == pw.align(b, a, mode=mode).score


def test_the_global_score_is_gemmis() -> None:
    scoring = gemmi.AlignmentScoring("b")  # BLOSUM62, open -10 and extend -1
    rng = np.random.default_rng(1)
    for _ in range(200):
        a = "".join(rng.choice(list(AMINO), rng.integers(1, 60)))
        b = "".join(rng.choice(list(AMINO), rng.integers(1, 60)))
        if rng.random() < 0.5:
            b = a[rng.integers(0, len(a)) :] + b[:10]
        theirs = gemmi.align_string_sequences([THREE[c] for c in a], [THREE[c] for c in b],
                                              [], scoring)  # fmt: skip
        ours = pw.align(a, b, mode=pw.GLOBAL, gap_open=10, gap_extend=1)
        assert ours.score == theirs.score, (a, b)


def test_identity_counts_identical_aligned_pairs() -> None:
    result = pw.align("ACDEFGHIK", "ACDEFGHIK", mode=pw.GLOBAL)
    assert result.identity == 1.0 and result.aligned == 9
