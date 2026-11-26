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
"""The substitution matrices: pinned, symmetric, NCBI's values, and the letter rule."""

from __future__ import annotations

import hashlib

import numpy as np
import pytest

from chemsuite.macromolecule.align import matrices


def test_the_shipped_matrices_are_the_pinned_ones() -> None:
    assert hashlib.sha256(matrices.DATA.read_bytes()).hexdigest() == matrices.SHA256


@pytest.mark.parametrize("name", matrices.NAMES)
def test_every_matrix_is_symmetric_with_its_gap_costs(name: str) -> None:
    table = matrices.matrix(name)
    assert table.scores.shape == (24, 24)
    assert np.array_equal(table.scores, table.scores.T)
    assert name in matrices.GAPS


def test_blosum62_has_ncbis_values() -> None:
    m = matrices.matrix("BLOSUM62")
    assert m.pair_scores("W", "W")[0, 0] == 11
    assert m.pair_scores("C", "C")[0, 0] == 9
    assert m.pair_scores("A", "R")[0, 0] == -1
    assert m.pair_scores("W", "*")[0, 0] == -4


def test_unusual_letters_follow_the_rule() -> None:
    m = matrices.matrix("BLOSUM62")
    index = {ch: i for i, ch in enumerate(m.letters)}
    assert list(m.encode("uoJ?a")) == [index["C"], index["K"], index["X"], index["X"], index["A"]]
