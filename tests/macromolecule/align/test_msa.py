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
"""Multiple alignment: every residue placed once and in order, and what it should get right."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.macromolecule.align import msa, tree

AMINO = "ARNDCQEGHILKMFPSTWYV"


def _random_family(rng: np.random.Generator, count: int) -> list[str]:
    """Mutated copies of one ancestor: substitutions, insertions and deletions."""
    ancestor = "".join(rng.choice(list(AMINO), 80))
    out = []
    for _ in range(count):
        seq = []
        for ch in ancestor:
            roll = rng.random()
            if roll < 0.05:
                continue
            seq.append(rng.choice(list(AMINO)) if roll < 0.3 else ch)
            if rng.random() < 0.03:
                seq.append("".join(rng.choice(list(AMINO), rng.integers(1, 4))))
        out.append("".join(seq))
    return out


def test_every_residue_is_placed_once_in_order() -> None:
    rng = np.random.default_rng(5)
    for count in (2, 3, 7):
        sequences = _random_family(rng, count)
        rows = msa.align(sequences)
        assert rows.shape[0] == count
        for seq, row in zip(sequences, rows, strict=True):
            assert list(row[row >= 0]) == list(range(len(seq)))
        assert not (rows < 0).all(axis=0).any()


def test_identical_sequences_need_no_gaps() -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    assert msa.align_strings([seq] * 4) == [seq] * 4


def test_a_shorter_copy_lines_up_with_its_parent() -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    aligned = msa.align_strings([seq, seq[5:25], seq])
    assert aligned[1] == "-" * 5 + seq[5:25] + "-" * (len(seq) - 25)


def test_empty_and_single_inputs() -> None:
    assert msa.align([]).shape == (0, 0)
    assert msa.align_strings(["ACD"]) == ["ACD"]
    assert msa.align_strings(["ACD", ""]) == ["ACD", "---"]


def test_cancel_stops_the_alignment() -> None:
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        msa.align(_random_family(np.random.default_rng(1), 4), cancel=token)


def test_near_copies_share_their_weight() -> None:
    from scipy.cluster.hierarchy import linkage

    # three near-copies and one distinct sequence: the distinct one weighs most
    condensed = np.array([0.01, 0.01, 0.9, 0.01, 0.9, 0.9])
    weights = tree.clustal_weights(linkage(condensed, method="average"), 4)
    assert weights[3] > 2 * weights[0]
    assert weights.sum() == pytest.approx(4.0)
