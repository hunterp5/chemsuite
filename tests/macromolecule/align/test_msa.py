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


@pytest.mark.parametrize("library", [0.0, msa.LIBRARY])
def test_every_residue_is_placed_once_in_order(library: float) -> None:
    rng = np.random.default_rng(5)
    for count in (2, 3, 7):
        sequences = _random_family(rng, count)
        rows = msa.align(sequences, library=library)
        assert rows.shape[0] == count
        for seq, row in zip(sequences, rows, strict=True):
            assert list(row[row >= 0]) == list(range(len(seq)))
        assert not (rows < 0).all(axis=0).any()


@pytest.mark.parametrize("library", [0.0, msa.LIBRARY])
def test_identical_sequences_need_no_gaps(library: float) -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    assert msa.align_strings([seq] * 4, library=library) == [seq] * 4


def test_consistency_recovers_a_mutated_family() -> None:
    """Mutated copies of one ancestor, aligned without knowing it: the library keeps more of
    the true columns than substitution alone."""
    rng = np.random.default_rng(9)
    ancestor = "".join(rng.choice(list(AMINO), 120))
    sequences, truth = [], []
    for _ in range(8):
        seq, origin = [], []
        for k, ch in enumerate(ancestor):
            roll = rng.random()
            if roll < 0.08:
                continue
            seq.append(rng.choice(list(AMINO)) if roll < 0.55 else ch)
            origin.append(k)
        sequences.append("".join(seq))
        truth.append(origin)

    def right(rows: np.ndarray) -> float:
        good = total = 0
        for s in range(len(rows)):
            for t in range(s + 1, len(rows)):
                where_t = {int(p): c for c, p in enumerate(rows[t]) if p >= 0}
                origin_t = {k: i for i, k in enumerate(truth[t])}
                for c, p in enumerate(rows[s]):
                    if p < 0 or truth[s][p] not in origin_t:
                        continue
                    total += 1
                    good += where_t.get(origin_t[truth[s][p]]) == c
        return good / total

    assert right(msa.align(sequences, library=msa.LIBRARY)) >= right(
        msa.align(sequences, library=0)
    )


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
