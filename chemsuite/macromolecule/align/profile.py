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
"""Profiles: a group of aligned sequences as weighted residue frequencies per column.

Two profiles are aligned on their sum-of-pairs column scores (``F_a S F_b^T``: every pair of
residues across the two columns, weighted) with gap costs that follow the gaps already there,
as ClustalW and MUSCLE do: opening a gap where most sequences already have one costs only what
the rest pay, and gaps at the ends cost a fraction of the rest.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

#: what end gaps cost, as a share of inner ones: opening (MUSCLE halves it) and extending
TERMINAL = 0.5
TERMINAL_EXTEND = 0.5


@dataclass(frozen=True, eq=False)
class Profile:
    freqs: np.ndarray  # (L, letters): weighted residue frequencies, weights summing to 1
    gaps: np.ndarray  # (L,): the weighted share of sequences with a gap in each column

    @property
    def width(self) -> int:
        return len(self.gaps)


def build(
    codes: Sequence[np.ndarray], rows: np.ndarray, members: Sequence[int], weights: np.ndarray,
    letters: int,
) -> Profile:  # fmt: skip
    """The profile of ``rows`` (members x L positions, -1 gaps) of the sequences ``members``."""
    width = rows.shape[1]
    w = np.asarray([weights[s] for s in members], dtype=np.float64)
    w = w / w.sum() if w.sum() > 0 else np.full(len(members), 1.0 / max(len(members), 1))
    freqs = np.zeros((width, letters))
    gaps = np.zeros(width)
    for k, s in enumerate(members):
        row = rows[k]
        present = row >= 0
        np.add.at(freqs, (np.flatnonzero(present), codes[s][row[present]]), w[k])
        gaps[~present] += w[k]
    return Profile(freqs, gaps)


def column_scores(a: Profile, b: Profile, substitution: np.ndarray) -> np.ndarray:
    """(La, Lb) sum-of-pairs scores between every column of ``a`` and of ``b``."""
    return a.freqs @ substitution @ b.freqs.T


def gap_costs(p: Profile, gap_open: float, gap_extend: float) -> tuple[np.ndarray, np.ndarray]:
    """The open and extend costs of a gap at each boundary of ``p`` (L + 1 of each)."""
    width = p.width
    padded = np.concatenate(([0.0], p.gaps, [0.0]))
    # the share of sequences already in a gap beside the boundary only extend theirs
    beside = np.maximum(padded[:-1], padded[1:])
    opening = gap_open * (1.0 - beside)
    extending = np.full(width + 1, float(gap_extend))
    opening[0] *= TERMINAL
    opening[-1] *= TERMINAL
    extending[0] *= TERMINAL_EXTEND
    extending[-1] *= TERMINAL_EXTEND
    return opening, extending
