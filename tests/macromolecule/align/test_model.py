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
"""Saved alignments: round trip, and what is dropped instead of raised."""

from __future__ import annotations

import numpy as np

from chemsuite.macromolecule.align.model import MSA, AlignedRow, Alignment


def test_an_alignment_round_trips() -> None:
    rows = (AlignedRow("1atp", "E", np.array([0, 1, -1, 2], np.int32)),
            AlignedRow("af-p17612", "A", np.array([-1, 0, 1, 2], np.int32)))  # fmt: skip
    saved = Alignment(rows, MSA, {"matrix": "BLOSUM62"}, "PKA").to_json()
    back = Alignment.from_json(saved)
    assert back is not None
    assert back.to_json() == saved
    assert back.width == 4 and back.proteins() == {"1atp", "af-p17612"}


def test_unreadable_alignments_are_dropped() -> None:
    assert Alignment.from_json(None) is None
    assert Alignment.from_json({"rows": []}) is None
    ragged = {"rows": [{"protein": "a", "chain": "A", "positions": [0, 1]},
                       {"protein": "b", "chain": "A", "positions": [0]}]}  # fmt: skip
    assert Alignment.from_json(ragged) is None
    one_bad = {"rows": [{"protein": "a", "chain": "A", "positions": [0, 1]},
                        {"protein": 3, "chain": "A", "positions": [0, 1]}]}  # fmt: skip
    kept = Alignment.from_json(one_bad)
    assert kept is not None and len(kept.rows) == 1
