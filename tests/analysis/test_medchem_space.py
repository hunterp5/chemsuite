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

from __future__ import annotations

import math

import numpy as np

from chemsuite.analysis.medchem_space import (
    BBB,
    GUIDES,
    HIA,
    INSIDE,
    OUTSIDE,
    boiled_egg,
    golden_triangle,
    lipe,
)


def test_boiled_egg_regions() -> None:
    # caffeine (61.8, −1.03): absorbed, not brain-penetrant; aspirin (63.6, 1.31): yolk;
    # sucrose (189.5, −5.4): outside; a blank value: no class.
    tpsa = np.array([61.8, 63.6, 189.5, np.nan])
    wlogp = np.array([-1.03, 1.31, -5.4, 2.0])
    assert boiled_egg(tpsa, wlogp) == [HIA, BBB, OUTSIDE, ""]


def test_golden_triangle_and_lipe() -> None:
    assert golden_triangle([1.5, 1.5, 6.0, math.nan], [300, 460, 300, 300]) == [
        INSIDE,
        OUTSIDE,
        OUTSIDE,
        "",
    ]
    assert lipe([8.0], [3.0]).tolist() == [5.0]


def test_guides_are_closed_shapes_or_lines() -> None:
    for guide in GUIDES.values():
        assert guide.regions or guide.lines
        for region in guide.regions:
            assert len(region.polygon) >= 3
