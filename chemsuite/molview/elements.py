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
"""Element colours and radii.

Colours are the Jmol scheme, the one 3Dmol.js uses by default (``elementColors.Jmol``;
3Dmol.js is BSD-3-Clause, Copyright (c) 2014, University of Pittsburgh and contributors).
Van der Waals radii come from RDKit's periodic table.
"""

from __future__ import annotations

from functools import cache

import numpy as np
from rdkit import Chem

#: atomic number -> 0xRRGGBB (Jmol)
JMOL = {
    1: 0xFFFFFF, 2: 0xD9FFFF, 3: 0xCC80FF, 4: 0xC2FF00, 5: 0xFFB5B5, 6: 0x909090,
    7: 0x3050F8, 8: 0xFF0D0D, 9: 0x90E050, 10: 0xB3E3F5, 11: 0xAB5CF2, 12: 0x8AFF00,
    13: 0xBFA6A6, 14: 0xF0C8A0, 15: 0xFF8000, 16: 0xFFFF30, 17: 0x1FF01F, 18: 0x80D1E3,
    19: 0x8F40D4, 20: 0x3DFF00, 21: 0xE6E6E6, 22: 0xBFC2C7, 23: 0xA6A6AB, 24: 0x8A99C7,
    25: 0x9C7AC7, 26: 0xE06633, 27: 0xF090A0, 28: 0x50D050, 29: 0xC88033, 30: 0x7D80B0,
    31: 0xC28F8F, 32: 0x668F8F, 33: 0xBD80E3, 34: 0xFFA100, 35: 0xA62929, 36: 0x5CB8D1,
    37: 0x702EB0, 38: 0x00FF00, 39: 0x94FFFF, 40: 0x94E0E0, 41: 0x73C2C9, 42: 0x54B5B5,
    43: 0x3B9E9E, 44: 0x248F8F, 45: 0x0A7D8C, 46: 0x006985, 47: 0xC0C0C0, 48: 0xFFD98F,
    49: 0xA67573, 50: 0x668080, 51: 0x9E63B5, 52: 0xD47A00, 53: 0x940094, 54: 0x429EB0,
    55: 0x57178F, 56: 0x00C900, 78: 0xD0D0E0, 79: 0xFFD123, 80: 0xB8B8D0, 81: 0xA6544D,
    82: 0x575961, 83: 0x9E4FB5,
}  # fmt: skip
#: elements not in the table (Jmol's default)
UNKNOWN = 0xFF1493


def rgb(value: int) -> tuple[float, float, float]:
    return ((value >> 16) / 255.0, ((value >> 8) & 255) / 255.0, (value & 255) / 255.0)


@cache
def _colour_table() -> np.ndarray:
    table = np.tile(np.array(rgb(UNKNOWN), np.float32), (119, 1))
    for z, value in JMOL.items():
        table[z] = rgb(value)
    return table


@cache
def _radius_table() -> np.ndarray:
    pt = Chem.GetPeriodicTable()
    return np.array([pt.GetRvdw(z) if z else 1.5 for z in range(119)], np.float32)


def colours(atomic_numbers: np.ndarray) -> np.ndarray:
    """RGB in 0..1 for each atom, shape (n, 3)."""
    return _colour_table()[np.clip(atomic_numbers, 0, 118)]


def vdw_radii(atomic_numbers: np.ndarray) -> np.ndarray:
    return _radius_table()[np.clip(atomic_numbers, 0, 118)]
