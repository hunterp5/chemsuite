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
"""Which atom a ray from the eye hits first."""

from __future__ import annotations

import numpy as np


def pick_sphere(
    origin: np.ndarray, direction: np.ndarray, centres: np.ndarray, radii: np.ndarray
) -> int | None:
    """The index of the nearest sphere the ray enters, or None."""
    if not len(centres):
        return None
    oc = centres - origin
    b = oc @ direction
    h = b * b - np.einsum("ij,ij->i", oc, oc) + radii * radii
    hit = h >= 0
    t = np.where(hit, b - np.sqrt(np.where(hit, h, 0.0)), np.inf)
    t = np.where(t > 0, t, np.inf)
    best = int(np.argmin(t))
    return best if np.isfinite(t[best]) else None
