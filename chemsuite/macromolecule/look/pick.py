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
"""Which atom is under the pointer, fast enough to ask on every mouse move of a big protein.

The spheres' centres are projected to the screen once per camera position and put in a 2D
k-d tree. A pointer then only ray-tests the few spheres whose screen discs could cover it
(:func:`chemsuite.molview.pick.pick_sphere`, exact), instead of every atom. No GPU readback.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.molview.camera import Camera
from chemsuite.molview.pick import pick_sphere

#: the largest screen radius looked around (pixels); bigger discs are found by their centre
_MAX_REACH = 400.0


class ScreenPicker:
    def __init__(
        self,
        camera: Camera,
        centres: np.ndarray,
        radii: np.ndarray,
        width: float,
        height: float,
        slab: tuple[float, float] = (0.0, np.inf),
    ) -> None:
        self.camera = camera
        self.width, self.height = width, height
        xy, depth = (
            camera.project(centres, width, height)
            if len(centres)
            else (np.zeros((0, 2)), np.zeros(0))
        )
        keep = (depth > slab[0]) & (depth < slab[1])
        self._index = np.flatnonzero(keep)
        self._centres = np.asarray(centres, float)[keep]
        self._radii = np.asarray(radii, float)[keep]
        self._xy = xy[keep]
        pixels = camera.pixels_per_unit(1.0, height) / np.maximum(depth[keep], 1e-6)
        self._reach = float(min((self._radii * pixels).max(), _MAX_REACH)) if len(pixels) else 0
        self._tree = cKDTree(self._xy) if len(self._xy) else None

    def pick(self, x: float, y: float) -> int:
        """The index (into the centres given) of the nearest sphere under (x, y), or -1."""
        if self._tree is None:
            return -1
        near = self._tree.query_ball_point((x, y), self._reach + 1.0)
        if not near:
            return -1
        near = np.asarray(near)
        origin, direction = self.camera.ray(x, y, self.width, self.height)
        hit = pick_sphere(origin, direction, self._centres[near], self._radii[near])
        return -1 if hit is None else int(self._index[near[hit]])
