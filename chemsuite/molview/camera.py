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
"""The 3D Viewer's camera: turn, pan, zoom and fit, and the matrices the renderer uses.

The camera circles a centre point at a distance; ``rotation`` takes world directions to view
directions (view space looks down -z, y up), so the eye sits at ``centre + rotation.T @ (0, 0,
distance)``. The mouse does what it does in 3Dmol.js: drag turns, right- or middle-drag pans,
the wheel zooms.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

FOV = 30.0  # vertical field of view, degrees
TURN_PER_PIXEL = 0.01  # radians
MIN_RADIUS = 2.0  # Å: a single atom still gets some room
MARGIN = 1.08  # a little room round a fitted scene
FIT_PASSES = 6  # refinements of a fit in perspective


def _rot(axis: np.ndarray, angle: float) -> np.ndarray:
    x, y, z = axis / np.linalg.norm(axis)
    c, s = math.cos(angle), math.sin(angle)
    return np.array(
        [
            [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
        ]
    )


@dataclass
class Camera:
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))
    centre: np.ndarray = field(default_factory=lambda: np.zeros(3))
    distance: float = 20.0
    radius: float = MIN_RADIUS  # of the scene, for the near and far planes
    fov: float = FOV

    # -- placing ----------------------------------------------------------------------------
    def fit(self, positions: np.ndarray, radii: np.ndarray | None = None, aspect: float = 1.0):
        """Look at *positions* from the current direction, just far enough to see them all."""
        if not len(positions):
            return
        lo, hi = positions.min(axis=0), positions.max(axis=0)
        self.centre = ((lo + hi) / 2).astype(float)
        pad = np.zeros(len(positions)) if radii is None else np.asarray(radii, float)
        offset = (positions - self.centre) @ self.rotation.T  # view directions, centred
        self.radius = max(float((np.linalg.norm(offset, axis=1) + pad).max()), MIN_RADIUS)
        tan_y = math.tan(math.radians(self.fov) / 2)
        tan_x = tan_y * max(aspect, 1e-6)
        # each atom's edge must lie inside the view cone from an eye at (0, 0, distance)
        need_y = (np.abs(offset[:, 1]) + pad) / tan_y + offset[:, 2] + pad
        need_x = (np.abs(offset[:, 0]) + pad) / tan_x + offset[:, 2] + pad
        self.distance = max(float(max(need_x.max(), need_y.max())) * MARGIN, self.radius * 1.2)
        self._refine_fit(positions, pad, aspect)

    def _refine_fit(self, positions: np.ndarray, pad: np.ndarray, aspect: float) -> None:
        """Centre the atoms' drawn outline on the screen and make it fill the view (less the
        margin) in its tighter direction, as drawn in perspective; a few passes converge."""
        tan_y = math.tan(math.radians(self.fov) / 2)
        for _ in range(FIT_PASSES):
            xy, depth = self.project(positions, aspect, 1.0)  # in view heights
            size = pad / (2 * np.maximum(depth, 1e-6) * tan_y)
            lo = (xy - size[:, None]).min(axis=0)
            hi = (xy + size[:, None]).max(axis=0)
            shift = (lo + hi) / 2 - np.array([aspect / 2, 0.5])
            per_height = 2 * self.distance * tan_y  # Å per view height at the centre's depth
            self.centre = self.centre + self.rotation.T @ np.array(
                [shift[0] * per_height, -shift[1] * per_height, 0.0]
            )
            fill = max((hi[0] - lo[0]) / aspect, hi[1] - lo[1])
            self.distance *= max(fill * MARGIN, 1e-3)
            nearest = float(((positions - self.centre) @ self.rotation.T)[:, 2].max() + pad.max())
            self.distance = max(self.distance, nearest * 1.05 + 0.5)  # the eye stays outside
        offset = positions - self.centre
        self.radius = max(float((np.linalg.norm(offset, axis=1) + pad).max()), MIN_RADIUS)

    def orient(self, positions: np.ndarray, aspect: float = 1.0) -> None:
        """Turn to the scene's principal axes: its longest along the view's longer side, its
        thinnest pointing at the eye (as PyMOL's ``orient``)."""
        if len(positions) < 2:
            return
        centred = np.asarray(positions, float) - np.mean(positions, axis=0)
        _values, vectors = np.linalg.eigh(centred.T @ centred)  # ascending
        longest, middle, thinnest = vectors[:, 2], vectors[:, 1], vectors[:, 0]
        across, up = (longest, middle) if aspect >= 1.0 else (middle, longest)
        rotation = np.stack([across, up, thinnest])
        if np.linalg.det(rotation) < 0:
            rotation[2] = -rotation[2]
        self.rotation = rotation

    def turn(self, dx: float, dy: float) -> None:
        """Turn the scene as a mouse drag of (dx, dy) pixels does."""
        spin = _rot(np.array([0.0, 1.0, 0.0]), dx * TURN_PER_PIXEL)
        tilt = _rot(np.array([1.0, 0.0, 0.0]), dy * TURN_PER_PIXEL)
        self.rotation = tilt @ spin @ self.rotation

    def pan(self, dx: float, dy: float, height: float) -> None:
        """Move the scene with the mouse: (dx, dy) pixels in a view *height* pixels high."""
        per_pixel = 2 * self.distance * math.tan(math.radians(self.fov) / 2) / max(height, 1.0)
        self.centre = self.centre - self.rotation.T @ np.array([dx * per_pixel, -dy * per_pixel, 0])

    def zoom(self, factor: float) -> None:
        """Move closer (factor < 1) or further away."""
        self.distance = float(np.clip(self.distance * factor, 0.5, 2000.0))

    # -- matrices ---------------------------------------------------------------------------
    def view(self) -> np.ndarray:
        m = np.eye(4)
        m[:3, :3] = self.rotation
        m[:3, 3] = -self.rotation @ self.centre - np.array([0.0, 0.0, self.distance])
        return m

    def clip_range(self) -> tuple[float, float]:
        near = max(self.distance - self.radius * 1.5, self.distance * 0.02, 0.05)
        far = self.distance + self.radius * 1.5
        return near, far

    def projection(self, aspect: float) -> np.ndarray:
        """OpenGL-style perspective (clip z from -1 to 1); the renderer corrects it per API."""
        near, far = self.clip_range()
        f = 1.0 / math.tan(math.radians(self.fov) / 2)
        m = np.zeros((4, 4))
        m[0, 0] = f / max(aspect, 1e-6)
        m[1, 1] = f
        m[2, 2] = (far + near) / (near - far)
        m[2, 3] = 2 * far * near / (near - far)
        m[3, 2] = -1.0
        return m

    def eye(self) -> np.ndarray:
        return self.centre + self.rotation.T @ np.array([0.0, 0.0, self.distance])

    # -- screen <-> world -------------------------------------------------------------------
    def ray(self, x: float, y: float, width: float, height: float):
        """The eye and the unit direction through pixel (x, y) of a width x height view."""
        t = math.tan(math.radians(self.fov) / 2)
        aspect = width / max(height, 1.0)
        d = np.array([(2 * x / width - 1) * t * aspect, (1 - 2 * y / height) * t, -1.0])
        d = self.rotation.T @ (d / np.linalg.norm(d))
        return self.eye(), d

    def project(self, points: np.ndarray, width: float, height: float):
        """Pixel positions (n, 2) and distances in front of the eye (n,) of world *points*."""
        v = (np.asarray(points, float) - self.centre) @ self.rotation.T
        v[:, 2] -= self.distance
        depth = -v[:, 2]
        t = math.tan(math.radians(self.fov) / 2)
        aspect = width / max(height, 1.0)
        safe = np.where(depth > 1e-6, depth, 1e-6)
        x = (v[:, 0] / (safe * t * aspect) + 1) * width / 2
        y = (1 - v[:, 1] / (safe * t)) * height / 2
        return np.stack([x, y], axis=1), depth

    def pixels_per_unit(self, depth: float, height: float) -> float:
        """How many pixels one Å spans at *depth* in front of the eye."""
        return height / (2 * max(depth, 1e-6) * math.tan(math.radians(self.fov) / 2))

    # -- saving -----------------------------------------------------------------------------
    def state(self) -> dict:
        return {
            "rotation": [float(v) for v in self.rotation.ravel()],
            "centre": [float(v) for v in self.centre],
            "distance": float(self.distance),
            "radius": float(self.radius),
        }

    @classmethod
    def from_state(cls, state: dict) -> Camera:
        camera = cls()
        try:
            camera.rotation = np.array(state["rotation"], float).reshape(3, 3)
            camera.centre = np.array(state["centre"], float).reshape(3)
            camera.distance = float(state["distance"])
            camera.radius = float(state.get("radius", MIN_RADIUS))
        except (KeyError, TypeError, ValueError):
            return cls()
        return camera
