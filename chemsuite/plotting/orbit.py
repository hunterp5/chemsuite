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
"""The camera of the 3D plots: it orbits a centre, z up, as pyqtgraph's GLViewWidget did.

The maths matches pyqtgraph's (horizontal field of view, distance, elevation and azimuth in
degrees), so the cameras saved with documents before the plots moved to QRhi look the same.
Matrices are row-major numpy arrays with OpenGL conventions; the renderer adds the graphics
backend's clip-space correction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .geometry import project

FOV = 60.0  # horizontal field of view, degrees
DEFAULT = {"distance": 4.6, "elevation": 22.0, "azimuth": -55.0}
WHEEL_BASE = 0.999  # distance factor per wheel unit (120 units a notch), as pyqtgraph


def _translate(x: float, y: float, z: float) -> np.ndarray:
    m = np.eye(4)
    m[:3, 3] = (x, y, z)
    return m


def _rotate_x(degrees: float) -> np.ndarray:
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    m = np.eye(4)
    m[1:3, 1:3] = ((c, -s), (s, c))
    return m


def _rotate_z(degrees: float) -> np.ndarray:
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    m = np.eye(4)
    m[0:2, 0:2] = ((c, -s), (s, c))
    return m


def _number(value: object) -> bool:
    return isinstance(value, int | float) and not isinstance(value, bool) and math.isfinite(value)


@dataclass
class OrbitCamera:
    distance: float = DEFAULT["distance"]
    elevation: float = DEFAULT["elevation"]
    azimuth: float = DEFAULT["azimuth"]
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))

    # -- matrices -------------------------------------------------------------------------------
    def view(self) -> np.ndarray:
        return (
            _translate(0.0, 0.0, -self.distance)
            @ _rotate_x(self.elevation - 90.0)
            @ _rotate_z(-(self.azimuth + 90.0))
            @ _translate(*(-self.center))
        )

    def clip_range(self) -> tuple[float, float]:
        """Near and far planes around the ±1 cube (tighter than pyqtgraph's, for depth
        precision; the picture is the same)."""
        reach = math.sqrt(3.0) + float(np.linalg.norm(self.center)) + 1.0
        return max(self.distance - reach, self.distance * 0.01), self.distance + reach

    def projection(self, aspect: float) -> np.ndarray:
        """Perspective for a view *aspect* (width / height) wide."""
        near, far = self.clip_range()
        r = near * math.tan(math.radians(FOV) / 2)
        t = r / max(aspect, 1e-6)
        return np.array(
            [
                [near / r, 0.0, 0.0, 0.0],
                [0.0, near / t, 0.0, 0.0],
                [0.0, 0.0, -(far + near) / (far - near), -2 * far * near / (far - near)],
                [0.0, 0.0, -1.0, 0.0],
            ]
        )

    def matrix(self, width: float, height: float) -> np.ndarray:
        return self.projection(width / max(height, 1.0)) @ self.view()

    def project(self, points: np.ndarray, width: float, height: float):
        """Where *points* are in a *width* × *height* view, in pixels from the top left, and
        whether each is in front of the camera."""
        return project(points, self.matrix(width, height), width, height)

    def eye_direction(self) -> np.ndarray:
        """The unit vector from the centre towards the camera."""
        e, a = math.radians(self.elevation), math.radians(self.azimuth)
        return np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])

    # -- moving ---------------------------------------------------------------------------------
    def turn(self, dx: float, dy: float) -> None:
        """A drag of (dx, dy) pixels: a degree a pixel; the elevation stops at the poles."""
        self.azimuth -= dx
        self.elevation = min(max(self.elevation + dy, -90.0), 90.0)

    def pan(self, dx: float, dy: float, width: float) -> None:
        """Slide the centre so the scene follows a drag of (dx, dy) pixels; up stays up."""
        scale = 2 * self.distance * math.tan(math.radians(FOV) / 2) / max(width, 1.0)
        a = math.radians(self.azimuth)
        right = np.array([-math.sin(a), math.cos(a), 0.0])
        self.center = self.center - right * dx * scale + np.array([0.0, 0.0, dy * scale])

    def zoom(self, wheel_units: float) -> None:
        """A wheel turn (120 units a notch; away from you zooms in)."""
        self.distance *= WHEEL_BASE**wheel_units

    # -- saving ---------------------------------------------------------------------------------
    def state(self) -> dict:
        return {
            "distance": float(self.distance),
            "elevation": float(self.elevation),
            "azimuth": float(self.azimuth),
            "center": [float(v) for v in self.center],
        }

    @classmethod
    def from_state(cls, state: object) -> OrbitCamera:
        """A camera from :meth:`state`; missing or bad values take the default."""
        camera = cls()
        if not isinstance(state, dict):
            return camera
        for key in DEFAULT:
            value = state.get(key)
            if _number(value) and (key != "distance" or value > 0):
                setattr(camera, key, float(value))
        center = state.get("center")
        if isinstance(center, list) and len(center) == 3 and all(map(_number, center)):
            camera.center = np.array(center, float)
        camera.elevation = min(max(camera.elevation, -90.0), 90.0)
        return camera
