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
"""The 3D plots' orbit camera: the same picture pyqtgraph drew, and moving it."""

from __future__ import annotations

import math

import numpy as np
import pytest

from chemsuite.plotting.orbit import DEFAULT, FOV, OrbitCamera

CUBE = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], float)


def _pyqtgraph_matrix(distance, elevation, azimuth, width, height) -> np.ndarray:
    """pyqtgraph's GLViewWidget.projectionMatrix() × viewMatrix(), as it computed them."""
    from PySide6.QtGui import QMatrix4x4

    near, far = distance * 0.001, distance * 1000.0
    r = near * math.tan(0.5 * math.radians(FOV))
    t = r * height / width
    projection = QMatrix4x4()
    projection.frustum(-r, r, -t, t, near, far)
    view = QMatrix4x4()
    view.translate(0.0, 0.0, -distance)
    view.rotate(elevation - 90, 1, 0, 0)
    view.rotate(azimuth + 90, 0, 0, -1)
    return np.array((projection * view).data(), float).reshape(4, 4).T


@pytest.mark.parametrize(
    ("distance", "elevation", "azimuth", "size"),
    [(4.6, 22.0, -55.0, (640, 480)), (3.0, -40.0, 130.0, (300, 700)), (8.0, 90.0, 0.0, (500, 500))],
)
def test_saved_views_look_as_they_did_with_pyqtgraph(distance, elevation, azimuth, size) -> None:
    camera = OrbitCamera(distance=distance, elevation=elevation, azimuth=azimuth)
    ours, front = camera.project(CUBE, *size)
    theirs = _pyqtgraph_matrix(distance, elevation, azimuth, *size)
    from chemsuite.plotting.geometry import project

    expected, _ = project(CUBE, theirs, *size)
    assert front.all()
    np.testing.assert_allclose(ours, expected, atol=1e-3)  # Qt computes in float32


def test_the_cube_stays_between_the_clip_planes() -> None:
    for distance in (2.5, 4.6, 40.0):
        camera = OrbitCamera(distance=distance)
        clip = np.column_stack([CUBE, np.ones(8)]) @ camera.matrix(400, 300).T
        depth = clip[:, 2] / clip[:, 3]
        assert (np.abs(depth) < 1).all()


def test_turning_stops_at_the_poles() -> None:
    camera = OrbitCamera()
    camera.turn(10, 500)
    assert camera.elevation == 90.0 and camera.azimuth == DEFAULT["azimuth"] - 10


def test_panning_moves_the_scene_with_the_pointer() -> None:
    camera = OrbitCamera(elevation=0.0)
    (before,), _ = camera.project(np.zeros((1, 3)), 400, 300)
    camera.pan(30, 20, 400)
    (after,), _ = camera.project(np.zeros((1, 3)), 400, 300)
    np.testing.assert_allclose(after - before, (30, 20), atol=1e-6)


def test_the_wheel_zooms_like_pyqtgraph() -> None:
    camera = OrbitCamera(distance=4.0)
    camera.zoom(120)
    assert camera.distance == pytest.approx(4.0 * 0.999**120)


def test_state_round_trips_and_bad_values_fall_back() -> None:
    camera = OrbitCamera(distance=3.0, elevation=10.0, azimuth=5.0, center=np.array([0.1, 0, 0]))
    again = OrbitCamera.from_state(camera.state())
    assert again.state() == camera.state()
    old = OrbitCamera.from_state({"distance": 5.0, "elevation": 30.0, "azimuth": 10.0})
    assert old.distance == 5.0 and list(old.center) == [0, 0, 0]  # saved before "center"
    bad = OrbitCamera.from_state({"distance": -1, "elevation": "x", "center": [1, None, 2]})
    assert bad.state() == OrbitCamera().state()
    assert OrbitCamera.from_state(None).state() == OrbitCamera().state()
