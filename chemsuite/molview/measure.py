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
"""Distances, angles and dihedrals between atoms (Å and degrees)."""

from __future__ import annotations

import math

import numpy as np


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(b, float) - np.asarray(a, float)))


def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """The angle at *b*."""
    u = np.asarray(a, float) - b
    v = np.asarray(c, float) - b
    cos = u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-12)
    return math.degrees(math.acos(float(np.clip(cos, -1.0, 1.0))))


def dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    """The torsion a-b-c-d, from -180 to 180."""
    b0 = np.asarray(a, float) - b
    b1 = np.asarray(c, float) - b
    b2 = np.asarray(d, float) - c
    b1 = b1 / max(np.linalg.norm(b1), 1e-12)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return math.degrees(math.atan2(float(np.cross(b1, v) @ w), float(v @ w)))


def describe(points: list[np.ndarray]) -> str:
    """The measurement two, three or four picked atoms make."""
    if len(points) == 2:
        return f"{distance(*points):.2f} Å"
    if len(points) == 3:
        return f"{angle(*points):.1f}°"
    if len(points) == 4:
        return f"{dihedral(*points):.1f}°"
    return ""
