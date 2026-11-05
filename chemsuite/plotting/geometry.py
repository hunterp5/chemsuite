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
"""Which points a selection gesture covers (vectorised; 200k points in milliseconds)."""

from __future__ import annotations

import numpy as np


def box_mask(
    x: np.ndarray, y: np.ndarray, x0: float, x1: float, y0: float, y1: float
) -> np.ndarray:
    lo_x, hi_x = sorted((x0, x1))
    lo_y, hi_y = sorted((y0, y1))
    return (x >= lo_x) & (x <= hi_x) & (y >= lo_y) & (y <= hi_y)


def lasso_mask(x: np.ndarray, y: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    """Points inside *polygon* (``k × 2`` vertices, implicitly closed), by ray casting.

    One pass per edge over all points, so the cost is O(edges × points) in numpy.
    """
    polygon = np.asarray(polygon, dtype=np.float64)
    inside = np.zeros(len(x), dtype=np.bool_)
    if len(polygon) < 3 or not len(x):
        return inside
    # Only test points in the polygon's bounding box.
    candidate = box_mask(
        x, y, polygon[:, 0].min(), polygon[:, 0].max(), polygon[:, 1].min(), polygon[:, 1].max()
    )
    idx = np.flatnonzero(candidate)
    px, py = x[idx], y[idx]
    hit = np.zeros(len(idx), dtype=np.bool_)
    xj, yj = polygon[-1]
    for xi, yi in polygon:
        crosses = (yi > py) != (yj > py)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_cross = (xj - xi) * (py - yi) / (yj - yi) + xi
        hit ^= crosses & (px < x_cross)
        xj, yj = xi, yi
    inside[idx] = hit
    return inside


def nearest(
    x: np.ndarray, y: np.ndarray, px: float, py: float, sx: float, sy: float, radius: float
) -> int | None:
    """Index of the point nearest (px, py) within *radius* pixels; *sx*, *sy* are data units
    per pixel (so the distance is measured on screen, not in data units)."""
    if not len(x) or sx <= 0 or sy <= 0:
        return None
    d2 = ((x - px) / sx) ** 2 + ((y - py) / sy) ** 2
    i = int(np.argmin(d2))
    return i if d2[i] <= radius * radius else None


def jitter(oids: np.ndarray, width: float) -> np.ndarray:
    """A stable offset in ``[-width/2, width/2]`` per row, so points do not jump on redraw."""
    mixed = (np.asarray(oids, dtype=np.uint64) * np.uint64(2654435761)) % np.uint64(2**32)
    return (mixed.astype(np.float64) / 2**32 - 0.5) * width
