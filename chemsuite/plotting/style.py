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
"""Colours for numeric scales and categories, and marker sizes."""

from __future__ import annotations

import numpy as np

# Viridis anchors (matplotlib's colormap, sampled at 9 points; perceptually uniform and
# readable for the common forms of colour blindness).
_VIRIDIS = np.array(
    [
        (68, 1, 84),
        (71, 44, 122),
        (59, 81, 139),
        (44, 113, 142),
        (33, 144, 141),
        (39, 173, 129),
        (92, 200, 99),
        (170, 220, 50),
        (253, 231, 37),
    ],
    dtype=np.float64,
)
#: Okabe-Ito palette, distinguishable with colour blindness; grey reserved for Other/blank.
CATEGORY_COLORS = [
    (0, 114, 178),
    (230, 159, 0),
    (0, 158, 115),
    (204, 121, 167),
    (86, 180, 233),
    (213, 94, 0),
    (240, 228, 66),
    (0, 0, 0),
    (117, 112, 179),
    (102, 166, 30),
    (166, 118, 29),
    (231, 41, 138),
]
GREY = (160, 160, 160)
DEFAULT_COLOR = (31, 119, 180)


def scale_colors(steps: int) -> list[tuple[int, int, int]]:
    """*steps* colours evenly spaced along the scale."""
    positions = np.linspace(0.0, len(_VIRIDIS) - 1, steps)
    lower = np.floor(positions).astype(int).clip(0, len(_VIRIDIS) - 2)
    frac = (positions - lower)[:, None]
    rgb = _VIRIDIS[lower] * (1 - frac) + _VIRIDIS[lower + 1] * frac
    return [tuple(int(round(c)) for c in row) for row in rgb]


def scale_bins(values: np.ndarray, steps: int, lo: float | None = None, hi: float | None = None):
    """``(bin per value, lo, hi)``: bin 0..steps-1 along [lo, hi]; -1 for NaN."""
    finite = values[np.isfinite(values)]
    if lo is None:
        lo = float(finite.min()) if len(finite) else 0.0
    if hi is None:
        hi = float(finite.max()) if len(finite) else 1.0
    span = hi - lo if hi > lo else 1.0
    with np.errstate(invalid="ignore"):
        bins = np.floor((values - lo) / span * steps)
    bins = np.where(np.isfinite(bins), np.clip(bins, 0, steps - 1), -1).astype(np.int64)
    return bins, lo, hi


def category_color(index: int, label: str, *, special: tuple[str, ...]) -> tuple[int, int, int]:
    if label in special:
        return GREY
    return CATEGORY_COLORS[index % len(CATEGORY_COLORS)]


def marker_sizes(values: np.ndarray, small: float = 4.0, large: float = 16.0) -> np.ndarray:
    """Pixel sizes scaled linearly between *small* and *large*; NaN gets *small*."""
    finite = values[np.isfinite(values)]
    if not len(finite):
        return np.full(len(values), small)
    lo, hi = float(finite.min()), float(finite.max())
    span = hi - lo if hi > lo else 1.0
    with np.errstate(invalid="ignore"):
        sizes = small + (values - lo) / span * (large - small)
    return np.where(np.isfinite(sizes), sizes, small)
