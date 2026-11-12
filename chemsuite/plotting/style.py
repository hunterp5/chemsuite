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
"""Colours for numeric scales (colour maps) and categories, and marker sizes."""

from __future__ import annotations

import numpy as np

#: Colour maps for number scales: key -> (title, anchor colours, low to high), interpolated
#: linearly. Viridis, Plasma, Inferno, Magma and Cividis are matplotlib's (CC0); Turbo is
#: Google's (Mikhailov 2019, Apache-2.0); Rainbow (R1) and Blue–Red (D1) are Peter Kovesi's
#: CET maps (CC-BY 4.0, colorcet.com). All are sampled at 17 points from the data pyqtgraph
#: ships. Jet is MATLAB's classic layout; Red–Yellow follows the constellation plots of
#: Naveja & Medina-Franco (2019).
COLOR_MAPS: dict[str, tuple[str, list[tuple[int, int, int]]]] = {
    "viridis": (
        "Viridis",
        [
            (68, 1, 84),
            (72, 24, 106),
            (71, 45, 123),
            (66, 64, 134),
            (59, 82, 139),
            (51, 99, 141),
            (44, 114, 142),
            (38, 130, 142),
            (33, 145, 140),
            (31, 159, 136),
            (39, 173, 129),
            (61, 188, 116),
            (92, 200, 99),
            (129, 211, 77),
            (170, 220, 50),
            (213, 226, 26),
            (253, 231, 37),
        ],
    ),
    "plasma": (
        "Plasma",
        [
            (13, 8, 135),
            (49, 5, 151),
            (76, 2, 161),
            (102, 0, 167),
            (126, 3, 168),
            (149, 17, 161),
            (170, 35, 149),
            (188, 53, 135),
            (204, 71, 120),
            (217, 88, 106),
            (229, 107, 93),
            (240, 127, 79),
            (248, 148, 65),
            (253, 171, 51),
            (253, 195, 40),
            (249, 221, 37),
            (240, 249, 33),
        ],
    ),
    "inferno": (
        "Inferno",
        [
            (0, 0, 4),
            (11, 7, 36),
            (33, 12, 74),
            (61, 9, 101),
            (87, 16, 110),
            (113, 25, 110),
            (138, 34, 106),
            (163, 44, 97),
            (188, 55, 84),
            (208, 69, 69),
            (227, 89, 51),
            (241, 113, 31),
            (249, 140, 10),
            (252, 170, 15),
            (249, 201, 50),
            (242, 232, 101),
            (252, 255, 164),
        ],
    ),
    "magma": (
        "Magma",
        [
            (0, 0, 4),
            (10, 8, 34),
            (29, 17, 71),
            (54, 16, 107),
            (81, 18, 124),
            (106, 28, 129),
            (131, 38, 129),
            (156, 46, 127),
            (183, 55, 121),
            (207, 64, 112),
            (229, 80, 100),
            (244, 105, 92),
            (251, 135, 97),
            (254, 165, 113),
            (254, 194, 135),
            (253, 224, 161),
            (252, 253, 191),
        ],
    ),
    "cividis": (
        "Cividis",
        [
            (0, 32, 77),
            (0, 44, 105),
            (6, 54, 110),
            (45, 66, 108),
            (65, 77, 107),
            (82, 89, 108),
            (97, 100, 111),
            (111, 112, 115),
            (125, 124, 120),
            (139, 135, 121),
            (155, 148, 119),
            (172, 161, 116),
            (188, 175, 111),
            (206, 188, 104),
            (224, 203, 94),
            (242, 218, 80),
            (255, 234, 70),
        ],
    ),
    "turbo": (
        "Turbo",
        [
            (48, 18, 59),
            (64, 64, 162),
            (70, 107, 227),
            (66, 148, 255),
            (40, 188, 235),
            (24, 221, 194),
            (50, 242, 152),
            (109, 254, 98),
            (164, 252, 60),
            (203, 237, 52),
            (236, 209, 58),
            (253, 174, 53),
            (251, 129, 34),
            (236, 83, 15),
            (210, 49, 5),
            (172, 23, 1),
            (122, 4, 3),
        ],
    ),
    "rainbow": (
        "Rainbow",
        [
            (0, 48, 246),
            (25, 94, 190),
            (39, 122, 139),
            (62, 143, 90),
            (71, 163, 37),
            (114, 177, 17),
            (161, 188, 24),
            (204, 198, 31),
            (241, 202, 36),
            (250, 185, 34),
            (248, 160, 29),
            (245, 133, 24),
            (241, 105, 21),
            (242, 86, 47),
            (252, 101, 111),
            (255, 123, 179),
            (253, 146, 250),
        ],
    ),
    "jet": (
        "Jet",
        [
            (0, 0, 128),
            (0, 0, 255),
            (0, 128, 255),
            (0, 255, 255),
            (128, 255, 128),
            (255, 255, 0),
            (255, 128, 0),
            (255, 0, 0),
            (128, 0, 0),
        ],
    ),
    "red_yellow": (
        "Red–Yellow",
        [(200, 0, 0), (235, 60, 0), (250, 130, 0), (255, 200, 20), (255, 245, 60)],
    ),
    "coolwarm": (
        "Blue–Red",
        [
            (33, 81, 219),
            (82, 98, 222),
            (112, 117, 226),
            (138, 136, 229),
            (160, 156, 232),
            (182, 177, 234),
            (202, 197, 237),
            (221, 219, 239),
            (238, 231, 233),
            (242, 215, 208),
            (240, 190, 178),
            (237, 165, 148),
            (231, 140, 119),
            (224, 115, 91),
            (215, 88, 63),
            (204, 57, 36),
            (192, 2, 6),
        ],
    ),
    "greys": ("Greys", [(250, 250, 250), (20, 20, 20)]),
}
DEFAULT_MAP = "viridis"

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


def scale_colors(
    steps: int, cmap: str = DEFAULT_MAP, reverse: bool = False
) -> list[tuple[int, int, int]]:
    """*steps* colours evenly spaced along colour map *cmap* (high to low when *reverse*)."""
    anchors = np.array(COLOR_MAPS.get(cmap, COLOR_MAPS[DEFAULT_MAP])[1], dtype=np.float64)
    if reverse:
        anchors = anchors[::-1]
    positions = np.linspace(0.0, len(anchors) - 1, max(steps, 1))
    lower = np.floor(positions).astype(int).clip(0, len(anchors) - 2)
    frac = (positions - lower)[:, None]
    rgb = anchors[lower] * (1 - frac) + anchors[lower + 1] * frac
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


def marker_sizes(
    values: np.ndarray,
    small: float = 4.0,
    large: float = 16.0,
    lo: float | None = None,
    hi: float | None = None,
) -> np.ndarray:
    """Pixel sizes scaled linearly between *small* and *large* along [lo, hi] (the values'
    own range where not given; values beyond it get the end sizes); NaN gets *small*."""
    finite = values[np.isfinite(values)]
    if not len(finite):
        return np.full(len(values), small)
    lo = float(finite.min()) if lo is None else lo
    hi = float(finite.max()) if hi is None else hi
    span = hi - lo if hi > lo else 1.0
    with np.errstate(invalid="ignore"):
        sizes = small + np.clip((values - lo) / span, 0.0, 1.0) * (large - small)
    return np.where(np.isfinite(sizes), sizes, small)


def parse_range(text: str) -> float | None:
    """A range end typed by hand: a number, or None for empty text or anything else (Auto)."""
    try:
        value = float(text.strip().replace(",", "."))
    except ValueError:
        return None
    return value if np.isfinite(value) else None
