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
"""What a plot draws (:class:`PlotSpec`) and what is under the mouse (:class:`Hover`)."""

from __future__ import annotations

from dataclasses import dataclass

from chemsuite.plotting.style import DEFAULT_MAP

Range = tuple[float | None, float | None]


@dataclass(frozen=True)
class PlotSpec:
    kind: str = "scatter"  # or "histogram"
    x: str | None = None
    y: str | None = None
    color: str | None = None
    size: str | None = None
    bins: int | None = None  # None: automatic
    fit: bool = False  # scatter: least-squares line(s)
    curve: str | None = None  # histogram: "normal" or "lognormal" density overlay
    guide: str | None = None  # scatter: a key of GUIDES (BOILED-Egg, Golden Triangle, LipE)
    aggregate: str = "mean"  # bar: a key of AGGREGATES
    errors: bool = False  # bar: standard deviation whiskers (with "mean")
    columns: tuple[str, ...] = ()  # radar: the axes
    z: str | None = None  # 3D: the third axis (the landscape's height)
    drop: bool = False  # 3D scatter: lines from each point down to the floor
    cmap: str = DEFAULT_MAP  # number colours: a key of COLOR_MAPS
    reverse: bool = False  # the colour map runs high to low
    color_range: Range = (None, None)  # number colours span this (None: the data's own end)
    size_range: Range = (None, None)  # sizes span these values (None: the data's own end)
    point_sizes: tuple[float, float] = (4.0, 18.0)  # smallest and largest marker, in pixels
    structure: str | None = None  # constellation: the structure column
    summary: str = "mean"  # constellation: a key of SUMMARIES, over each core's compounds
    layout: str = "tsne"  # constellation: a key of MAPS
    min_compounds: int = 2  # constellation: cores shared by fewer are left out
    core_percent: int = 67  # constellation: a core keeps at least this % of the heavy atoms
    labels: bool = True  # constellation: name each analog series


@dataclass(frozen=True)
class Hover:
    """What is under the mouse: one row (``oid``) or a group of rows (a histogram bar)."""

    oid: int | None
    title: str
    lines: tuple[str, ...]
    smiles: str | None = None  # a structure to draw on the card (a constellation's core)
