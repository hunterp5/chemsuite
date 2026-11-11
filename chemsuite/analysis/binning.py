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
"""Assign numbers to bins: a range label ("[5, 6)") and a class name ("Medium") per row."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise

import numpy as np

MAX_BINS = 256
METHODS = {"equal_width": "Equal width", "equal_count": "Equal count", "custom": "Custom edges"}
_PRESETS: dict[int, tuple[str, ...]] = {
    2: ("Low", "High"),
    3: ("Low", "Medium", "High"),
    4: ("Low", "Medium", "High", "Very high"),
}


@dataclass(frozen=True)
class Interval:
    lo: float
    hi: float
    right_closed: bool = False  # every bin includes its low edge

    def label(self) -> str:
        if math.isinf(self.lo) and not math.isinf(self.hi):
            return f"< {_edge(self.hi)}"
        if math.isinf(self.hi) and not math.isinf(self.lo):
            return f"≥ {_edge(self.lo)}"
        if self.lo == self.hi:
            return _edge(self.lo)
        return f"[{_edge(self.lo)}, {_edge(self.hi)}{']' if self.right_closed else ')'}"


@dataclass(frozen=True)
class BinSpec:
    method: str = "equal_width"
    count: int = 5
    width: float | None = None  # equal width: a bin width instead of a count
    minimum: float | None = None  # equal width: limits (values outside get no bin)
    maximum: float | None = None
    edges: tuple[float, ...] = ()  # custom
    open_ends: bool = True  # custom: add "< first" and "≥ last" bins
    names: tuple[str, ...] = ()  # one class name per bin, low to high; empty: defaults


@dataclass(frozen=True)
class Binned:
    index: np.ndarray  # int per value: its bin, -1 for none
    intervals: tuple[Interval, ...]
    class_names: tuple[str, ...]

    @property
    def range_names(self) -> tuple[str, ...]:
        return tuple(i.label() for i in self.intervals)

    def labels(self, names: Sequence[str]) -> list[str]:
        lookup = [*names, ""]  # index -1 → ""
        return [lookup[i] for i in self.index.tolist()]

    def counts(self) -> np.ndarray:
        return np.bincount(self.index[self.index >= 0], minlength=len(self.intervals))


def _edge(value: float) -> str:
    text = f"{value:.6g}"
    return "0" if text == "-0" else text


def default_names(count: int) -> tuple[str, ...]:
    return _PRESETS.get(count, tuple(f"Bin {i + 1}" for i in range(count)))


def parse_numbers(text: str) -> tuple[float, ...]:
    """Edges from "5, 6, 7" (commas, semicolons or spaces); must increase."""
    parts = text.replace(";", " ").replace(",", " ").split()
    if not parts:
        raise ValueError("Enter bin edges, for example 5, 6, 7.")
    try:
        edges = tuple(float(p) for p in parts)
    except ValueError as err:
        raise ValueError(f"Could not read the edges: {err}") from None
    if any(not math.isfinite(e) for e in edges):
        raise ValueError("Edges must be finite numbers.")
    if any(b <= a for a, b in pairwise(edges)):
        raise ValueError("Edges must increase from left to right.")
    return edges


def parse_names(text: str) -> tuple[str, ...]:
    parts = [p.strip() for p in text.split(",")] if text.strip() else []
    if any(not p for p in parts):
        raise ValueError("Bin names cannot be blank; separate them with commas.")
    return tuple(parts)


def _closed(edges: Sequence[float]) -> tuple[Interval, ...]:
    """Bins [e0, e1), [e1, e2), …, [e(n-1), en]: the last one includes its top edge."""
    n = len(edges) - 1
    return tuple(Interval(edges[i], edges[i + 1], i == n - 1) for i in range(n))


def intervals_for(values: np.ndarray, spec: BinSpec) -> tuple[Interval, ...]:  # noqa: C901 - one rule per binning mode
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x)]
    if spec.method == "custom":
        edges = list(spec.edges)
        if any(b <= a for a, b in pairwise(edges)):
            raise ValueError("Edges must increase from left to right.")
        if spec.open_ends:
            if not edges:
                raise ValueError("Enter at least one edge.")
            inner = [Interval(a, b) for a, b in pairwise(edges)]
            out = (Interval(-math.inf, edges[0]), *inner, Interval(edges[-1], math.inf))
        else:
            if len(edges) < 2:
                raise ValueError("Enter at least two edges, or include values beyond them.")
            out = _closed(edges)
    elif not len(x):
        raise ValueError("There are no numbers to bin.")
    elif spec.method == "equal_count":
        out = _equal_count(x, spec.count)
    elif spec.method == "equal_width":
        lo = float(x.min()) if spec.minimum is None else spec.minimum
        hi = float(x.max()) if spec.maximum is None else spec.maximum
        if hi < lo:
            raise ValueError("The minimum is above the maximum.")
        if spec.width is not None:
            if not spec.width > 0:
                raise ValueError("The bin width must be above 0.")
            if spec.minimum is None:  # align to multiples of the width, like a histogram
                lo = math.floor(lo / spec.width) * spec.width
            steps = max(1, math.ceil((hi - lo) / spec.width - 1e-9))
            if steps > MAX_BINS:
                raise ValueError(f"That width makes more than {MAX_BINS} bins.")
            edges = [lo + i * spec.width for i in range(steps)] + [
                hi if spec.maximum is not None else lo + steps * spec.width
            ]
            out = _closed(edges)
        elif hi == lo:
            out = (Interval(lo, hi, True),)
        else:
            _check_count(spec.count)
            out = _closed(np.linspace(lo, hi, spec.count + 1).tolist())
    else:
        raise ValueError(f"unknown binning method: {spec.method!r}")
    if len(out) > MAX_BINS:
        raise ValueError(f"Use at most {MAX_BINS} bins.")
    return out


def _check_count(count: int) -> None:
    if not 2 <= count <= MAX_BINS:
        raise ValueError(f"Use between 2 and {MAX_BINS} bins.")


def _equal_count(x: np.ndarray, count: int) -> tuple[Interval, ...]:
    """Groups of nearly equal size; equal values always share a bin (so ties can mean fewer)."""
    _check_count(count)
    ordered = np.sort(x)
    n = len(ordered)
    group = np.minimum(count - 1, np.arange(n) * count // n)
    uniques, first = np.unique(ordered, return_index=True)
    group_of_unique = group[first]  # a value's group is where it first appears
    starts = np.flatnonzero(np.diff(group_of_unique, prepend=-1))
    lows = uniques[starts].tolist()
    edges = [*lows, float(uniques[-1])]
    if len(lows) == 1:
        return (Interval(lows[0], edges[-1], True),)
    return _closed(edges)


def assign(values: np.ndarray, intervals: Sequence[Interval]) -> np.ndarray:
    """Bin index per value (-1: missing or outside every bin)."""
    x = np.asarray(values, dtype=np.float64)
    lows = np.array([i.lo for i in intervals])
    index = np.searchsorted(lows, x, side="right") - 1  # last bin whose low edge ≤ x
    index = np.where(np.isfinite(x) & (index >= 0), index, -1)
    safe = np.clip(index, 0, None)
    highs = np.array([i.hi for i in intervals])[safe]
    closed = np.array([i.right_closed for i in intervals])[safe]
    lows_of = lows[safe]
    inside = np.where(closed, x <= highs, x < highs) | ((lows_of == highs) & (x == lows_of))
    return np.where((index >= 0) & inside, index, -1)


def bin_values(values: np.ndarray, spec: BinSpec) -> Binned:
    intervals = intervals_for(values, spec)
    names = spec.names or default_names(len(intervals))
    if len(names) != len(intervals):
        ranges = ", ".join(i.label() for i in intervals[:8]) + ("…" if len(intervals) > 8 else "")
        raise ValueError(
            f"These settings make {len(intervals)} bins ({ranges}); "
            f"enter {len(intervals)} names, low to high, or leave the names blank."
        )
    return Binned(assign(values, intervals), intervals, tuple(names))
