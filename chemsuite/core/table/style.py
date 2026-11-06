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

"""How a column is shown: number precision, the log10 marker and cell colouring.

A :class:`ColumnStyle` belongs to its column object, so it moves with renames, comes back with
an undone delete and is saved in the document. Changing it is an undoable table edit
(:meth:`Table.set_style <chemsuite.core.table.table.Table.set_style>`).

* ``precision``: every number in the column is shown with this many decimals. The stored
  values keep their full precision, so sorting, filtering and tools are unaffected, and
  ``None`` brings back each cell's own format.
* ``log10``: the values were converted with "Logarithmic" (``log10(x)``) and converting back
  gives ``10**x``; the header says so.
* ``colour``: a :class:`ColourRule` for cell backgrounds.
"""

from __future__ import annotations

import math
import zlib
from dataclasses import asdict, dataclass, replace
from typing import Any

GRADIENT = "gradient"  # low -> high
GRADIENT3 = "gradient3"  # low -> mid -> high
CATEGORIES = "categories"  # one colour per distinct value
COLOUR_MODES = (GRADIENT, GRADIENT3, CATEGORIES)

#: blue, amber, red: the low/mid/high default; reads on light and dark backgrounds
LOW, MID, HIGH = 0x3077F2, 0xF5D154, 0xEC4949


@dataclass(frozen=True)
class ColourRule:
    """Cell background colours. Bounds left at ``None`` follow the column's data."""

    mode: str = GRADIENT
    low: int = LOW  # 0xRRGGBB
    mid: int = MID
    high: int = HIGH
    minimum: float | None = None  # None: the column's smallest value
    middle: float | None = None  # None: halfway between minimum and maximum
    maximum: float | None = None  # None: the column's largest value
    opacity: int = 96  # 0-255

    def __post_init__(self) -> None:
        if self.mode not in COLOUR_MODES:
            raise ValueError(f"unknown colour mode {self.mode!r}")

    @property
    def numeric(self) -> bool:
        return self.mode in (GRADIENT, GRADIENT3)

    def bounds(self, data_range: tuple[float, float] | None) -> tuple[float, float, float] | None:
        """``(minimum, middle, maximum)`` with the data filling in what the rule leaves open."""
        lo = self.minimum if self.minimum is not None else (data_range or (None, None))[0]
        hi = self.maximum if self.maximum is not None else (data_range or (None, None))[1]
        if lo is None or hi is None:
            return None
        mid = self.middle if self.middle is not None else (lo + hi) / 2
        return lo, min(max(mid, lo), hi), hi

    def number_colour(self, value: float, bounds: tuple[float, float, float]) -> int:
        """0xAARRGGBB for *value*; values outside the bounds take the end colour."""
        lo, mid, hi = bounds
        if self.mode == GRADIENT3:
            if value <= mid:
                rgb = _mix(self.low, self.mid, _fraction(value, lo, mid))
            else:
                rgb = _mix(self.mid, self.high, _fraction(value, mid, hi))
        else:
            rgb = _mix(self.low, self.high, _fraction(value, lo, hi))
        return (self.opacity << 24) | rgb

    def text_colour(self, text: str) -> int:
        """0xAARRGGBB for a category: the same text always gets the same colour."""
        hue = zlib.crc32(text.encode("utf-8")) % 360
        return (self.opacity << 24) | _hsl(hue, 0.55, 0.62)

    def to_json(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> ColourRule:
        fields = cls.__dataclass_fields__
        return cls(**{k: v for k, v in data.items() if k in fields})


@dataclass(frozen=True)
class ColumnStyle:
    precision: int | None = None
    log10: bool = False
    colour: ColourRule | None = None

    def __post_init__(self) -> None:
        if self.precision is not None and not 0 <= self.precision <= MAX_PRECISION:
            raise ValueError(f"precision must be 0-{MAX_PRECISION}")

    @property
    def is_default(self) -> bool:
        return self == DEFAULT_STYLE

    def for_kind(self, kind: str) -> ColumnStyle:
        """The parts of this style that still mean something after a change of kind."""
        if kind == "number":
            return self
        colour = self.colour if self.colour and not self.colour.numeric else None
        return ColumnStyle(colour=colour) if kind == "text" else DEFAULT_STYLE

    def with_(self, **changes: Any) -> ColumnStyle:
        return replace(self, **changes)

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        if self.precision is not None:
            out["precision"] = self.precision
        if self.log10:
            out["log10"] = True
        if self.colour is not None:
            out["colour"] = self.colour.to_json()
        return out

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> ColumnStyle:
        colour = data.get("colour")
        return cls(
            precision=data.get("precision"),
            log10=bool(data.get("log10", False)),
            colour=ColourRule.from_json(colour) if isinstance(colour, dict) else None,
        )


MAX_PRECISION = 12
DEFAULT_STYLE = ColumnStyle()


def format_fixed(value: float, decimals: int) -> str:
    """*value* with exactly *decimals* decimals; never ``-0.00``."""
    text = f"{value:.{decimals}f}"
    if text.startswith("-") and not text.strip("-0."):
        return text[1:]
    return text


def _fraction(value: float, lo: float, hi: float) -> float:
    if not math.isfinite(value) or hi <= lo:
        return 0.5
    return min(1.0, max(0.0, (value - lo) / (hi - lo)))


def _mix(a: int, b: int, t: float) -> int:
    out = 0
    for shift in (16, 8, 0):
        ca, cb = (a >> shift) & 255, (b >> shift) & 255
        out |= round(ca + (cb - ca) * t) << shift
    return out


def _hsl(hue: float, saturation: float, lightness: float) -> int:
    chroma = (1 - abs(2 * lightness - 1)) * saturation
    x = chroma * (1 - abs((hue / 60) % 2 - 1))
    m = lightness - chroma / 2
    r, g, b = [
        (chroma, x, 0),
        (x, chroma, 0),
        (0, chroma, x),
        (0, x, chroma),
        (x, 0, chroma),
        (chroma, 0, x),
    ][int(hue // 60) % 6]
    return (round((r + m) * 255) << 16) | (round((g + m) * 255) << 8) | round((b + m) * 255)
