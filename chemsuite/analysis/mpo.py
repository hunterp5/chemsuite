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
"""Multi-parameter optimisation (MPO): desirability scores combined into one score.

Each property maps to a desirability between 0 (unwanted) and 1 (ideal); the overall score is
their weighted geometric mean (Derringer-Suich: one unacceptable property sinks the compound)
or weighted arithmetic mean (properties trade off against each other).
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass

import numpy as np

SHAPES = {
    "maximize": "Higher is better",
    "minimize": "Lower is better",
    "target": "Within a range, best at a target",
    "gaussian": "Near a centre (bell curve)",
}
COMBINE = {"geometric": "Geometric mean", "arithmetic": "Arithmetic mean"}


@dataclass(frozen=True)
class Desirability:
    """How one column scores.

    * ``maximize``: 0 at or below *low*, rising linearly to 1 at *high*.
    * ``minimize``: 1 at or below *low*, falling linearly to 0 at *high*.
    * ``target``: 0 outside [*low*, *high*], 1 at *target* (the midpoint by default).
    * ``gaussian``: exp(-z²/2) with z = (x − *target*) / *width*.

    *hard* makes maximize/minimize/target a step: 1 inside the acceptable side, else 0.
    """

    column: str
    shape: str = "maximize"
    low: float = 0.0
    high: float = 1.0
    target: float | None = None
    width: float = 1.0
    weight: float = 1.0
    hard: bool = False

    def __post_init__(self) -> None:
        if self.shape not in SHAPES:
            raise ValueError(f"unknown desirability shape: {self.shape!r}")
        if self.weight < 0 or not math.isfinite(self.weight):
            raise ValueError("weights must be zero or more")

    def describe(self) -> str:
        if self.shape == "gaussian":
            return f"best at {self.target or 0:g}, width {self.width:g}"
        lo, hi = sorted((self.low, self.high))
        if self.shape == "target":
            t = self.mid()
            return (
                f"0 outside {lo:g}–{hi:g}, 1 at {t:g}"
                if not self.hard
                else f"1 within {lo:g}–{hi:g}"
            )
        if self.hard:
            return f"1 at ≥ {hi:g}" if self.shape == "maximize" else f"1 at ≤ {lo:g}"
        better = "rising" if self.shape == "maximize" else "falling"
        return f"{better} from {lo:g} to {hi:g}"

    def mid(self) -> float:
        lo, hi = sorted((self.low, self.high))
        t = 0.5 * (lo + hi) if self.target is None else self.target
        return min(max(t, lo), hi)

    def score(self, values: np.ndarray) -> np.ndarray:
        """Desirability per value in [0, 1]; NaN where the value is missing."""
        x = np.asarray(values, dtype=np.float64)
        lo, hi = sorted((self.low, self.high))
        with np.errstate(invalid="ignore", divide="ignore"):
            if self.shape == "gaussian":
                width = abs(self.width)
                center = self.target or 0.0
                d = np.exp(-0.5 * ((x - center) / width) ** 2) if width > 0 else (x == center) * 1.0
            elif self.shape == "maximize":
                d = (x >= hi) * 1.0 if self.hard or hi == lo else (x - lo) / (hi - lo)
            elif self.shape == "minimize":
                d = (x <= lo) * 1.0 if self.hard or hi == lo else (hi - x) / (hi - lo)
            elif self.hard:
                d = ((x >= lo) & (x <= hi)) * 1.0
            else:
                t = self.mid()
                rising = (x - lo) / (t - lo) if t > lo else (x >= t) * 1.0
                falling = (hi - x) / (hi - t) if hi > t else (x <= t) * 1.0
                d = np.where(x <= t, rising, falling)
                d = np.where((x < lo) | (x > hi), 0.0, d)
            d = np.clip(d, 0.0, 1.0)
        return np.where(np.isfinite(x), d, np.nan)


def combine(
    scores: Sequence[np.ndarray],
    weights: Sequence[float],
    method: str = "geometric",
    *,
    skip_missing: bool = False,
) -> np.ndarray:
    """Weighted mean of desirabilities per row.

    A row missing any property is NaN, unless *skip_missing*: then it is scored on the
    properties it has (NaN only when it has none).
    """
    if method not in COMBINE:
        raise ValueError(f"unknown combine method: {method!r}")
    if not scores:
        raise ValueError("no properties to combine")
    d = np.column_stack([np.asarray(s, dtype=np.float64) for s in scores])
    w = np.asarray(weights, dtype=np.float64)[None, :] * np.ones_like(d)
    missing = ~np.isfinite(d)
    w = np.where(missing, 0.0, w)
    total = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        if method == "arithmetic":
            out = np.where(missing, 0.0, d * w).sum(axis=1) / total
        else:
            logs = np.log(np.where(missing, 1.0, d))  # log(0) = -inf gives 0, as it should
            out = np.exp(np.where(w > 0, logs * w, 0.0).sum(axis=1) / total)
    out = np.where(total > 0, out, np.nan)
    if not skip_missing:
        out = np.where(missing.any(axis=1), np.nan, out)
    return np.clip(out, 0.0, 1.0)


@dataclass(frozen=True)
class MpoProfile:
    """A complete scoring scheme, saved to and loaded from JSON."""

    properties: tuple[Desirability, ...]
    method: str = "geometric"
    skip_missing: bool = False

    def score(self, columns: dict[str, np.ndarray]) -> tuple[np.ndarray, list[np.ndarray]]:
        """``(overall score, desirability per property)``; *columns* maps name → values."""
        parts = [p.score(columns[p.column]) for p in self.properties]
        overall = combine(
            parts, [p.weight for p in self.properties], self.method, skip_missing=self.skip_missing
        )
        return overall, parts

    def to_json(self) -> str:
        return json.dumps(
            {
                "format": "chemsuite-mpo",
                "version": 1,
                "method": self.method,
                "skip_missing": self.skip_missing,
                "properties": [asdict(p) for p in self.properties],
            },
            indent=2,
        )

    @classmethod
    def from_json(cls, text: str) -> MpoProfile:
        try:
            data = json.loads(text)
            if data.get("format") != "chemsuite-mpo":
                raise ValueError("not a chemsuite MPO profile")
            fields = Desirability.__dataclass_fields__
            properties = tuple(
                Desirability(**{k: v for k, v in item.items() if k in fields})
                for item in data["properties"]
            )
            return cls(properties, data.get("method", "geometric"), bool(data.get("skip_missing")))
        except (json.JSONDecodeError, KeyError, TypeError, AttributeError) as err:
            raise ValueError(f"not a valid MPO profile: {err}") from err
