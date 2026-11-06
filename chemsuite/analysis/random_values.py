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

"""Random numbers for a column, and random samples of rows (Data ▸ Utilities ▸ Random).

Every draw takes an optional seed: the same seed and settings give the same numbers or rows
again, so a split into training and test sets can be reproduced. Without a seed each draw is
different. numpy's PCG64 generator is used throughout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

UNIFORM, INTEGER, NORMAL, LOGNORMAL, ORDER = "uniform", "integer", "normal", "lognormal", "order"
DISTRIBUTIONS = {
    UNIFORM: "Uniform (any value from low to high)",
    INTEGER: "Whole numbers (low to high, both included)",
    NORMAL: "Normal (mean and standard deviation)",
    LOGNORMAL: "Log-normal (median and spread factor)",
    ORDER: "Random order (1 to the number of rows)",
}


@dataclass(frozen=True)
class RandomSpec:
    distribution: str = UNIFORM
    low: float = 0.0
    high: float = 1.0
    mean: float = 0.0  # normal: mean; log-normal: median
    sd: float = 1.0  # normal: standard deviation; log-normal: spread factor (> 1)
    clip: bool = False  # normal / log-normal: keep within low…high
    decimals: int | None = 3  # None: full precision
    seed: int | None = None

    def problem(self) -> str:
        """Why these settings cannot be drawn ("" when they can)."""
        if self.distribution not in DISTRIBUTIONS:
            return f"Unknown distribution {self.distribution!r}."
        ranged = self.distribution in (UNIFORM, INTEGER) or self.clip
        if ranged and self.low > self.high:
            return "Low must not be above high."
        if self.distribution == INTEGER and np.floor(self.high) < np.ceil(self.low):
            return "There is no whole number between low and high."
        if self.distribution == NORMAL and self.sd < 0:
            return "The standard deviation cannot be negative."
        if self.distribution == LOGNORMAL:
            if self.mean <= 0:
                return "The median of a log-normal distribution must be above 0."
            if self.sd < 1:
                return "The spread factor must be at least 1 (1 means no spread)."
        return ""


def random_numbers(spec: RandomSpec, n: int) -> np.ndarray:
    """*n* numbers drawn as *spec* says (float64; whole numbers for INTEGER and ORDER)."""
    problem = spec.problem()
    if problem:
        raise ValueError(problem)
    rng = np.random.default_rng(spec.seed)
    kind = spec.distribution
    if kind == ORDER:
        return rng.permutation(n).astype(np.float64) + 1
    if kind == INTEGER:
        lo, hi = int(np.ceil(spec.low)), int(np.floor(spec.high))
        return rng.integers(lo, hi, size=n, endpoint=True).astype(np.float64)
    if kind == UNIFORM:
        values = rng.uniform(spec.low, spec.high, size=n)
    elif kind == NORMAL:
        values = rng.normal(spec.mean, spec.sd, size=n)
    else:  # log-normal: median m and spread factor s, i.e. log(x) ~ N(log m, log s)
        values = rng.lognormal(np.log(spec.mean), np.log(spec.sd), size=n)
    if spec.clip and kind in (NORMAL, LOGNORMAL):
        values = np.clip(values, spec.low, spec.high)
    if spec.decimals is not None:
        values = np.round(values, spec.decimals)
        if kind == UNIFORM:  # rounding must not step outside the range
            values = np.clip(values, spec.low, spec.high)
    return values


def random_sample(oids: np.ndarray, count: int, seed: int | None = None) -> np.ndarray:
    """*count* of *oids* chosen at random without repeats (all of them if fewer), ascending."""
    oids = np.asarray(oids, dtype=np.int64)
    count = max(0, min(int(count), len(oids)))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(oids, size=count, replace=False))


def sample_size(total: int, amount: float, percent: bool) -> int:
    """How many rows *amount* means: a count, or a percentage of *total* (rounded)."""
    if percent:
        return max(0, min(total, int(round(total * amount / 100))))
    return max(0, min(total, int(amount)))
