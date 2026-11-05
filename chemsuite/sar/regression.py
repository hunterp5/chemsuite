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
"""Ordinary least squares with the statistics a medicinal chemist expects.

Shared by Free-Wilson and Hansch. Columns that the data cannot separate (a substituent that
always occurs with another) are reported, not silently given arbitrary values.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Fit:
    names: tuple[str, ...]  # one per column of X (the intercept is separate)
    coefficients: np.ndarray  # NaN for a column that could not be estimated
    errors: np.ndarray  # standard errors (NaN when not estimable or no residual freedom)
    intercept: float
    intercept_error: float
    predicted: np.ndarray
    residuals: np.ndarray
    n: int
    parameters: int  # estimated parameters, intercept included
    r2: float
    q2: float  # leave-one-out cross-validated R²
    rmse: float
    s: float  # residual standard deviation (n − parameters degrees of freedom)
    f: float

    @property
    def estimable(self) -> np.ndarray:
        return np.isfinite(self.coefficients)


def _estimable_columns(x: np.ndarray) -> list[int]:
    """Greedy: keep a column when it adds rank to the intercept and the columns kept so far."""
    kept: list[int] = []
    basis = np.ones((len(x), 1))
    rank = 1
    for j in range(x.shape[1]):
        trial = np.column_stack([basis, x[:, j]])
        if np.linalg.matrix_rank(trial) > rank:
            basis, rank = trial, rank + 1
            kept.append(j)
    return kept


def least_squares(x: np.ndarray, y: np.ndarray, names: Sequence[str]) -> Fit:
    x = np.asarray(x, dtype=np.float64).reshape(len(y), -1)
    y = np.asarray(y, dtype=np.float64)
    n = len(y)
    kept = _estimable_columns(x)
    design = np.column_stack([np.ones(n), x[:, kept]])
    p = design.shape[1]
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    predicted = design @ beta
    residuals = y - predicted
    rss = float(residuals @ residuals)
    tss = float(((y - y.mean()) ** 2).sum())
    dof = n - p
    s = math.sqrt(rss / dof) if dof > 0 else math.nan
    xtx_inv = np.linalg.pinv(design.T @ design)
    errors_all = np.sqrt(np.diag(xtx_inv)) * s if dof > 0 else np.full(p, math.nan)
    hat = np.einsum("ij,jk,ik->i", design, xtx_inv, design)
    with np.errstate(divide="ignore", invalid="ignore"):
        loo = residuals / (1 - hat)
    press = float((loo**2).sum()) if np.all(hat < 1 - 1e-10) else math.nan
    coefficients = np.full(x.shape[1], math.nan)
    errors = np.full(x.shape[1], math.nan)
    coefficients[kept] = beta[1:]
    errors[kept] = errors_all[1:]
    r2 = 1 - rss / tss if tss > 0 else math.nan
    f = ((tss - rss) / (p - 1)) / (rss / dof) if p > 1 and dof > 0 and rss > 0 else math.nan
    return Fit(
        tuple(names),
        coefficients,
        errors,
        float(beta[0]),
        float(errors_all[0]),
        predicted,
        residuals,
        n,
        p,
        r2,
        1 - press / tss if tss > 0 and math.isfinite(press) else math.nan,
        math.sqrt(rss / n) if n else math.nan,
        s,
        f,
    )
