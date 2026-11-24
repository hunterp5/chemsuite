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
"""Filling the gaps and putting the columns on one scale, as arrays rather than estimators.

A descriptor fails on some molecules, so the pool has holes; a distance-based or penalized
model needs columns of comparable size. Both are learned from a training fold and replayed as
arithmetic, which is what lets :mod:`chemsuite.qsar.apply` work on numpy alone — scikit-learn's
``SimpleImputer`` and ``StandardScaler`` would have to be pickled to travel.

The order a model replays is always: **fill** every column of the pool, **take** the columns
selection kept, then **scale** those. :class:`Prepare` holds both halves.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: A column whose spread is below this is treated as flat: dividing by it would blow small
#: numerical noise up into a feature.
MIN_SCALE = 1e-12


@dataclass(frozen=True)
class Prepare:
    """What a fold learned about its columns, as plain arrays.

    *fill* has one value per column of the whole pool; *center* and *scale* have one per
    **kept** column, or are None when the candidate wants the columns as they are (a tree does
    not care, and scaling fingerprint counts only costs accuracy in the archive).
    """

    fill: np.ndarray
    center: np.ndarray | None = None
    scale: np.ndarray | None = None

    @classmethod
    def fit(cls, x: np.ndarray) -> Prepare:
        """The fill values of a training fold: each column's median, and 0 for a column with
        no value at all (it carries nothing, and selection drops it)."""
        x = np.asarray(x, dtype=np.float64)
        usable = np.where(np.isfinite(x), x, np.nan)
        fill = np.zeros(x.shape[1], dtype=np.float64)
        has_value = np.isfinite(usable).any(axis=0)  # nanmedian warns on an all-NaN column
        if has_value.any():
            fill[has_value] = np.nanmedian(usable[:, has_value], axis=0)
        return cls(fill)

    def filled(self, x: np.ndarray) -> np.ndarray:
        """*x* with every hole — NaN or an infinity — filled from the training fold."""
        x = np.asarray(x, dtype=np.float64)
        if x.shape[1] != len(self.fill):
            raise ValueError(f"{x.shape[1]} columns, but the fit saw {len(self.fill)}")
        out = np.where(np.isfinite(x), x, np.broadcast_to(self.fill, x.shape))
        return np.ascontiguousarray(out)

    def with_scaling(self, x: np.ndarray, kept: np.ndarray, scaled: bool) -> Prepare:
        """This fit, plus the center and spread of the kept columns of *x* (already filled).

        *scaled* False keeps the columns as they are, so the archive stores no scaling.
        """
        if not scaled:
            return Prepare(self.fill)
        taken = np.asarray(x, dtype=np.float64)[:, kept]
        center = taken.mean(axis=0)
        spread = taken.std(axis=0)
        return Prepare(self.fill, center, np.where(spread > MIN_SCALE, spread, 1.0))

    def design(self, x: np.ndarray, kept: np.ndarray) -> np.ndarray:
        """Raw pool columns in, the matrix the estimator sees out: fill, take, scale."""
        taken = self.filled(x)[:, np.asarray(kept, dtype=np.int64)]
        if self.center is None or self.scale is None:
            return taken
        return (taken - self.center) / self.scale

    def as_arrays(self) -> dict[str, np.ndarray]:
        """For the model archive: only the parts that exist."""
        out = {"prepare_fill": self.fill}
        if self.center is not None and self.scale is not None:
            out["prepare_center"] = self.center
            out["prepare_scale"] = self.scale
        return out

    @classmethod
    def from_arrays(cls, arrays: dict[str, np.ndarray]) -> Prepare:
        if "prepare_fill" not in arrays:
            raise ValueError("the model archive has no fill values")
        center = arrays.get("prepare_center")
        scale = arrays.get("prepare_scale")
        return cls(np.asarray(arrays["prepare_fill"], dtype=np.float64), center, scale)
