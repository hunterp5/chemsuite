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
"""Selectivity of a primary assay against off-target assays, as log units (1 = 10-fold)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

_ASSAY_HINTS = ("pic50", "pec50", "pki", "pkd", "ic50", "ec50", "ki", "kd", "activity", "potency")


def looks_like_assay(name: str) -> bool:
    """A column name that looks like an activity readout (pIC50, Ki…), to preselect it."""
    flat = name.casefold().replace(" ", "").replace("_", "").replace("-", "")
    words = name.casefold().replace("-", " ").replace("_", " ").split()
    return any(key in flat for key in _ASSAY_HINTS[:6]) or any(w in ("ki", "kd") for w in words)


def selectivity_name(off_target: str) -> str:
    return f"logSel vs {off_target}"


MIN_NAME = "Min logSel"
WINDOW_NAME = "Narrow window"


@dataclass(frozen=True)
class SelectivityResult:
    by_target: dict[str, np.ndarray]  # off-target name → log selectivity per row
    minimum: np.ndarray  # the smallest measured log selectivity per row (NaN: none)
    narrow: list[str]  # "yes", "no" or "" per row


def selectivity(
    primary: np.ndarray,
    off_targets: dict[str, np.ndarray],
    *,
    concentrations: bool = False,
    potency_cutoff: float | None = None,
    window: float = 1.0,
) -> SelectivityResult:
    """Log selectivity of *primary* over each off-target.

    With p-values (pIC50, pKi: higher is more potent) it is primary − off-target. With
    *concentrations* (IC50, Ki in one unit: lower is more potent) it is
    log10(off-target / primary). Non-positive concentrations count as missing.

    A row is in a "narrow window" when its primary potency reaches *potency_cutoff* (a
    p-value, or a concentration at most that when *concentrations*) and its smallest
    selectivity is below *window*. No cutoff: every potent-enough row qualifies.
    """

    def p_scale(v: np.ndarray) -> np.ndarray:
        v = np.asarray(v, dtype=np.float64)
        if not concentrations:
            return v
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(v > 0, -np.log10(np.where(v > 0, v, 1.0)), np.nan)

    p = p_scale(primary)
    by_target = {name: p - p_scale(values) for name, values in off_targets.items()}
    stacked = np.column_stack(list(by_target.values())) if by_target else np.empty((len(p), 0))
    measured = np.isfinite(stacked).any(axis=1) if stacked.size else np.zeros(len(p), bool)
    with np.errstate(invalid="ignore"):
        minimum = np.where(
            measured,
            np.nanmin(np.where(np.isfinite(stacked), stacked, np.inf), axis=1)
            if stacked.size
            else np.nan,
            np.nan,
        )
        cutoff = (
            np.ones(len(p), bool)
            if potency_cutoff is None
            else p >= (-np.log10(potency_cutoff) if concentrations else potency_cutoff)
        )
    narrow = [
        "" if not judged else ("yes" if potent and low < window else "no")
        for judged, potent, low in zip(measured & np.isfinite(p), cutoff, minimum, strict=True)
    ]
    return SelectivityResult(by_target, minimum, narrow)


def off_target_candidates(names: Sequence[str], primary: str) -> list[str]:
    return [n for n in names if n != primary and looks_like_assay(n)]
