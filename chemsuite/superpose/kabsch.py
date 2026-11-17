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
"""Least-squares rigid fits of paired points (Kabsch, Acta Cryst. A 1976, 32, 922)."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

#: symmetric matches tried per fit (a core with many symmetric groups has very many)
MAX_MATCHES = 1000


def fit(probe: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, float]:
    """The 4×4 rigid transform that best moves the *probe* points onto the *reference* points
    (same length, paired in order), and the RMSD after it."""
    p_mean = probe.mean(axis=0)
    q_mean = reference.mean(axis=0)
    p = probe - p_mean
    q = reference - q_mean
    u, _s, vt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    transform = np.eye(4)
    transform[:3, :3] = rotation
    transform[:3, 3] = q_mean - rotation @ p_mean
    moved = p @ rotation.T
    rmsd = float(np.sqrt(((moved - q) ** 2).sum(axis=1).mean()))
    return transform, rmsd


def best_fit(
    probe: np.ndarray, reference: np.ndarray, probe_matches: Sequence[Sequence[int]],
    reference_match: Sequence[int],
) -> tuple[np.ndarray, float]:  # fmt: skip
    """The fit over the probe match (of several symmetric ones) that gives the lowest RMSD:
    ``probe[match[k]]`` pairs with ``reference[reference_match[k]]``."""
    target = reference[list(reference_match)]
    found: tuple[np.ndarray, float] | None = None
    for match in probe_matches[:MAX_MATCHES]:
        candidate = fit(probe[list(match)], target)
        if found is None or candidate[1] < found[1]:
            found = candidate
    if found is None:
        raise ValueError("no atom pairs to fit")
    return found


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))
