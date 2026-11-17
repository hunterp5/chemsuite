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
"""Fingerprint similarity against several references (2D pharmacophores included: they are a
fingerprint kind), fused per row as the best or the mean over the references."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
from rdkit import Chem, rdBase

from chemsuite.analysis.fingerprints import (
    KINDS,
    Fingerprints,
    fingerprint_chunk,
    max_similarity,
    mean_similarity,
)
from chemsuite.core.jobs import CancelToken

from .model import Method, Reference


def reference_smiles(references: Sequence[Reference]) -> list[str]:
    """The references as SMILES without explicit hydrogens, so a 3D reference matches its 2D
    form (and a reference travels to a worker process as text)."""
    with rdBase.BlockLogs():
        return [Chem.MolToSmiles(Chem.RemoveHs(r.mol, sanitize=False)) for r in references]


def reference_fingerprints(references: Sequence[Reference], key: str) -> Fingerprints:
    """The references' fingerprints, as the searched molecules get theirs; a reference that
    fails is flagged."""
    packed, valid = fingerprint_chunk((key, reference_smiles(references)))
    return Fingerprints(KINDS[key], packed, valid)


def self_rows(references: Sequence[Reference], rows: int) -> np.ndarray:
    """Per searched row, the reference that is that row (-1: none), so it is left out."""
    out = np.full(rows, -1, dtype=np.int64)
    for index, reference in enumerate(references):
        if 0 <= reference.row < rows:
            out[reference.row] = index
    return out


def score(
    fps: Fingerprints,
    references: Sequence[Reference],
    method: Method,
    fuse: str,
    *,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(score, nearest reference)`` per row of *fps* (NaN and -1 for an invalid row). On
    cancel, raises ``Cancelled`` whose partial is that pair for the first rows, or None."""
    refs = reference_fingerprints(references, method.fingerprint)
    scan = max_similarity if fuse == "max" else mean_similarity
    return scan(fps, refs, method.metric, exclude_self=self_rows(references, len(fps)),
                alpha=method.alpha, beta=method.beta, progress=progress,
                cancel=cancel)  # fmt: skip
