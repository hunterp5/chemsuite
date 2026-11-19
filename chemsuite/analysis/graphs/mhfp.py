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
"""MHFP6: the MinHash fingerprint TMAP was built on.

Probst, D. and Reymond, J.-L. (2018). "A probabilistic molecular fingerprint for big data
settings." *J. Cheminf.* 10, 66. Each molecule's circular substructures up to radius 3 (as
SMILES strings, its "shingles") are hashed; for each of *perms* random hash functions the
fingerprint keeps the smallest value over the shingles. Two molecules agree at a position with
probability equal to the Jaccard similarity of their shingle sets, so the share of agreeing
positions estimates it, on the same 0–1 scale as Tanimoto.

The values come from RDKit's own encoder (``rdMHFPFingerprint``); unlike the bit fingerprints
they are not folded, so :mod:`.knn` puts them straight into its prefix forest.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass

import numpy as np
from rdkit import rdBase
from rdkit.Chem import rdMHFPFingerprint

from chemsuite.analysis.fingerprints import Cell, cell_mol
from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map

#: hash functions per molecule (each a ``uint32``)
PERMS = 256
RADIUS = 3  # MHFP6: substructures up to six bonds across
SEED = 42  # the encoder's: the same molecule always gives the same values
CHUNK_ROWS = 2000


@dataclass(frozen=True)
class MinHashes:
    """One MHFP per row (``n × perms`` ``uint32``), with which rows had a valid structure."""

    hashes: np.ndarray
    valid: np.ndarray

    def __len__(self) -> int:
        return len(self.valid)

    def take(self, rows: np.ndarray) -> MinHashes:
        return MinHashes(self.hashes[rows], self.valid[rows])


def agreement(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The share of positions where MinHash rows agree (broadcast over leading axes): the
    estimated Jaccard similarity of the molecules' shingle sets."""
    return (a == b).mean(axis=-1)


def mhfp_chunk(args: tuple[int, Sequence[Cell]]) -> tuple[np.ndarray, np.ndarray]:
    """``(hashes, valid flags)`` for one chunk (a picklable worker entry point)."""
    perms, cells = args
    encoder = rdMHFPFingerprint.MHFPEncoder(perms, SEED)
    out = np.zeros((len(cells), perms), dtype=np.uint32)
    valid = np.zeros(len(cells), dtype=np.bool_)
    with rdBase.BlockLogs():
        for i, cell in enumerate(cells):
            try:
                mol = cell_mol(cell)
                if mol is None or not mol.GetNumAtoms():
                    continue
                out[i] = np.asarray(encoder.EncodeMol(mol, radius=RADIUS), dtype=np.uint32)
                valid[i] = True
            except (ValueError, RuntimeError):
                pass
    return out, valid


def _chunks(perms: int, cells: Sequence[Cell], size: int) -> Iterator[tuple[int, list[Cell]]]:
    for start in range(0, len(cells), size):
        yield perms, list(cells[start : start + size])


def compute_mhfp(
    cells: Sequence[Cell],
    perms: int = PERMS,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> MinHashes:
    """The MHFP6 of every cell (blank and invalid cells are not valid)."""
    parts: list[np.ndarray] = []
    flags: list[np.ndarray] = []
    done = 0
    for hashes, valid in ordered_map(
        mhfp_chunk, _chunks(perms, cells, chunk_rows), executor, cancel=cancel
    ):
        check(cancel)
        parts.append(hashes)
        flags.append(valid)
        done += len(valid)
        if progress is not None:
            progress(done, len(cells))
    if not parts:
        return MinHashes(np.zeros((0, perms), np.uint32), np.zeros(0, bool))
    return MinHashes(np.concatenate(parts), np.concatenate(flags))
