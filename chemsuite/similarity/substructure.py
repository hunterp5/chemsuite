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
"""Similarity by common substructure: the maximum common substructure's share of the two
molecules' heavy atoms (an atom Tanimoto), or whether two molecules share their Bemis-Murcko
scaffold (exactly, or as a generic scaffold: every atom carbon, every bond single)."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric

from chemsuite.analysis.fingerprints import Cell, cell_mol
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.descriptors.scaffold import murcko

from .fusion import fuse
from .model import GENERIC_SCAFFOLD, MCS, SCAFFOLD, Method, Reference

#: per-row work fails with these (a molecule RDKit cannot handle), never with anything else
_FAILURES = (ValueError, RuntimeError)
NAN = float("nan")


def mcs_tanimoto(a: Chem.Mol, b: Chem.Mol, timeout: int) -> float:
    """Atoms in the largest common substructure over the atoms in either molecule (heavy atoms;
    rings match only rings and only whole). Past *timeout* seconds, the largest found so far."""
    params = rdFMCS.MCSParameters()
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareOrder
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    params.Timeout = timeout
    common = rdFMCS.FindMCS([a, b], params).numAtoms
    total = a.GetNumAtoms() + b.GetNumAtoms() - common
    return common / total if total else NAN


def scaffold_key(mol: Chem.Mol, generic: bool) -> str:
    """Canonical SMILES of *mol*'s scaffold ("" when it has no ring)."""
    scaffold = murcko(mol)
    if not scaffold.GetNumAtoms():
        return ""
    if generic:
        scaffold = MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(scaffold)


def score_chunk(
    args: tuple[str, int, Sequence[str], Sequence[int], Sequence[Cell]],
) -> list[tuple[float, ...]]:
    """Per cell, its score against each reference (a picklable worker entry point).

    *args* is ``(method, MCS timeout, reference SMILES, the reference each cell is (-1: none),
    cells)``. NaN for a blank or invalid cell, for its own reference, and where RDKit fails.
    """
    kind, timeout, ref_smiles, own, cells = args
    if kind not in (MCS, SCAFFOLD, GENERIC_SCAFFOLD):
        raise ValueError(f"not a substructure method: {kind!r}")
    blank = tuple(NAN for _ in ref_smiles)
    with rdBase.BlockLogs():
        refs = [cell_mol(s) for s in ref_smiles]
        keys = None
        if kind in (SCAFFOLD, GENERIC_SCAFFOLD):
            generic = kind == GENERIC_SCAFFOLD
            keys = [scaffold_key(r, generic) if r is not None else None for r in refs]
        out: list[tuple[float, ...]] = []
        for cell, skip in zip(cells, own, strict=True):
            mol = cell_mol(cell)
            if mol is None:
                out.append(blank)
                continue
            try:
                if keys is not None:
                    mine = scaffold_key(mol, kind == GENERIC_SCAFFOLD)
                    row = [NAN if k is None else float(bool(mine) and mine == k) for k in keys]
                else:
                    row = [NAN if r is None else mcs_tanimoto(mol, r, timeout) for r in refs]
            except _FAILURES:
                out.append(blank)
                continue
            if skip >= 0:
                row[skip] = NAN
            out.append(tuple(row))
    return out


#: cells per worker task: the common substructure of a pair takes far longer than a scaffold
CHUNK_ROWS = {MCS: 4, SCAFFOLD: 500, GENERIC_SCAFFOLD: 500}


def _chunks(
    kind: str,
    timeout: int,
    smiles: Sequence[str],
    own: np.ndarray,
    cells: Sequence[Cell],
    size: int,
) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        stop = start + size
        yield kind, timeout, list(smiles), own[start:stop].tolist(), list(cells[start:stop])


def score(
    cells: Sequence[Cell],
    references: Sequence[Reference],
    method: Method,
    fuse_how: str,
    own: np.ndarray,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(score, nearest reference)`` per cell by *method* (a substructure method).

    *own* gives, per cell, the reference that is that cell (-1: none), so it is not compared
    with itself. On cancel, raises ``Cancelled`` whose partial is that pair for the cells
    finished so far, or None if none were.
    """
    from .fingerprint import reference_smiles

    size = CHUNK_ROWS[method.kind]
    chunks = _chunks(
        method.kind, method.mcs_timeout, reference_smiles(references), own, cells, size
    )
    rows: list[tuple[float, ...]] = []
    done = [0]

    def finished(part: list) -> None:
        done[0] += len(part)
        if progress is not None:
            progress(done[0], len(cells))

    def fused(found: list[tuple[float, ...]]) -> tuple[np.ndarray, np.ndarray]:
        matrix = np.array(found, dtype=np.float64).reshape(len(found), len(references))
        return fuse(matrix, fuse_how)

    try:
        for part in ordered_map(score_chunk, chunks, executor, cancel=cancel, finished=finished):
            rows.extend(part)
    except Cancelled:
        raise Cancelled(fused(rows) if rows else None) from None
    return fused(rows)
