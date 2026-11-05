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
"""Conformers with RDKit's ETKDG, optimised with MMFF94 (or UFF where MMFF has no parameters).

Each molecule's conformers are stored together in one RDKit molecule, explicit hydrogens
included, lowest energy first, so a structure cell holds a whole conformer ensemble. Energies
are kcal/mol from the force field used; they compare conformers of one molecule only.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdDistGeom, rdMolAlign

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text, mol_to_blob, smiles_of

#: A structure cell: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
CHUNK_ROWS = 16  # molecules per worker task: each takes tens to hundreds of milliseconds
FORCE_FIELDS = ("mmff", "uff")
_FAILURES = (ValueError, RuntimeError, ArithmeticError)


@dataclass(frozen=True)
class ConformerSettings:
    """How to generate conformers. *energy_window* (kcal/mol above the lowest) drops strained
    conformers; ``None`` keeps them all. *optimize* is ``"mmff"``, ``"uff"`` or ``None``."""

    count: int = 10
    prune_rms: float = 0.5
    optimize: str | None = "mmff"
    max_iters: int = 500
    energy_window: float | None = 10.0
    seed: int = 42
    timeout: int = 60  # seconds of embedding per molecule before giving up

    def __post_init__(self) -> None:
        if self.count < 1:
            raise ValueError("count must be at least 1")
        if self.optimize not in (*FORCE_FIELDS, None):
            raise ValueError(f"unknown force field: {self.optimize!r}")


@dataclass
class Conformers:
    """One molecule's ensemble: *mol* holds the conformers in *energies* order (lowest first)."""

    mol: Chem.Mol | None
    energies: list[float] = field(default_factory=list)
    force_field: str = ""  # "MMFF94", "UFF" or "" (not optimised)
    error: str = ""

    @property
    def count(self) -> int:
        return self.mol.GetNumConformers() if self.mol is not None else 0


def _params(settings: ConformerSettings, random_coords: bool) -> rdDistGeom.EmbedParameters:
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = settings.seed
    params.pruneRmsThresh = settings.prune_rms
    params.numThreads = 1  # parallelism comes from worker processes
    params.useRandomCoords = random_coords
    params.timeout = settings.timeout
    return params


def _optimize(mol: Chem.Mol, settings: ConformerSettings) -> tuple[list[float], str]:
    """Optimise every conformer in place; ``(energies, force field name)``."""
    if settings.optimize is None:
        return [np.nan] * mol.GetNumConformers(), ""
    if settings.optimize == "mmff" and AllChem.MMFFHasAllMoleculeParams(mol):
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=settings.max_iters)
        name = "MMFF94"
    elif AllChem.UFFHasAllMoleculeParams(mol):
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=settings.max_iters)
        name = "UFF"
    else:
        return [np.nan] * mol.GetNumConformers(), ""
    return [float(energy) for _not_converged, energy in results], name


def _distinct(mol: Chem.Mol, order: list[int], prune_rms: float) -> list[int]:
    """*order* without conformers that optimised onto an earlier (lower-energy) one."""
    if prune_rms <= 0 or len(order) < 2:
        return order
    heavy = Chem.RemoveHs(mol)  # heavy-atom RMS with symmetry, on a copy (aligning moves atoms)
    ids = [conf.GetId() for conf in mol.GetConformers()]
    kept: list[int] = []
    for index in order:
        if all(
            rdMolAlign.GetBestRMS(heavy, heavy, ids[k], ids[index], maxMatches=100) > prune_rms
            for k in kept
        ):
            kept.append(index)
    return kept


def generate(mol: Chem.Mol, settings: ConformerSettings | None = None) -> Conformers:
    """Conformers of *mol* (which is not changed). Failures come back with ``mol=None``."""
    settings = settings or ConformerSettings()
    with rdBase.BlockLogs():
        try:
            work = Chem.AddHs(Chem.Mol(mol))
            work.RemoveAllConformers()
            ids = list(
                rdDistGeom.EmbedMultipleConfs(work, settings.count, _params(settings, False))
            )
            if not ids:  # hard cases (big rings, cages) often embed from random coordinates
                ids = list(
                    rdDistGeom.EmbedMultipleConfs(work, settings.count, _params(settings, True))
                )
            if not ids:
                return Conformers(None, error="could not embed in 3D")
            energies, force_field = _optimize(work, settings)
        except _FAILURES as err:
            return Conformers(None, error=str(err) or type(err).__name__)
        order = list(range(len(ids)))
        if force_field:
            order.sort(key=lambda i: energies[i])
            if settings.energy_window is not None:
                lowest = energies[order[0]]
                order = [i for i in order if energies[i] - lowest <= settings.energy_window]
            order = _distinct(work, order, settings.prune_rms)
        out = Chem.Mol(work)
        out.RemoveAllConformers()
        conformers = work.GetConformers()
        for new_id, index in enumerate(order):
            conformer = Chem.Conformer(conformers[index])
            conformer.SetId(new_id)
            out.AddConformer(conformer, assignId=False)
    return Conformers(out, [energies[i] for i in order], force_field)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def generate_chunk(args: tuple[ConformerSettings, Sequence[Cell]]) -> list[tuple]:
    """``(blob, display, count, energy span)`` per cell (a picklable worker entry point).

    A blank or invalid input, or a molecule that would not embed, gives ``(b"", "", 0, nan)``
    (count is NaN for a blank input, so failures stand out).
    """
    settings, cells = args
    out = []
    for cell in cells:
        mol = _mol(cell)
        if mol is None:
            out.append((b"", "", np.nan, np.nan))
            continue
        result = generate(mol, settings)
        if result.mol is None:
            out.append((b"", "", 0.0, np.nan))
            continue
        with rdBase.BlockLogs():
            display = smiles_of(Chem.RemoveHs(result.mol))
        energies = result.energies
        span = max(energies) - min(energies) if result.force_field else np.nan
        out.append((mol_to_blob(result.mol), display, float(result.count), float(span)))
    return out


def _chunks(settings: ConformerSettings, cells: Sequence[Cell], size: int) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        yield settings, list(cells[start : start + size])


@dataclass
class ConformerBatch:
    """Results cell for cell: blobs (``b""`` for none), display SMILES, counts, energy spans."""

    blobs: list[bytes]
    displays: list[str]
    counts: np.ndarray
    spans: np.ndarray

    @property
    def failed(self) -> int:
        """Structures that were given but would not embed."""
        return int(np.sum(self.counts == 0))


def generate_conformers(
    cells: Sequence[Cell],
    settings: ConformerSettings | None = None,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> ConformerBatch:
    settings = settings or ConformerSettings()
    total = len(cells)
    rows: list[tuple] = []
    for part in ordered_map(generate_chunk, _chunks(settings, cells, chunk_rows), executor):
        check(cancel)
        rows.extend(part)
        if progress is not None:
            progress(len(rows), total)
    if not rows:
        return ConformerBatch([], [], np.empty(0), np.empty(0))
    blobs, displays, counts, spans = zip(*rows, strict=True)
    return ConformerBatch(
        list(blobs),
        list(displays),
        np.array(counts, dtype=np.float64),
        np.array(spans, dtype=np.float64),
    )
