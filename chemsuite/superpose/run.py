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
"""Superpose many rows, on the process pool, keeping finished rows on cancel.

Each row's stored conformers are all tried and the best fit is kept, moved into the
reference's frame; a row with only 2D gets conformers first (CONFORGE where CDPKit is
installed, else ETKDG).
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign

from chemsuite.conformers import conforge
from chemsuite.conformers.energies import delta_of
from chemsuite.conformers.etkdg import Cell, ConformerSettings, generate
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text, mol_to_blob, smiles_of

from . import anchor, field, o3a, shape, substructure
from .model import (
    COLOR,
    CORE,
    FIELD,
    MAX_RMSD,
    MCS,
    MEAN_RMSD,
    O3A,
    SHAPE,
    Fit,
    Settings,
    apply_transform,
    best,
    blob_with_reference,
    first_conformer,
    has_3d,
)

#: rows per worker task: a row takes milliseconds by shape, far more with conformers to make
CHUNK_ROWS = 4
ALIGNERS: dict[str, Callable[[Chem.Mol, Chem.Mol, Settings], list[Fit]]] = {
    MCS: substructure.align,
    CORE: substructure.align,
    O3A: o3a.align,
    SHAPE: shape.align,
    COLOR: shape.align,
    FIELD: field.align,
}
_FAILURES = (ValueError, RuntimeError, ArithmeticError)


def fits_of(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> list[Fit]:
    """Every probe conformer's fit: on the anchor when there is one and the probe has it."""
    if settings.anchored:
        fits = anchor.align(probe, reference, settings)
        if fits is not None:
            return fits
    return ALIGNERS[settings.method](probe, reference, settings)


def prepare_reference(mol: Chem.Mol) -> Chem.Mol:
    """The reference as the aligners take it: its first conformer, hydrogens added in place
    (a protein's ligand comes without them; MMFF typing and charges need them)."""
    if not has_3d(mol):
        raise ValueError("the reference has no 3D coordinates")
    with rdBase.BlockLogs():
        return Chem.AddHs(first_conformer(mol), addCoords=True)


def conformer_settings(count: int) -> ConformerSettings:
    """As Generate Conformers' default: CONFORGE where installed, else ETKDG with no window."""
    if conforge.available():
        return ConformerSettings(count=count, method="conforge")
    return ConformerSettings(count=count, energy_window=None)


def mol_of_cell(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def probe_of(mol: Chem.Mol, settings: Settings) -> Chem.Mol | None:
    """*mol* ready to fit: its stored 3D conformers (hydrogens added in place if it has none),
    or new conformers when it has only 2D (None if it would not embed)."""
    with rdBase.BlockLogs():
        if has_3d(mol):
            if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
                return mol
            return Chem.AddHs(mol, addCoords=True)
        return generate(mol, conformer_settings(settings.conformers)).mol


def superpose_one(mol: Chem.Mol, reference: Chem.Mol, settings: Settings) -> tuple[Chem.Mol, Fit]:
    """*mol*'s best-fitting conformer moved onto *reference* (prepared), with its fit."""
    probe = probe_of(mol, settings)
    if probe is None:
        raise ValueError("could not embed in 3D")
    found = best(fits_of(probe, reference, settings), settings.method)
    if found is None:
        raise ValueError("no fit")
    return placed(probe, found)


def placed(probe: Chem.Mol, fit: Fit) -> tuple[Chem.Mol, Fit]:
    """*probe*'s conformer that *fit* names, alone and moved into place, and the fit renumbered
    (``conf_id`` the conformer's position among the probe's: 0 is the first)."""
    out = Chem.Mol(probe)
    keep = Chem.Conformer(probe.GetConformer(fit.conf_id))
    out.RemoveAllConformers()
    keep.SetId(0)
    apply_transform(keep, fit.transform)
    out.AddConformer(keep, assignId=False)
    conf_index = [c.GetId() for c in probe.GetConformers()].index(fit.conf_id)
    return out, Fit(conf_index, fit.transform, fit.scores)


def fits_over_references(
    mol: Chem.Mol, references: Sequence[Chem.Mol], settings: Settings
) -> list[Fit | None] | None:
    """*mol*'s best fit onto each of *references* (each already prepared), in their order.

    The probe's conformers are generated once and tried against every reference, so scoring a
    row against several references costs one embedding, not one each. An entry is None where
    that reference gives no fit; the result itself is None when the molecule will not embed.
    """
    probe = probe_of(mol, settings)
    if probe is None:
        return None
    return probe_fits(probe, references, settings)


def probe_fits(
    probe: Chem.Mol, references: Sequence[Chem.Mol], settings: Settings
) -> list[Fit | None]:
    """*probe*'s best fit onto each of *references*, in their order (see :func:`probe_of`)."""
    out: list[Fit | None] = []
    for reference in references:
        try:
            out.append(best(fits_of(probe, reference, settings), settings.method))
        except _FAILURES:
            out.append(None)
    return out


def best_over_references(
    mol: Chem.Mol, references: Sequence[Chem.Mol], settings: Settings
) -> tuple[int, Fit] | None:
    """Which of *references* *mol* fits best, and that fit (None: none of them, or no embed)."""
    fits = fits_over_references(mol, references, settings)
    if fits is None:
        return None
    found = [(i, f) for i, f in enumerate(fits) if f is not None]
    if not found:
        return None
    return max(found, key=lambda pair: pair[1].key(settings.method))


def superpose_chunk(args: tuple[Settings, bytes, str, Sequence[Cell]]) -> list[tuple]:
    """``(blob, display, conformer number, scores, ΔE)`` per cell (a picklable worker entry
    point); a blank, invalid or failed cell gives ``(b"", "", nan, nans, nan)``. ΔE is the kept
    conformer's stored energy above its ensemble's lowest (NaN when it has none)."""
    settings, reference_blob, reference_id, cells = args
    reference = mol_from_blob(reference_blob)
    blank = tuple(np.nan for _ in settings.scores)
    out = []
    for cell in cells:
        mol = mol_of_cell(cell)
        if mol is None:
            out.append((b"", "", np.nan, blank, np.nan))
            continue
        try:
            moved, fit = superpose_one(mol, reference, settings)
        except _FAILURES:
            out.append((b"", "", np.nan, blank, np.nan))
            continue
        with rdBase.BlockLogs():
            display = smiles_of(Chem.RemoveHs(moved))
        scores = tuple(float(fit.scores.get(name, np.nan)) for name in settings.scores)
        out.append((blob_with_reference(moved, reference_id), display, fit.conf_id + 1.0, scores,
                    delta_of(moved)))  # fmt: skip
    return out


@dataclass
class SuperposeBatch:
    """Results cell for cell: blobs (``b""`` for none), display SMILES, the conformer used
    (1 = the first), one array per score (named as in ``settings.scores``) and the kept
    conformer's ΔE above its ensemble's lowest (kcal/mol; NaN where unknown)."""

    blobs: list[bytes]
    displays: list[str]
    conformers: np.ndarray
    scores: dict[str, np.ndarray]
    delta_e: np.ndarray

    @property
    def failed(self) -> int:
        return sum(1 for b in self.blobs if not b)


def _batch(rows: list[tuple], names: Sequence[str]) -> SuperposeBatch:
    blobs = [r[0] for r in rows]
    displays = [r[1] for r in rows]
    conformers = np.array([r[2] for r in rows], dtype=np.float64)
    table = np.array([r[3] for r in rows], dtype=np.float64).reshape(len(rows), len(names))
    delta_e = np.array([r[4] if len(r) > 4 else np.nan for r in rows], dtype=np.float64)
    return SuperposeBatch(blobs, displays, conformers,
                          {n: table[:, i] for i, n in enumerate(names)}, delta_e)  # fmt: skip


def _run(
    fn: Callable[[tuple], list[tuple]],
    chunks: Iterator[tuple],
    total: int,
    executor: Executor | None,
    progress: Callable[[int, int], None] | None,
    cancel: CancelToken | None,
) -> list[tuple]:
    rows: list[tuple] = []
    done = [0]

    def finished(part: list) -> None:
        done[0] += len(part)
        if progress is not None:
            progress(done[0], total)

    try:
        for part in ordered_map(fn, chunks, executor, cancel=cancel, finished=finished):
            rows.extend(part)
    except Cancelled:
        raise Cancelled(rows or None) from None
    return rows


def superpose_rows(
    cells: Sequence[Cell],
    reference: Chem.Mol,
    reference_id: str,
    settings: Settings,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> SuperposeBatch:
    """Every cell superposed onto *reference* (a molecule with 3D coordinates; its first
    conformer is used). Each result names *reference_id* (see
    :func:`~.model.reference_id_of`).

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` is the batch for the
    cells finished so far (the first ``len(batch.blobs)``), or None if none were.
    """
    prepared = mol_to_blob(prepare_reference(reference))
    chunks = ((settings, prepared, reference_id, list(cells[i : i + chunk_rows]))
              for i in range(0, len(cells), chunk_rows))  # fmt: skip
    try:
        rows = _run(superpose_chunk, chunks, len(cells), executor, progress, cancel)
    except Cancelled as stop:
        partial = stop.partial
        raise Cancelled(_batch(partial, settings.scores) if partial else None) from None
    return _batch(rows, settings.scores)


# -- each row's conformers onto its first ------------------------------------------------------
def conformers_onto_first(mol: Chem.Mol, heavy_only: bool = True) -> tuple[Chem.Mol, list[float]]:
    """A copy of *mol* with every conformer fitted onto the first, and each one's RMSD to it
    after the fit (the first's is 0)."""
    out = Chem.Mol(mol)
    if out.GetNumConformers() < 2:
        return out, [0.0] * out.GetNumConformers()
    heavy = [a.GetIdx() for a in out.GetAtoms() if a.GetAtomicNum() > 1]
    atoms = heavy if heavy_only and heavy else list(range(out.GetNumAtoms()))
    rms: list[float] = []
    rdMolAlign.AlignMolConformers(out, atomIds=atoms, RMSlist=rms)
    return out, [0.0, *rms]


def conformers_chunk(args: tuple[bool, Sequence[Cell]]) -> list[tuple]:
    """``(blob, display, mean RMSD, max RMSD)`` per cell; a cell without 3D gives blanks."""
    heavy_only, cells = args
    out = []
    for cell in cells:
        mol = mol_of_cell(cell)
        if not has_3d(mol):
            out.append((b"", "", np.nan, (np.nan, np.nan)))
            continue
        moved, rms = conformers_onto_first(mol, heavy_only)
        others = rms[1:] or [0.0]
        out.append((mol_to_blob(moved), smiles_of(Chem.RemoveHs(moved)), np.nan,
                    (float(np.mean(others)), float(np.max(others)))))  # fmt: skip
    return out


def align_conformers_rows(
    cells: Sequence[Cell],
    heavy_only: bool = True,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = 16,
) -> SuperposeBatch:
    """Each cell's conformers fitted onto its first; scores are the mean and the largest RMSD
    of the others (:data:`~.model.MEAN_RMSD`, :data:`~.model.MAX_RMSD`)."""
    names = (MEAN_RMSD, MAX_RMSD)
    chunks = ((heavy_only, list(cells[i : i + chunk_rows]))
              for i in range(0, len(cells), chunk_rows))  # fmt: skip
    try:
        rows = _run(conformers_chunk, chunks, len(cells), executor, progress, cancel)
    except Cancelled as stop:
        partial = stop.partial
        raise Cancelled(_batch(partial, names) if partial else None) from None
    return _batch(rows, names)
