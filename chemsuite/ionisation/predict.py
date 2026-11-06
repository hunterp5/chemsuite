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

"""pKa prediction for many structures: chunked, optionally in worker processes."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

from rdkit import Chem, rdBase
from rdkit.Chem import Crippen

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text, mol_to_blob, smiles_of

from .molgpka import ACID
from .profile import Profile, profile

#: A structure cell: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
CHUNK_ROWS = 200  # per task in a worker process
INLINE_CHUNK_ROWS = 20  # in the job's own thread, so a cancel is noticed quickly

_FAILURES = (ValueError, RuntimeError, ArithmeticError, KeyError, IndexError)


@dataclass
class PkaResults:
    """One entry per row, in order. NaN (or "") where a row has no value."""

    most_acidic: list[float] = field(default_factory=list)
    most_basic: list[float] = field(default_factory=list)
    sites: list[str] = field(default_factory=list)  # "9.05 (base), 13.77 (acid)"
    charge: list[float] = field(default_factory=list)
    neutral: list[float] = field(default_factory=list)  # percent
    logd: list[float] = field(default_factory=list)
    species_blobs: list[bytes] = field(default_factory=list)  # b"" with a "" text: none
    species_texts: list[str] = field(default_factory=list)
    failed: int = 0

    def __len__(self) -> int:
        return len(self.charge)

    def extend(self, other: PkaResults) -> None:
        for name in (
            "most_acidic",
            "most_basic",
            "sites",
            "charge",
            "neutral",
            "logd",
            "species_blobs",
            "species_texts",
        ):
            getattr(self, name).extend(getattr(other, name))
        self.failed += other.failed


def sites_text(found: Profile) -> str:
    """Every site's pKa, lowest first, each marked acid or base."""
    ordered = sorted(found.sites, key=lambda site: site.pka)
    return ", ".join(
        f"{site.pka:.2f} ({'acid' if site.kind == ACID else 'base'})" for site in ordered
    )


def _none(value: float | None) -> float:
    return math.nan if value is None else value


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def predict_chunk(args: tuple[float, Sequence[Cell]]) -> PkaResults:
    """Results for one chunk of cells at pH *ph* (a picklable worker entry point)."""
    ph, cells = args
    out = PkaResults()
    with rdBase.BlockLogs():
        for cell in cells:
            found = None
            try:
                mol = _mol(cell)
                found = profile(mol) if mol is not None else None
                if found is not None:
                    species = found.major_species(ph)
                    row = (
                        _none(found.most_acidic),
                        _none(found.most_basic),
                        sites_text(found),
                        found.charge(ph),
                        100.0 * found.fraction_neutral(ph),
                        found.logd(float(Crippen.MolLogP(found.mol)), ph),
                        mol_to_blob(species),
                        smiles_of(species),
                    )
            except _FAILURES:
                found = None
            if found is None:
                out.failed += cell is not None
                row = (math.nan, math.nan, "", math.nan, math.nan, math.nan, b"", "")
            out.most_acidic.append(row[0])
            out.most_basic.append(row[1])
            out.sites.append(row[2])
            out.charge.append(row[3])
            out.neutral.append(row[4])
            out.logd.append(row[5])
            out.species_blobs.append(row[6])
            out.species_texts.append(row[7])
    return out


def predict_pka(
    cells: Sequence[Cell],
    ph: float = 7.4,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> PkaResults:
    """MolGpKa's pKa values for each cell, and what they mean at pH *ph*.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` holds the results for the
    rows finished so far (the first rows, in order), or None if none were.
    """
    total = len(cells)
    size = CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS
    chunks = ((ph, list(cells[i : i + size])) for i in range(0, total, size))
    results = PkaResults()
    try:
        for part in ordered_map(predict_chunk, chunks, executor, cancel=cancel):
            results.extend(part)
            if progress is not None:
                progress(len(results), total)
    except Cancelled:
        raise Cancelled(results if len(results) else None) from None
    return results
