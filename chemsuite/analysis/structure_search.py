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
"""Structure search: substructure, exact match and similarity over a structure column.

Substructure and exact searches screen first with RDKit pattern fingerprints, which can only
rule rows *out*: a row whose fingerprint lacks one of the query's bits cannot contain it. The
screen is a numpy bit test over the cached fingerprints (milliseconds for 200k rows). Only the
rows that pass are matched atom by atom, in chunks that a process pool can share.

Similarity searches compare fingerprints of the chosen kind and keep rows at or above the
threshold.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass
from typing import Literal

import numpy as np
from rdkit import Chem, rdBase

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_text

from .fingerprints import (
    DEFAULT_KIND,
    Cell,
    Fingerprints,
    cell_mol,
    mol_bits,
    pack,
    query_fingerprint,
    similarity,
)
from .fingerprints import KINDS as FINGERPRINT_KINDS

SearchMode = Literal["substructure", "exact", "similarity"]
MODES: dict[str, str] = {
    "substructure": "Substructure",
    "exact": "Exact match",
    "similarity": "Similarity",
}
SCREEN_KIND = "pattern"
CHUNK_ROWS = 1000
POOL_MIN_ROWS = 3000  # fewer candidates are matched faster in the job's own thread


class QueryError(ValueError):
    """The query text is not a structure (or not a SMARTS pattern, for substructures)."""


@dataclass(frozen=True)
class StructureQuery:
    """What to look for.

    ``text`` is SMILES, a molfile block or InChI; substructure searches also take SMARTS.
    ``stereo`` makes substructure and exact searches respect stereochemistry. ``threshold``,
    ``fingerprint`` and ``metric`` apply to similarity searches.
    """

    text: str
    mode: SearchMode = "substructure"
    stereo: bool = False
    threshold: float = 0.7
    fingerprint: str = DEFAULT_KIND
    metric: str = "tanimoto"

    @property
    def screen_kind(self) -> str:
        """The fingerprints the search needs for every row."""
        return self.fingerprint if self.mode == "similarity" else SCREEN_KIND

    def describe(self) -> str:
        """A short label, e.g. ``substructure c1ccccc1`` or ``similarity ≥ 0.70 to CCO``."""
        text = self.text.strip()
        if "\n" in text:
            text = "drawn structure"
        elif len(text) > 40:
            text = text[:39] + "…"
        if self.mode == "similarity":
            return f"similarity ≥ {self.threshold:.2f} to {text}"
        return f"{MODES[self.mode].lower()} {text}"


@dataclass(frozen=True)
class ParsedQuery:
    mol: Chem.Mol
    is_smarts: bool


def parse_query(query: StructureQuery) -> ParsedQuery:
    """The query molecule; raises :class:`QueryError` with a reason when there is none."""
    text = query.text
    if not text.strip():
        raise QueryError("Enter a structure.")
    mol = mol_from_text(text)
    if mol is not None:
        return ParsedQuery(mol, False)
    if query.mode == "substructure" and "\n" not in text.strip():
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmarts(text.strip())
        if mol is not None:
            return ParsedQuery(mol, True)
        raise QueryError("Not valid SMILES or SMARTS.")
    raise QueryError("Not a valid structure (SMARTS works only for substructure searches).")


def exact_key(mol: Chem.Mol, stereo: bool) -> str:
    """Canonical SMILES that two identical structures share (no explicit hydrogens)."""
    with rdBase.BlockLogs():
        try:
            return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=stereo)
        except (RuntimeError, ValueError):
            return ""


def query_screen(parsed: ParsedQuery) -> Fingerprints:
    """The query's pattern fingerprint, packed like the rows' fingerprints."""
    kind = FINGERPRINT_KINDS[SCREEN_KIND]
    with rdBase.BlockLogs():
        mol = parsed.mol
        if not parsed.is_smarts:
            mol = Chem.RemoveHs(mol)
        bits = mol_bits(mol, SCREEN_KIND)
    return Fingerprints(kind, pack(bits, kind), np.ones(1, dtype=np.bool_))


def screen(fps: Fingerprints, query: Fingerprints, *, exact: bool = False) -> np.ndarray:
    """``bool`` per row: rows that may contain the query (*exact*: may equal it)."""
    q = query.words[0]
    words = fps.words
    possible = np.all((words & q) == q, axis=1)
    if exact:
        possible &= fps.counts == query.counts[0]
    return possible & fps.valid


# -- matching (runs in worker processes) -------------------------------------------------------
def match_chunk(args: tuple[str, str, bool, Sequence[Cell]]) -> np.ndarray:
    """``bool`` per cell: does it match? A picklable worker entry point.

    The query travels as text and is parsed again here: query molecules from SMARTS do not
    survive pickling in every RDKit build.
    """
    text, mode, stereo, cells = args
    parsed = parse_query(StructureQuery(text, mode, stereo))  # type: ignore[arg-type]
    out = np.zeros(len(cells), dtype=np.bool_)
    with rdBase.BlockLogs():
        if mode == "exact":
            wanted = exact_key(parsed.mol, stereo)
            heavy = parsed.mol.GetNumHeavyAtoms()
            for i, cell in enumerate(cells):
                try:
                    mol = cell_mol(cell)
                except (RuntimeError, ValueError):
                    continue
                if mol is not None and mol.GetNumHeavyAtoms() == heavy:
                    out[i] = exact_key(mol, stereo) == wanted
            return out
        params = Chem.SubstructMatchParameters()
        params.useChirality = stereo
        params.maxMatches = 1
        for i, cell in enumerate(cells):
            try:
                mol = cell_mol(cell)
                out[i] = mol is not None and mol.HasSubstructMatch(parsed.mol, params)
            except (RuntimeError, ValueError):
                pass
    return out


def _chunks(query: StructureQuery, cells: Sequence[Cell], size: int) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        yield query.text, query.mode, query.stereo, list(cells[start : start + size])


@dataclass(frozen=True)
class SearchResult:
    hits: np.ndarray  # bool per searched row
    scores: np.ndarray | None = None  # similarity per row (NaN for invalid rows)
    candidates: int = 0  # rows left after the fingerprint screen

    @property
    def count(self) -> int:
        return int(self.hits.sum())


def search(
    query: StructureQuery,
    cells: Sequence[Cell] | Callable[[np.ndarray], Sequence[Cell]],
    fps: Fingerprints,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> SearchResult:
    """Search rows whose fingerprints (of ``query.screen_kind``) are *fps*.

    *cells* are the rows' structures in the same order, or a function that returns the cells at
    given row indices (so only screened candidates are read). *executor* is used only when
    enough candidates pass the screen to be worth the hand-off.
    """
    parsed = parse_query(query)
    n = len(fps)
    if query.mode == "similarity":
        scores = similarity(fps, query_fingerprint(parsed.mol, query.fingerprint), query.metric)
        hits = np.nan_to_num(scores, nan=-1.0) >= query.threshold
        return SearchResult(hits, scores, n)
    candidates = np.flatnonzero(screen(fps, query_screen(parsed), exact=query.mode == "exact"))
    hits = np.zeros(n, dtype=np.bool_)
    if not len(candidates):
        return SearchResult(hits, None, 0)
    cells_at = cells if callable(cells) else (lambda rows: [cells[i] for i in rows.tolist()])
    chosen = list(cells_at(candidates))
    use_pool = executor is not None and len(chosen) >= POOL_MIN_ROWS
    found: list[np.ndarray] = []
    done = 0
    for part in ordered_map(
        match_chunk,
        _chunks(query, chosen, chunk_rows),
        executor if use_pool else None,
        cancel=cancel,
    ):
        check(cancel)
        found.append(part)
        done += len(part)
        if progress is not None:
            progress(done, len(chosen))
    hits[candidates] = np.concatenate(found)
    return SearchResult(hits, None, len(candidates))


def match_atoms(mol: Chem.Mol, query: Chem.Mol, stereo: bool = False) -> tuple[list, list]:
    """``(atoms, bonds)`` of *mol* covered by any match of *query*, for highlighting."""
    params = Chem.SubstructMatchParameters()
    params.useChirality = stereo
    params.maxMatches = 100
    atoms: set[int] = set()
    bonds: set[int] = set()
    with rdBase.BlockLogs():
        try:
            matches = mol.GetSubstructMatches(query, params)
        except (RuntimeError, ValueError):
            return [], []
    for match in matches:
        atoms.update(match)
        for bond in query.GetBonds():
            b = mol.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()])
            if b is not None:
                bonds.add(b.GetIdx())
    return sorted(atoms), sorted(bonds)
