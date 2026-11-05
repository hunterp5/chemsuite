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
"""Fingerprints as packed bit matrices, with fast similarity arithmetic.

One fingerprint per row is a row of a ``uint8`` matrix (2048 bits = 256 bytes), so 200k
compounds take about 50 MB. Similarities come from bit counts: one query against every row is
a numpy popcount (~10 ms for 200k rows); many against many is a matrix product of 0/1 tiles,
which BLAS does far faster than pairwise loops while holding only a tile in memory.
"""

from __future__ import annotations

import math
import weakref
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text

if TYPE_CHECKING:
    from chemsuite.core.table import Table

Cell = bytes | str | None
CHUNK_ROWS = 2000
TILE_ROWS = 2048  # rows per BLAS tile: 2048 × 2048 bits as float32 is 16 MB
METRICS = {"tanimoto": "Tanimoto", "dice": "Dice", "cosine": "Cosine"}


@dataclass(frozen=True)
class FingerprintKind:
    key: str
    title: str
    bits: int
    help: str = ""


KINDS: dict[str, FingerprintKind] = {
    k.key: k
    for k in (
        FingerprintKind("morgan2", "Morgan r=2 (ECFP4), 2048 bits", 2048, "The usual default."),
        FingerprintKind("morgan3", "Morgan r=3 (ECFP6), 2048 bits", 2048),
        FingerprintKind(
            "fcfp4",
            "Feature Morgan r=2 (FCFP4), 2048 bits",
            2048,
            "Atoms by role (donor, acceptor, aromatic…): finds scaffold hops.",
        ),
        FingerprintKind("rdkit", "RDKit paths (up to 7 bonds), 2048 bits", 2048),
        FingerprintKind("atompair", "Atom pairs, 2048 bits", 2048, "Good for overall shape."),
        FingerprintKind("torsion", "Topological torsions, 2048 bits", 2048),
        FingerprintKind(
            "maccs", "MACCS keys (166)", 167, "Small, interpretable substructure keys."
        ),
    )
}
DEFAULT_KIND = "morgan2"


def _generator(key: str):
    if key in ("morgan2", "morgan3"):
        return rdFingerprintGenerator.GetMorganGenerator(
            radius=2 if key == "morgan2" else 3, fpSize=2048
        )
    if key == "fcfp4":
        return rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=2048,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        )
    if key == "rdkit":
        return rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=7, fpSize=2048)
    if key == "atompair":
        return rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
    if key == "torsion":
        return rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
    return None


def row_bytes(kind: FingerprintKind) -> int:
    """Bytes per packed row, padded to whole 64-bit words for fast popcounts."""
    return math.ceil(kind.bits / 64) * 8


def mol_bits(mol: Chem.Mol, key: str, generator=None) -> np.ndarray:
    """One fingerprint as a bool array of ``KINDS[key].bits``."""
    if key == "maccs":
        vector = MACCSkeys.GenMACCSKeys(mol)
        out = np.zeros(vector.GetNumBits(), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(vector, out)
        return out.astype(np.bool_)
    generator = generator or _generator(key)
    return generator.GetFingerprintAsNumPy(mol).astype(np.bool_)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def pack(bits: np.ndarray, kind: FingerprintKind) -> np.ndarray:
    """Bool rows → packed ``uint8`` rows of :func:`row_bytes` bytes."""
    bits = np.atleast_2d(bits)
    out = np.zeros((len(bits), row_bytes(kind)), dtype=np.uint8)
    packed = np.packbits(bits, axis=1)
    out[:, : packed.shape[1]] = packed
    return out


def fingerprint_chunk(args: tuple[str, Sequence[Cell]]) -> tuple[np.ndarray, np.ndarray]:
    """``(packed rows, valid flags)`` for one chunk (a picklable worker entry point)."""
    key, cells = args
    kind = KINDS[key]
    generator = _generator(key)
    bits = np.zeros((len(cells), kind.bits), dtype=np.bool_)
    valid = np.zeros(len(cells), dtype=np.bool_)
    with rdBase.BlockLogs():
        for i, cell in enumerate(cells):
            mol = _mol(cell)
            if mol is None:
                continue
            try:
                if mol.GetNumConformers() and mol.GetConformer().Is3D():
                    mol = Chem.RemoveHs(mol)  # a 3D cell matches its 2D source
                bits[i] = mol_bits(mol, key, generator)
                valid[i] = True
            except (ValueError, RuntimeError):
                pass
    return pack(bits, kind), valid


@dataclass
class Fingerprints:
    """Packed fingerprints, row for row; invalid rows are all zero and flagged."""

    kind: FingerprintKind
    packed: np.ndarray  # uint8, rows × row_bytes
    valid: np.ndarray  # bool per row
    _counts: np.ndarray | None = field(default=None, repr=False)

    def __len__(self) -> int:
        return len(self.packed)

    @property
    def words(self) -> np.ndarray:
        return self.packed.view(np.uint64)

    @property
    def counts(self) -> np.ndarray:
        """Bits set per row."""
        if self._counts is None:
            self._counts = np.bitwise_count(self.words).sum(axis=1, dtype=np.int32)
        return self._counts

    def take(self, rows: np.ndarray) -> Fingerprints:
        return Fingerprints(self.kind, self.packed[rows], self.valid[rows])

    def unpacked(self, start: int = 0, stop: int | None = None) -> np.ndarray:
        """0/1 float32 rows for BLAS (only the fingerprint's own bits)."""
        bits = np.unpackbits(self.packed[start:stop], axis=1, count=self.kind.bits)
        return bits.astype(np.float32)

    def nbytes(self) -> int:
        return self.packed.nbytes + self.valid.nbytes


def _chunks(key: str, cells: Sequence[Cell], size: int) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        yield key, list(cells[start : start + size])


def compute_fingerprints(
    cells: Sequence[Cell],
    key: str = DEFAULT_KIND,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> Fingerprints:
    kind = KINDS[key]
    parts: list[np.ndarray] = []
    flags: list[np.ndarray] = []
    done = 0
    for packed, valid in ordered_map(fingerprint_chunk, _chunks(key, cells, chunk_rows), executor):
        check(cancel)
        parts.append(packed)
        flags.append(valid)
        done += len(valid)
        if progress is not None:
            progress(done, len(cells))
    if not parts:
        return Fingerprints(kind, np.zeros((0, row_bytes(kind)), np.uint8), np.zeros(0, bool))
    return Fingerprints(kind, np.concatenate(parts), np.concatenate(flags))


def query_fingerprint(mol: Chem.Mol, key: str = DEFAULT_KIND) -> Fingerprints:
    packed, valid = fingerprint_chunk((key, [Chem.MolToSmiles(mol)]))
    if not valid[0]:
        raise ValueError("could not fingerprint the query structure")
    return Fingerprints(KINDS[key], packed, valid)


# -- similarity --------------------------------------------------------------------------------
def from_counts(common: np.ndarray, a: np.ndarray, b: np.ndarray, metric: str) -> np.ndarray:
    """Similarity from shared bits and each side's bit count (0 when both are empty)."""
    common = common.astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        if metric == "tanimoto":
            out = common / (a + b - common)
        elif metric == "dice":
            out = 2 * common / (a + b)
        elif metric == "cosine":
            out = common / np.sqrt(a.astype(np.float64) * b)
        else:
            raise ValueError(f"unknown similarity metric: {metric!r}")
    return np.nan_to_num(out, nan=0.0, posinf=0.0)


def at_least(
    common: np.ndarray, a: np.ndarray, b: np.ndarray, metric: str, threshold: float
) -> np.ndarray:
    """``similarity >= threshold`` without dividing: several times faster on big tiles.

    Tanimoto c/(a+b−c) ≥ t ⇔ c(1+t) ≥ t(a+b); Dice 2c/(a+b) ≥ t ⇔ 2c ≥ t(a+b);
    cosine c/√(ab) ≥ t ⇔ c ≥ t√(ab). Two empty fingerprints (similarity 0) never qualify.
    """
    t = np.float32(threshold)
    if metric == "tanimoto":
        hits = common * (1 + t) >= t * (a + b)
    elif metric == "dice":
        hits = 2 * common >= t * (a + b)
    elif metric == "cosine":
        hits = common >= t * np.sqrt(a * b)
    else:
        raise ValueError(f"unknown similarity metric: {metric!r}")
    return hits & (common > 0) if threshold <= 0 else hits


def similarity(fps: Fingerprints, query: Fingerprints, metric: str = "tanimoto") -> np.ndarray:
    """Similarity of every row to the first row of *query*; NaN for invalid rows."""
    q = query.words[0]
    common = np.bitwise_count(fps.words & q).sum(axis=1, dtype=np.int32)
    out = from_counts(common, fps.counts, query.counts[0], metric)
    return np.where(fps.valid, out, np.nan)


def _tiles(n: int, size: int = TILE_ROWS) -> list[tuple[int, int]]:
    return [(start, min(n, start + size)) for start in range(0, n, size)]


def max_similarity(
    fps: Fingerprints,
    refs: Fingerprints,
    metric: str = "tanimoto",
    *,
    exclude_self: np.ndarray | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(best similarity, index of the best ref)`` per row of *fps* against *refs*.

    *exclude_self* gives, per row of *fps*, the ref index to ignore (the row itself when the
    rows are also refs; -1 for none), so a selected compound is not its own nearest neighbour.
    Invalid rows get NaN and -1.
    """
    n = len(fps)
    best = np.full(n, -1.0)
    which = np.full(n, -1, dtype=np.int64)
    ref_ok = np.flatnonzero(refs.valid)
    if not len(ref_ok) or not n:
        return np.full(n, np.nan), which
    refs = refs.take(ref_ok)
    tiles = _tiles(n)
    for done, (start, stop) in enumerate(tiles):
        check(cancel)
        block = fps.unpacked(start, stop)
        for r_start, r_stop in _tiles(len(refs)):
            common = block @ refs.unpacked(r_start, r_stop).T
            sim = from_counts(
                common, fps.counts[start:stop, None], refs.counts[None, r_start:r_stop], metric
            )
            if exclude_self is not None:
                own = exclude_self[start:stop]
                mapped = np.searchsorted(ref_ok, own)
                hit = (own >= 0) & (mapped < len(ref_ok))
                hit &= ref_ok[np.minimum(mapped, len(ref_ok) - 1)] == own
                local = mapped - r_start
                rows = np.flatnonzero(hit & (local >= 0) & (local < r_stop - r_start))
                sim[rows, local[rows]] = -1.0
            arg = sim.argmax(axis=1)
            top = sim[np.arange(len(sim)), arg]
            better = top > best[start:stop]
            best[start:stop] = np.where(better, top, best[start:stop])
            which[start:stop] = np.where(better, ref_ok[arg + r_start], which[start:stop])
        if progress is not None:
            progress(done + 1, len(tiles))
    best = np.where(best < 0, np.nan, best)
    return np.where(fps.valid, best, np.nan), np.where(fps.valid & (which >= 0), which, -1)


def neighbours(
    fps: Fingerprints,
    threshold: float,
    metric: str = "tanimoto",
    *,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Every pair with similarity ≥ *threshold* as CSR ``(indptr, indices)``, self excluded.

    Invalid rows have no neighbours. Time grows with the square of the rows (about 6 s for
    50k Morgan fingerprints on a laptop); memory with the number of neighbour pairs.
    """
    n = len(fps)
    tiles = _tiles(n)
    counts = fps.counts.astype(np.float32)
    rows_out: list[np.ndarray] = []
    cols_out: list[np.ndarray] = []
    steps = len(tiles) * (len(tiles) + 1) // 2
    step = 0
    unpacked_cache: dict[int, np.ndarray] = {}
    for i, (a0, a1) in enumerate(tiles):
        block_a = fps.unpacked(a0, a1)
        for j in range(i, len(tiles)):
            check(cancel)
            b0, b1 = tiles[j]
            block_b = block_a if i == j else unpacked_cache.get(j)
            if block_b is None:
                block_b = fps.unpacked(b0, b1)
                if len(tiles) <= 8:  # small tables: keep tiles, skip re-unpacking
                    unpacked_cache[j] = block_b
            common = block_a @ block_b.T
            hits = at_least(common, counts[a0:a1, None], counts[None, b0:b1], metric, threshold)
            r, c = np.nonzero(hits)
            r, c = r + a0, c + b0
            keep = r < c  # each pair once (upper triangle), mirrored below
            rows_out.append(r[keep])
            cols_out.append(c[keep])
            step += 1
            if progress is not None:
                progress(step, steps)
    r = np.concatenate(rows_out) if rows_out else np.empty(0, np.int64)
    c = np.concatenate(cols_out) if cols_out else np.empty(0, np.int64)
    ok = fps.valid[r] & fps.valid[c]
    r, c = r[ok], c[ok]
    src = np.concatenate([r, c])
    dst = np.concatenate([c, r])
    order = np.lexsort((dst, src))
    src, dst = src[order], dst[order]
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.add.at(indptr, src + 1, 1)
    return np.cumsum(indptr), dst.astype(np.int64)


# -- a per-table cache -------------------------------------------------------------------------
class FingerprintCache:
    """Fingerprints of a structure column, kept until its cells change.

    Keyed by (column, kind). Rows are matched by oid, so sorting or deleting rows keeps the
    cache; editing the column, adding rows or renaming drops it. :meth:`generation` lets a
    background job tell whether the column changed while it computed.
    """

    MAX_ENTRIES = 3

    def __init__(self, table: Table) -> None:
        from chemsuite.core.table.events import (
            CellsChanged,
            ColumnKindChanged,
            ColumnRemoved,
            ColumnRenamed,
            RowsInserted,
        )

        self._entries: dict[tuple[str, str], tuple[np.ndarray, Fingerprints]] = {}
        self._generations: dict[str, int] = {}
        self._epoch = 0  # bumped when every column is invalidated (rows added)
        kinds = (CellsChanged, ColumnKindChanged, ColumnRemoved)
        ref = weakref.ref(self)

        def listener(event: object) -> None:
            self_ = ref()
            if self_ is None or not hasattr(event, "changes"):
                return
            for change in event.changes:
                if isinstance(change, kinds):
                    self_.drop(change.column if isinstance(change, CellsChanged) else change.name)
                elif isinstance(change, ColumnRenamed):
                    self_.drop(change.old)
                    self_.drop(change.new)
                elif isinstance(change, RowsInserted):
                    self_.drop(None)

        self.unsubscribe = table.subscribe(listener)

    def generation(self, column: str) -> tuple[int, int]:
        return self._epoch, self._generations.get(column, 0)

    def drop(self, column: str | None) -> None:
        """Forget *column* (None: every column)."""
        for key in [k for k in self._entries if column is None or k[0] == column]:
            del self._entries[key]
        if column is None:
            self._epoch += 1
        else:
            self._generations[column] = self._generations.get(column, 0) + 1

    def get(self, column: str, key: str, oids: np.ndarray) -> Fingerprints | None:
        """Fingerprints for *oids* in that order, or None unless every oid is cached."""
        entry = self._entries.get((column, key))
        if entry is None:
            return None
        cached_oids, fps = entry
        pos = np.searchsorted(cached_oids, oids)
        pos = np.minimum(pos, len(cached_oids) - 1) if len(cached_oids) else pos
        if not len(cached_oids) or not np.array_equal(cached_oids[pos], oids):
            return None
        self._entries[(column, key)] = self._entries.pop((column, key))  # most recently used
        return fps.take(pos)

    def put(
        self,
        column: str,
        key: str,
        oids: np.ndarray,
        fps: Fingerprints,
        generation: tuple[int, int],
    ) -> None:
        """Store fingerprints computed while the column was at *generation* (else ignore)."""
        if generation != self.generation(column):
            return
        order = np.argsort(oids, kind="stable")
        old = self._entries.pop((column, key), None)
        sorted_oids, sorted_fps = oids[order], fps.take(order)
        if old is not None:  # merge with what was cached before
            keep = ~np.isin(old[0], sorted_oids)
            merged = np.concatenate([old[0][keep], sorted_oids])
            packed = np.concatenate([old[1].packed[keep], sorted_fps.packed])
            valid = np.concatenate([old[1].valid[keep], sorted_fps.valid])
            order = np.argsort(merged, kind="stable")
            sorted_oids = merged[order]
            sorted_fps = Fingerprints(fps.kind, packed[order], valid[order])
        self._entries[(column, key)] = (sorted_oids, sorted_fps)
        while len(self._entries) > self.MAX_ENTRIES:
            del self._entries[next(iter(self._entries))]


_CACHES: weakref.WeakKeyDictionary = weakref.WeakKeyDictionary()


def cache_for(table: Table) -> FingerprintCache:
    """The table's fingerprint cache (made on first use, gone with the table)."""
    cache = _CACHES.get(table)
    if cache is None:
        cache = _CACHES[table] = FingerprintCache(table)
    return cache
