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

import functools
import math
import weakref
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from chemsuite.core.jobs import Cancelled, CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text

if TYPE_CHECKING:
    from chemsuite.core.table import Table

Cell = bytes | str | None
CHUNK_ROWS = 2000
TILE_ROWS = 2048  # rows per BLAS tile: 2048 × 2048 bits as float32 is 16 MB
METRICS = {"tanimoto": "Tanimoto", "dice": "Dice", "cosine": "Cosine"}
#: metrics where the order matters (row against reference): offered only by Similarity Search
ASYMMETRIC_METRICS = {"tversky": "Tversky"}
#: Gobbi 2D pharmacophore keys are folded onto this many bits
PHARM2D_BITS = 4096


@dataclass(frozen=True)
class FingerprintKind:
    key: str
    title: str
    bits: int
    help: str = ""
    listed: bool = True  # offered in the similarity tools' fingerprint menus
    short: str = ""  # the column name when written to the table (default: the title)
    pool_min_rows: int = 5000  # fewer rows are fingerprinted quicker in the job's own thread

    @property
    def column(self) -> str:
        return self.short or self.title


KINDS: dict[str, FingerprintKind] = {
    k.key: k
    for k in (
        FingerprintKind(
            "morgan2", "Morgan r=2 (ECFP4), 2048 bits", 2048, "The usual default.", short="ECFP4"
        ),
        FingerprintKind("morgan3", "Morgan r=3 (ECFP6), 2048 bits", 2048, short="ECFP6"),
        FingerprintKind(
            "fcfp4",
            "Feature Morgan r=2 (FCFP4), 2048 bits",
            2048,
            "Atoms by role (donor, acceptor, aromatic…): finds scaffold hops.",
            short="FCFP4",
        ),
        FingerprintKind(
            "rdkit", "RDKit paths (up to 7 bonds), 2048 bits", 2048, short="RDKit Fingerprint"
        ),
        FingerprintKind(
            "atompair", "Atom pairs, 2048 bits", 2048, "Good for overall shape.", short="Atom Pairs"
        ),
        FingerprintKind(
            "torsion", "Topological torsions, 2048 bits", 2048, short="Topological Torsions"
        ),
        FingerprintKind(
            "maccs",
            "MACCS keys (166)",
            167,
            "Small, interpretable substructure keys.",
            short="MACCS Keys",
        ),
        FingerprintKind(
            "pattern",
            "Substructure screen (RDKit pattern), 2048 bits",
            2048,
            "Screens out rows that cannot contain a query; not meant for similarity.",
            listed=False,
        ),
        FingerprintKind(
            "pharm2d",
            f"Gobbi 2D pharmacophore, {PHARM2D_BITS} bits",
            PHARM2D_BITS,
            "Pairs and triplets of pharmacophore features and their topological distances.",
            listed=False,
            short="Pharmacophore 2D",
            pool_min_rows=300,  # about 1 ms a molecule, so the pool pays off sooner
        ),
    )
}
DEFAULT_KIND = "morgan2"

# -- variants: another size, radius, or counts ---------------------------------------------------
#: kinds made by a hashing generator: their size can change, and they can count
SIZED = ("morgan2", "morgan3", "fcfp4", "rdkit", "atompair", "torsion")
#: kinds with a radius
RADIUS = {"morgan2": 2, "morgan3": 3, "fcfp4": 2}
SIZES = (512, 1024, 2048, 4096)
RADII = (1, 2, 3, 4)
DEFAULT_BITS = 2048


@dataclass(frozen=True)
class Variant:
    """A fingerprint kind with its size, radius and counting spelled out."""

    base: str
    bits: int
    radius: int | None
    counts: bool

    @property
    def key(self) -> str:
        """The base key when nothing differs from it, else ``base@b4096,r3,c``."""
        parts = []
        if self.bits != KINDS[self.base].bits:
            parts.append(f"b{self.bits}")
        if self.radius is not None and self.radius != RADIUS[self.base]:
            parts.append(f"r{self.radius}")
        if self.counts:
            parts.append("c")
        return f"{self.base}@{','.join(parts)}" if parts else self.base


def variant(base: str, bits: int | None = None, radius: int | None = None,
            counts: bool = False) -> str:  # fmt: skip
    """The key of *base* (a key of :data:`KINDS`) with another size, radius or counting.

    Only kinds in :data:`SIZED` take a size or counts, and only those in :data:`RADIUS` a
    radius; anything else raises ``ValueError``.
    """
    if base not in KINDS:
        raise ValueError(f"unknown fingerprint {base!r}")
    if (bits not in (None, KINDS[base].bits) or counts) and base not in SIZED:
        raise ValueError(f"{KINDS[base].title} has a fixed size and no counts")
    if radius is not None and base not in RADIUS:
        raise ValueError(f"{KINDS[base].title} has no radius")
    if bits is not None and bits not in SIZES:
        raise ValueError(f"a fingerprint has {', '.join(map(str, SIZES))} bits, not {bits}")
    if radius is not None and radius not in RADII:
        raise ValueError(f"a radius is 1 to {RADII[-1]}, not {radius}")
    return Variant(base, bits or KINDS[base].bits, radius, counts).key


@functools.cache
def parse(key: str) -> Variant:
    """The :class:`Variant` a key names (a plain key is its base kind)."""
    base, _, rest = key.partition("@")
    if base not in KINDS:
        raise ValueError(f"unknown fingerprint {key!r}")
    bits, radius, counts = KINDS[base].bits, RADIUS.get(base), False
    for part in filter(None, rest.split(",")):
        if part == "c":
            counts = True
        elif part[:1] == "b" and part[1:].isdigit():
            bits = int(part[1:])
        elif part[:1] == "r" and part[1:].isdigit():
            radius = int(part[1:])
        else:
            raise ValueError(f"unknown fingerprint {key!r}")
    if rest:  # check the values
        variant(base, bits, radius if base in RADIUS else None, counts)
    return Variant(base, bits, radius, counts)


@functools.cache
def kind_of(key: str) -> FingerprintKind:
    """The :class:`FingerprintKind` of any key, a variant's with its own title and bits."""
    if key in KINDS:
        return KINDS[key]
    v = parse(key)
    base = KINDS[v.base]
    if v.base in RADIUS:
        family = "FCFP" if v.base == "fcfp4" else "ECFP"
        name = f"{'Feature Morgan' if v.base == 'fcfp4' else 'Morgan'} r={v.radius}"
        short = f"{family}{2 * v.radius}"
        title = f"{name} ({short}), {v.bits} bits"
    else:
        short = base.column
        title = base.title.replace(f"{base.bits} bits", f"{v.bits} bits")
    if v.bits != base.bits:
        short += f" {v.bits}"
    if v.counts:
        title, short = f"{title}, counts", f"{short} counts"
    return FingerprintKind(key, title, v.bits, listed=False, short=short,
                           pool_min_rows=base.pool_min_rows)  # fmt: skip


def _generator(key: str):
    v = parse(key)
    if v.base in ("morgan2", "morgan3"):
        return rdFingerprintGenerator.GetMorganGenerator(radius=v.radius, fpSize=v.bits)
    if v.base == "fcfp4":
        return rdFingerprintGenerator.GetMorganGenerator(
            radius=v.radius,
            fpSize=v.bits,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        )
    if v.base == "rdkit":
        return rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=7, fpSize=v.bits)
    if v.base == "atompair":
        return rdFingerprintGenerator.GetAtomPairGenerator(fpSize=v.bits)
    if v.base == "torsion":
        return rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=v.bits)
    return None


def row_bytes(kind: FingerprintKind) -> int:
    """Bytes per packed row, padded to whole 64-bit words for fast popcounts."""
    return math.ceil(kind.bits / 64) * 8


def mol_bits(mol: Chem.Mol, key: str, generator=None) -> np.ndarray:
    """One fingerprint as a bool array of ``kind_of(key).bits``."""
    if key == "maccs":
        vector = MACCSkeys.GenMACCSKeys(mol)
        out = np.zeros(vector.GetNumBits(), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(vector, out)
        return out.astype(np.bool_)
    if key == "pattern":
        vector = Chem.PatternFingerprint(mol, fpSize=KINDS[key].bits)
        out = np.zeros(vector.GetNumBits(), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(vector, out)
        return out.astype(np.bool_)
    if key == "pharm2d":
        from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

        out = np.zeros(PHARM2D_BITS, dtype=np.bool_)
        keys = Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory).GetOnBits()
        out[np.fromiter(keys, np.int64) % PHARM2D_BITS] = True
        return out
    generator = generator or _generator(key)
    return generator.GetFingerprintAsNumPy(mol).astype(np.bool_)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def cell_mol(cell: Cell) -> Chem.Mol | None:
    """The molecule fingerprints and searches see: a 3D cell without its explicit hydrogens,
    so it matches its 2D source. None for a blank or invalid cell."""
    mol = _mol(cell)
    if mol is not None and mol.GetNumConformers() and mol.GetConformer().Is3D():
        mol = Chem.RemoveHs(mol)
    return mol


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
    kind = kind_of(key)
    generator = _generator(key)
    bits = np.zeros((len(cells), kind.bits), dtype=np.bool_)
    valid = np.zeros(len(cells), dtype=np.bool_)
    with rdBase.BlockLogs():
        for i, cell in enumerate(cells):
            try:
                mol = cell_mol(cell)
                if mol is None:
                    continue
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
    kind = kind_of(key)
    parts: list[np.ndarray] = []
    flags: list[np.ndarray] = []
    done = 0
    for packed, valid in ordered_map(
        fingerprint_chunk, _chunks(key, cells, chunk_rows), executor, cancel=cancel
    ):
        check(cancel)
        parts.append(packed)
        flags.append(valid)
        done += len(valid)
        if progress is not None:
            progress(done, len(cells))
    if not parts:
        return Fingerprints(kind, np.zeros((0, row_bytes(kind)), np.uint8), np.zeros(0, bool))
    return Fingerprints(kind, np.concatenate(parts), np.concatenate(flags))


def count_chunk(args: tuple[str, Sequence[Cell]]) -> tuple[np.ndarray, np.ndarray]:
    """``(counts, valid flags)`` for one chunk: how often each hashed feature occurs, capped
    at 255 (``uint8``), for a key whose base is in :data:`SIZED` (a picklable entry point)."""
    key, cells = args
    v = parse(key)
    if v.base not in SIZED:
        raise ValueError(f"{KINDS[v.base].title} has no counts")
    generator = _generator(key)
    out = np.zeros((len(cells), v.bits), dtype=np.uint8)
    valid = np.zeros(len(cells), dtype=np.bool_)
    with rdBase.BlockLogs():
        for i, cell in enumerate(cells):
            try:
                mol = cell_mol(cell)
                if mol is None:
                    continue
                counts = generator.GetCountFingerprintAsNumPy(mol)
                out[i] = np.minimum(counts, 255)
                valid[i] = True
            except (ValueError, RuntimeError):
                pass
    return out, valid


def compute_counts(
    cells: Sequence[Cell],
    key: str,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> tuple[np.ndarray, np.ndarray]:
    """``(counts, valid flags)`` for every cell: see :func:`count_chunk`."""
    parts: list[np.ndarray] = []
    flags: list[np.ndarray] = []
    done = 0
    for counts, valid in ordered_map(
        count_chunk, _chunks(key, cells, chunk_rows), executor, cancel=cancel
    ):
        check(cancel)
        parts.append(counts)
        flags.append(valid)
        done += len(valid)
        if progress is not None:
            progress(done, len(cells))
    if not parts:
        return np.zeros((0, kind_of(key).bits), np.uint8), np.zeros(0, bool)
    return np.concatenate(parts), np.concatenate(flags)


def as_hex(fps: Fingerprints) -> list[str]:
    """Each row's bits as hexadecimal text, bit 0 first (the high bit of the first byte); a
    row that could not be fingerprinted is ""."""
    width = math.ceil(fps.kind.bits / 8)
    return [
        row[:width].tobytes().hex() if ok else ""
        for row, ok in zip(fps.packed, fps.valid.tolist(), strict=True)
    ]


def query_fingerprint(mol: Chem.Mol, key: str = DEFAULT_KIND) -> Fingerprints:
    packed, valid = fingerprint_chunk((key, [Chem.MolToSmiles(mol)]))
    if not valid[0]:
        raise ValueError("could not fingerprint the query structure")
    return Fingerprints(kind_of(key), packed, valid)


# -- similarity --------------------------------------------------------------------------------
def from_counts(
    common: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    metric: str,
    alpha: float = 0.5,
    beta: float = 0.5,
) -> np.ndarray:
    """Similarity from shared bits and each side's bit count (0 when both are empty).

    Tversky weighs the bits only *b* has (the reference) by *alpha* and those only *a* has (the
    row) by *beta*: c / (alpha (b - c) + beta (a - c) + c). alpha = 1, beta = 0 is the share of
    the reference's bits the row has; 0.5 and 0.5 is Dice; 1 and 1 is Tanimoto.
    """
    common = common.astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        if metric == "tversky":
            out = common / (alpha * (b - common) + beta * (a - common) + common)
        elif metric == "tanimoto":
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


def similarity(
    fps: Fingerprints,
    query: Fingerprints,
    metric: str = "tanimoto",
    alpha: float = 0.5,
    beta: float = 0.5,
) -> np.ndarray:
    """Similarity of every row to the first row of *query*; NaN for invalid rows."""
    q = query.words[0]
    common = np.bitwise_count(fps.words & q).sum(axis=1, dtype=np.int32)
    out = from_counts(common, fps.counts, query.counts[0], metric, alpha, beta)
    return np.where(fps.valid, out, np.nan)


def _tiles(n: int, size: int | None = None) -> list[tuple[int, int]]:
    size = size or TILE_ROWS  # read at call time, so tests can make tiles small
    return [(start, min(n, start + size)) for start in range(0, n, size)]


def max_similarity(
    fps: Fingerprints,
    refs: Fingerprints,
    metric: str = "tanimoto",
    *,
    exclude_self: np.ndarray | None = None,
    alpha: float = 0.5,
    beta: float = 0.5,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(best similarity, index of the best ref)`` per row of *fps* against *refs*.

    *exclude_self* gives, per row of *fps*, the ref index to ignore (the row itself when the
    rows are also refs; -1 for none), so a selected compound is not its own nearest neighbour.
    Invalid rows get NaN and -1.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` is ``(best, which)`` for
    the first rows (those finished), or None if none were.
    """
    return _scan(fps, refs, metric, exclude_self, alpha, beta, False, progress, cancel)


def mean_similarity(
    fps: Fingerprints,
    refs: Fingerprints,
    metric: str = "tanimoto",
    *,
    exclude_self: np.ndarray | None = None,
    alpha: float = 0.5,
    beta: float = 0.5,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(mean similarity over the refs, index of the best ref)`` per row, as
    :func:`max_similarity` (a row's own ref, when excluded, is left out of its mean)."""
    return _scan(fps, refs, metric, exclude_self, alpha, beta, True, progress, cancel)


def _scan(
    fps: Fingerprints,
    refs: Fingerprints,
    metric: str,
    exclude_self: np.ndarray | None,
    alpha: float,
    beta: float,
    mean: bool,
    progress: Callable[[int, int], None] | None,
    cancel: CancelToken | None,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(fps)
    best = np.full(n, -1.0)
    which = np.full(n, -1, dtype=np.int64)
    total = np.zeros(n)
    count = np.zeros(n)
    ref_ok = np.flatnonzero(refs.valid)
    if not len(ref_ok) or not n:
        return np.full(n, np.nan), which

    def result(stop: int) -> tuple[np.ndarray, np.ndarray]:
        if mean:
            with np.errstate(invalid="ignore", divide="ignore"):
                score = np.where(count[:stop] > 0, total[:stop] / count[:stop], -1.0)
        else:
            score = best[:stop]
        return _best(fps.take(np.arange(stop)), score, which[:stop])

    refs = refs.take(ref_ok)
    tiles = _tiles(n)
    finished = 0  # rows done
    for done, (start, stop) in enumerate(tiles):
        if cancel is not None and cancel.cancelled:
            if not finished:
                raise Cancelled
            raise Cancelled(result(finished))
        block = fps.unpacked(start, stop)
        for r_start, r_stop in _tiles(len(refs)):
            common = block @ refs.unpacked(r_start, r_stop).T
            sim = from_counts(common, fps.counts[start:stop, None],
                              refs.counts[None, r_start:r_stop], metric, alpha, beta)  # fmt: skip
            own = np.zeros(sim.shape, dtype=np.bool_)
            if exclude_self is not None:
                mine = exclude_self[start:stop]
                mapped = np.searchsorted(ref_ok, mine)
                hit = (mine >= 0) & (mapped < len(ref_ok))
                hit &= ref_ok[np.minimum(mapped, len(ref_ok) - 1)] == mine
                local = mapped - r_start
                rows = np.flatnonzero(hit & (local >= 0) & (local < r_stop - r_start))
                own[rows, local[rows]] = True
                sim[own] = -1.0
            total[start:stop] += np.where(own, 0.0, sim).sum(axis=1)
            count[start:stop] += (~own).sum(axis=1)
            arg = sim.argmax(axis=1)
            top = sim[np.arange(len(sim)), arg]
            better = top > best[start:stop]
            best[start:stop] = np.where(better, top, best[start:stop])
            which[start:stop] = np.where(better, ref_ok[arg + r_start], which[start:stop])
        finished = stop
        if progress is not None:
            progress(done + 1, len(tiles))
    return result(n)


def _best(fps: Fingerprints, best: np.ndarray, which: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    best = np.where(best < 0, np.nan, best)
    return np.where(fps.valid, best, np.nan), np.where(fps.valid & (which >= 0), which, -1)


class TooManyPairs(ValueError):
    """More neighbour pairs than asked for (at least :attr:`count`)."""

    def __init__(self, count: int) -> None:
        super().__init__(f"at least {count:,} pairs")
        self.count = count


def neighbours(
    fps: Fingerprints,
    threshold: float,
    metric: str = "tanimoto",
    *,
    max_pairs: int | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Every pair with similarity ≥ *threshold* as CSR ``(indptr, indices)``, self excluded.

    Invalid rows have no neighbours. Time grows with the square of the rows (about 6 s for
    50k Morgan fingerprints on a laptop); memory with the number of neighbour pairs, so past
    *max_pairs* pairs it stops with :class:`TooManyPairs`.
    """
    n = len(fps)
    tiles = _tiles(n)
    counts = fps.counts.astype(np.float32)
    rows_out: list[np.ndarray] = []
    cols_out: list[np.ndarray] = []
    steps = len(tiles) * (len(tiles) + 1) // 2
    step = found = 0
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
            found += int(keep.sum())
            if max_pairs is not None and found > max_pairs:
                raise TooManyPairs(found)
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


def fingerprint_task(
    table: Table,
    column: str,
    key: str,
    oids: np.ndarray,
    pool: Callable[[], Executor],
) -> tuple[Callable[[object], Fingerprints], Callable[[Fingerprints], None]]:
    """``(compute, remember)`` for fingerprints of *oids* (in that order) through the cache.

    Call this on the GUI thread: it reads the cells now. *compute(job)* runs in the job (a
    cache hit is instant; big tables use *pool()*); *remember(fps)* stores the result once the
    job is done, unless the column changed meanwhile.
    """
    cache = cache_for(table)
    cached = cache.get(column, key, oids)
    generation = cache.generation(column)
    cells = _cells(table, column, oids) if cached is None else None
    use_pool = cells is not None and len(cells) >= kind_of(key).pool_min_rows

    def compute(job) -> Fingerprints:
        if cached is not None:
            return cached
        return compute_fingerprints(
            cells,
            key,
            executor=pool() if use_pool else None,
            progress=getattr(job, "report", None),
            cancel=getattr(job, "cancel_token", None),
        )

    def remember(fps: Fingerprints) -> None:
        if cached is None:
            cache.put(column, key, oids, fps, generation)

    return compute, remember


def _cells(table: Table, column: str, oids: np.ndarray) -> list[Cell]:
    from chemsuite.core.table.columns import flat_inputs

    return flat_inputs(table.column(column), table.slots_of(oids))


def counts_task(
    table: Table,
    column: str,
    key: str,
    oids: np.ndarray,
    pool: Callable[[], Executor],
) -> Callable[..., tuple[np.ndarray, np.ndarray]]:
    """*compute(progress, cancel)* -> ``(counts, valid)`` of *oids* (see :func:`count_chunk`);
    counts are not cached. Call this on the GUI thread: it reads the cells now."""
    cells = _cells(table, column, oids)
    use_pool = len(cells) >= kind_of(key).pool_min_rows

    def compute(
        progress: Callable[[int, int], None] | None = None, cancel: CancelToken | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        executor = pool() if use_pool else None
        return compute_counts(cells, key, executor=executor, progress=progress, cancel=cancel)

    return compute
