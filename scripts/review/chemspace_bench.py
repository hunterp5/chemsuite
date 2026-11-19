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

"""Speed and recall of the chemical-space graphs (ROADMAP §7.0d: TMAP and Chemical Space
Networks).

The sets are local data under ``samples/`` (not committed): BindingDB's 49,189 distinct
compounds up to that many, and above it the first *n* distinct compounds of ChEMBL 37's
``chembl_37_chemreps.txt.gz`` (CC BY-SA 3.0) in ``samples/chembl/``. Their SMILES are cached
under ``samples/chemspace/`` on first use. Steps::

    python -I scripts/review/chemspace_bench.py smiles             # cache BindingDB's SMILES
    python -I scripts/review/chemspace_bench.py knn [--n 50000]    # exact against LSH: recall
    python -I scripts/review/chemspace_bench.py tmap [--n 100000]  # the whole TMAP, timed
    python -I scripts/review/chemspace_bench.py csn [--n 5000]     # a network, timed

``knn`` and ``tmap`` take ``--fingerprint ecfp4`` (the default) or ``mhfp6``, and ``--perms``
(MinHash values: the forest's for ECFP4, MHFP6's own). ``tmap`` also writes
``samples/chemspace/tmap_<fingerprint>_<n>.png`` to look at, and reports peak memory and the
tree's weight.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase

from chemsuite.analysis.fingerprints import Fingerprints, compute_fingerprints
from chemsuite.core.jobs.parallel import default_workers

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "chemspace"
SOURCES = ("bindingdb_pubchem_unique.sdf", "bindingdb_subset_50k.sdf")
CHEMBL = ROOT / "samples" / "chembl" / "chembl_37_chemreps.txt.gz"


def peak_mb() -> float:
    """The process's peak memory so far (MB)."""
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]  # fmt: skip

        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        kernel = ctypes.windll.kernel32
        kernel.GetCurrentProcess.restype = wintypes.HANDLE  # a 64-bit -1, not an int
        info = kernel.K32GetProcessMemoryInfo
        info.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        info(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb)
        return counters.PeakWorkingSetSize / 1e6
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 1e6 if sys.platform == "darwin" else peak / 1e3


def smiles() -> list[str]:
    """Every SMILES of the largest available set, cached in ``all.smi``."""
    cache = DATA / "all.smi"
    if cache.exists():
        return cache.read_text(encoding="utf-8").split("\n")
    DATA.mkdir(parents=True, exist_ok=True)
    source = next(ROOT / "samples" / s for s in SOURCES if (ROOT / "samples" / s).exists())
    out: list[str] = []
    seen: set[str] = set()
    with rdBase.BlockLogs():
        for mol in Chem.ForwardSDMolSupplier(str(source)):
            if mol is None:
                continue
            text = Chem.MolToSmiles(mol)
            if text and text not in seen:
                seen.add(text)
                out.append(text)
    cache.write_text("\n".join(out), encoding="utf-8")
    print(f"{len(out):,} unique SMILES from {source.name}")
    return out


def chembl(n: int) -> list[str]:
    """The first *n* distinct SMILES of ChEMBL's chemreps file, cached in ``chembl_<n>.smi``."""
    cache = DATA / f"chembl_{n}.smi"
    if cache.exists():
        return cache.read_text(encoding="utf-8").split("\n")
    out: list[str] = []
    seen: set[str] = set()
    with gzip.open(CHEMBL, "rt", encoding="utf-8") as lines:
        next(lines)  # the header
        for line in lines:
            parts = line.split("\t")
            if len(parts) > 1 and parts[1] and parts[1] not in seen:
                seen.add(parts[1])
                out.append(parts[1])
                if len(out) == n:
                    break
    if len(out) < n:
        raise SystemExit(f"{CHEMBL.name} holds only {len(out):,} compounds")
    DATA.mkdir(parents=True, exist_ok=True)
    cache.write_text("\n".join(out), encoding="utf-8")
    return out


def cells_for(n: int) -> list[str]:
    """*n* SMILES: BindingDB's while it has enough, else ChEMBL's."""
    bindingdb = smiles()
    return bindingdb[:n] if n <= len(bindingdb) else chembl(n)


def _source(n: int) -> str:
    return "" if n <= len(smiles()) else "chembl_"


def fingerprints(n: int) -> Fingerprints:
    """Morgan (r=2) fingerprints of *n* compounds, cached as ``morgan2_<n>.npz``."""
    cache = DATA / f"{_source(n)}morgan2_{n}.npz"
    if cache.exists():
        saved = np.load(cache)
        from chemsuite.analysis.fingerprints import KINDS

        return Fingerprints(KINDS["morgan2"], saved["packed"], saved["valid"])
    from concurrent.futures import ProcessPoolExecutor

    cells = cells_for(n)
    with ProcessPoolExecutor(default_workers()) as pool:
        fps = compute_fingerprints(cells, "morgan2", executor=pool)
    np.savez(cache, packed=fps.packed, valid=fps.valid)
    return fps


def mhfp(n: int, perms: int):
    """MHFP6 of *n* compounds with *perms* values each, cached as ``mhfp6_<perms>_<n>.npz``."""
    from concurrent.futures import ProcessPoolExecutor

    from chemsuite.analysis.graphs.mhfp import MinHashes, compute_mhfp

    cache = DATA / f"{_source(n)}mhfp6_{perms}_{n}.npz"
    if cache.exists():
        saved = np.load(cache)
        return MinHashes(saved["hashes"], saved["valid"])
    cells = cells_for(n)
    with ProcessPoolExecutor(default_workers()) as pool:
        values = _timed("MHFP6", lambda: compute_mhfp(cells, perms, executor=pool))
    np.savez(cache, hashes=values.hashes, valid=values.valid)
    return values


def vectors(n: int, fingerprint: str, perms: int):
    """The valid rows of *n* compounds' fingerprints of either kind."""
    fps = fingerprints(n) if fingerprint == "ecfp4" else mhfp(n, perms)
    return fps.take(np.flatnonzero(fps.valid))


def _timed(label: str, fn):
    start = time.perf_counter()
    out = fn()
    print(f"{label}: {time.perf_counter() - start:.1f} s")
    return out


def run_knn(n: int, k: int, sample: int, trees: int, window: int, perms: int, refine: int,
            fingerprint: str = "ecfp4") -> None:  # fmt: skip
    from chemsuite.analysis.embedding import bit_neighbours
    from chemsuite.analysis.graphs.knn import lsh_neighbours

    if fingerprint == "mhfp6":
        run_knn_mhfp(n, k, sample, trees, window, perms, refine)
        return
    fps = fingerprints(n).take(np.flatnonzero(fingerprints(n).valid))
    print(f"{len(fps):,} valid fingerprints, k = {k}")
    found, sims = _timed(
        f"LSH (perms {perms}, trees {trees}, window {window}, refine {refine})",
        lambda: lsh_neighbours(fps, k, perms=perms, trees=trees, window=window, refine=refine),
    )
    rows = np.random.default_rng(0).choice(len(fps), min(sample, len(fps)), replace=False)
    rows.sort()
    exact_sims = np.empty((len(rows), k))
    for i, row in enumerate(rows):  # exact top-k of the sampled rows by brute force
        words = fps.words
        common = np.bitwise_count(words & words[row]).sum(axis=1)
        sim = common / np.maximum(fps.counts + fps.counts[row] - common, 1)
        sim[row] = -1
        exact_sims[i] = np.sort(sim)[::-1][:k]
    # recall by similarity, so ties between equally similar rows do not count as misses
    got = np.nan_to_num(sims[rows], nan=-1.0)
    recall = np.mean(got >= exact_sims[:, -1:] - 1e-9)
    print(f"recall@{k} on {len(rows):,} rows: {recall:.3f}")
    print(
        f"mean similarity of the k-th neighbour: exact {exact_sims[:, -1].mean():.3f}, "
        f"LSH {np.nanmean(sims[rows, -1]):.3f}"
    )
    if n <= 20_000:
        _timed("exact (bit_neighbours)", lambda: bit_neighbours(fps, k, "jaccard"))


def run_knn_mhfp(n: int, k: int, sample: int, trees: int, window: int, perms: int,
                 refine: int) -> None:  # fmt: skip
    """Recall of the MHFP6 forest against every pair's agreement, for sampled rows."""
    from chemsuite.analysis.graphs.knn import lsh_neighbours

    values = vectors(n, "mhfp6", perms)
    print(f"{len(values):,} valid MHFP6 ({perms} values), k = {k}")
    found, sims = _timed(
        f"LSH (trees {trees}, window {window}, refine {refine})",
        lambda: lsh_neighbours(values, k, trees=trees, window=window, refine=refine),
    )
    rows = np.sort(np.random.default_rng(0).choice(len(values), min(sample, len(values)),
                                                   replace=False))  # fmt: skip
    exact = np.empty((len(rows), k))
    hashes = values.hashes
    for i, row in enumerate(rows):
        sim = (hashes == hashes[row]).mean(axis=1)
        sim[row] = -1
        exact[i] = np.sort(sim)[::-1][:k]
    got = np.nan_to_num(sims[rows], nan=-1.0)
    print(f"recall@{k} on {len(rows):,} rows: {np.mean(got >= exact[:, -1:] - 1e-9):.3f}")
    print(f"mean agreement of the k-th neighbour: exact {exact[:, -1].mean():.3f}, "
          f"LSH {np.nanmean(sims[rows, -1]):.3f}")  # fmt: skip


def draw(path: Path, coords: np.ndarray, a: np.ndarray, b: np.ndarray, size: int = 2000) -> None:
    """The layout as a PNG: grey edges, dark points."""
    from PySide6.QtCore import QLineF, QPointF
    from PySide6.QtGui import QColor, QGuiApplication, QImage, QPainter, QPen

    _app = QGuiApplication.instance() or QGuiApplication([])
    ok = np.isfinite(coords[:, 0])
    lo, hi = coords[ok].min(axis=0), coords[ok].max(axis=0)
    scale = (size - 20) / max(float((hi - lo).max()), 1e-9)
    xy = (coords - lo) * scale + 10
    image = QImage(size, size, QImage.Format.Format_RGB32)
    image.fill(QColor("white"))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setPen(QPen(QColor(150, 150, 150), 0.6))
    painter.drawLines([QLineF(*xy[i], *xy[j]) for i, j in zip(a, b, strict=True)])
    painter.setPen(QPen(QColor(40, 90, 160), 2.0))
    painter.drawPoints([QPointF(*p) for p in xy[ok]])
    painter.end()
    image.save(str(path))
    print(f"wrote {path}")


def run_tmap(n: int, k: int, fingerprint: str = "ecfp4", perms: int = 256) -> None:
    from chemsuite.analysis.graphs.knn import knn, knn_edges
    from chemsuite.analysis.graphs.layout import layout
    from chemsuite.analysis.graphs.trees import components, spanning_forest

    fps = vectors(n, fingerprint, perms)
    print(f"{len(fps):,} valid {fingerprint}, k = {k}; peak so far {peak_mb():,.0f} MB")
    start = time.perf_counter()
    index, sim = _timed("neighbours", lambda: knn(fps, k))
    a, b, s = knn_edges(index, sim)
    ta, tb, tw = _timed("spanning forest", lambda: spanning_forest(len(fps), a, b, 1.0 - s))
    labels = components(len(fps), ta, tb)
    sizes = np.bincount(labels)
    print(f"{len(sizes):,} trees; largest {sizes[0]:,}; singletons {(sizes == 1).sum():,}")
    coords = _timed("layout", lambda: layout(len(fps), ta, tb, labels))
    print(f"total: {time.perf_counter() - start:.1f} s; peak memory {peak_mb():,.0f} MB")
    print(f"tree weight (sum of 1 - similarity): {float(np.sum(tw)):,.1f}")
    draw(DATA / f"tmap_{fingerprint}_{len(fps)}.png", coords, ta, tb)


def run_csn(n: int, threshold: float) -> None:
    from chemsuite.analysis.graphs.csn import network, similarity_edges

    fps = fingerprints(max(n, 10_000)).take(np.arange(n))
    start = time.perf_counter()
    a, b, s = _timed(f"edges at {threshold}", lambda: similarity_edges(fps, threshold))
    net = _timed("layout", lambda: network(n, a, b, s, valid=fps.valid))
    sizes = np.bincount(net.component[net.component >= 0])
    print(
        f"{len(a):,} edges; {len(sizes):,} parts; largest {sizes.max():,}; "
        f"singletons {(sizes == 1).sum():,}; total {time.perf_counter() - start:.1f} s"
    )
    draw(DATA / f"csn_{n}.png", net.coords, net.a, net.b)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="step", required=True)
    sub.add_parser("smiles")
    p = sub.add_parser("knn")
    p.add_argument("--n", type=int, default=50_000)
    p.add_argument("--k", type=int, default=10)
    p.add_argument("--sample", type=int, default=1000)
    p.add_argument("--trees", type=int, default=16)
    p.add_argument("--window", type=int, default=10)
    p.add_argument("--perms", type=int, default=128)
    p.add_argument("--refine", type=int, default=1)
    p.add_argument("--fingerprint", choices=("ecfp4", "mhfp6"), default="ecfp4")
    p = sub.add_parser("tmap")
    p.add_argument("--n", type=int, default=50_000)
    p.add_argument("--k", type=int, default=10)
    p.add_argument("--fingerprint", choices=("ecfp4", "mhfp6"), default="ecfp4")
    p.add_argument("--perms", type=int, default=256)
    p = sub.add_parser("csn")
    p.add_argument("--n", type=int, default=5000)
    p.add_argument("--threshold", type=float, default=0.65)
    args = parser.parse_args()
    if args.step == "smiles":
        print(f"{len(smiles()):,} SMILES")
    elif args.step == "knn":
        run_knn(
            args.n,
            args.k,
            args.sample,
            args.trees,
            args.window,
            args.perms,
            args.refine,
            args.fingerprint,
        )
    elif args.step == "tmap":
        run_tmap(args.n, args.k, args.fingerprint, args.perms)
    elif args.step == "csn":
        run_csn(args.n, args.threshold)


if __name__ == "__main__":
    main()
