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
"""chemsuite's TMAP beside the reference implementation, ``tmap-silicon`` (ROADMAP §7.0d).

Both read the same MHFP6 values (``samples/chemspace/mhfp6_256_<n>.npz``, made by
``chemspace_bench.py knn --fingerprint mhfp6 --n <n>``), so only the neighbour search, the
tree and the layout differ. ``tmap-silicon`` is never a chemsuite dependency: it runs in a
virtual environment of its own (``pip install tmap-silicon numpy``)::

    <tmap venv>/python -I scripts/review/tmap_compare.py tmap --n 49189      # its TMAP
    .venv/python -I scripts/review/tmap_compare.py report --n 49189          # both, compared

``report`` scores both against every pair's MinHash agreement for sampled rows (recall of the
ten nearest), compares the trees' total weight (the sum of 1 − agreement over their edges: the
lower, the closer the neighbours each tree joins) and writes both pictures under
``samples/chemspace/``.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "chemspace"
PERMS = 256
K = 10


def values(n: int) -> tuple[np.ndarray, np.ndarray]:
    saved = np.load(DATA / f"mhfp6_{PERMS}_{n}.npz")
    return saved["hashes"], saved["valid"]


def run_tmap(n: int, trees: int) -> None:
    """tmap-silicon's TMAP of the valid rows: its k-NN graph, tree and layout, saved."""
    import tmap as tm

    hashes, valid = values(n)
    hashes = hashes[valid]
    start = time.perf_counter()
    forest = tm.LSHForest(PERMS, trees)
    forest.batch_add([tm.VectorUint(row.tolist()) for row in hashes])
    forest.index()
    a, b, w = tm.VectorUint(), tm.VectorUint(), tm.VectorFloat()
    forest.get_knn_graph(a, b, w, K + 1, 10)  # it counts each compound among its own
    knn_time = time.perf_counter() - start
    config = tm.LayoutConfiguration()
    config.k, config.kc = K + 1, 10  # so k + 1 gives k others
    x, y, s, t, _ = tm.layout_from_lsh_forest(forest, config)
    total = time.perf_counter() - start
    np.savez(
        DATA / f"tmap_silicon_{n}.npz",
        knn_a=np.array(a), knn_b=np.array(b), x=np.array(x), y=np.array(y),
        s=np.array(s), t=np.array(t), knn_time=knn_time, total_time=total, trees=trees,
    )  # fmt: skip
    print(f"tmap-silicon (l = {trees}): k-NN {knn_time:.1f} s, total {total:.1f} s")


def _recall(index_of: list[np.ndarray], hashes: np.ndarray, rows: np.ndarray) -> float:
    """The share of each sampled row's ten nearest (by agreement) found among *index_of*."""
    hits = []
    for row in rows:
        sim = (hashes == hashes[row]).mean(axis=1)
        sim[row] = -1
        kth = np.sort(sim)[::-1][K - 1]
        found = index_of[row]
        found = found[found != row]
        hits.append(min(int((sim[found] >= kth - 1e-9).sum()), K) / K)
    return float(np.mean(hits))


def _weight(hashes: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(1.0 - (hashes[a] == hashes[b]).mean(axis=1)))


def report(n: int, sample: int) -> None:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from chemspace_bench import draw

    from chemsuite.analysis.graphs.knn import knn
    from chemsuite.analysis.graphs.mhfp import MinHashes
    from chemsuite.analysis.graphs.tmap import tmap
    from chemsuite.analysis.graphs.trees import components

    hashes, valid = values(n)
    hashes = hashes[valid]
    m = len(hashes)
    theirs = np.load(DATA / f"tmap_silicon_{n}.npz")
    rows = np.random.default_rng(0).choice(m, min(sample, m), replace=False)

    mine = MinHashes(hashes, np.ones(m, bool))
    start = time.perf_counter()
    index, _sim = knn(mine, K)
    knn_time = time.perf_counter() - start
    start = time.perf_counter()
    tree = tmap(mine, K)  # the neighbours again, then the tree and the layout
    total = time.perf_counter() - start

    their_lists: list[list[int]] = [[] for _ in range(m)]
    for i, j in zip(theirs["knn_a"].tolist(), theirs["knn_b"].tolist(), strict=True):
        their_lists[i].append(j)
    their_index = [np.array(x, dtype=np.int64) for x in their_lists]
    my_index = [r[r >= 0] for r in index]

    s, t = theirs["s"].astype(np.int64), theirs["t"].astype(np.int64)
    their_trees = int(components(m, s, t).max()) + 1
    my_trees = int(tree.component.max()) + 1
    print(f"{m:,} compounds, MHFP6 with {PERMS} values, k = {K}; recall on {len(rows)} rows")
    print("| | chemsuite | tmap-silicon |")
    print("|---|---|---|")
    print(f"| recall of the ten nearest | {_recall(my_index, hashes, rows):.3f} | "
          f"{_recall(their_index, hashes, rows):.3f} |")  # fmt: skip
    print(f"| tree edges | {len(tree.a):,} | {len(s):,} |")
    print(f"| separate trees | {my_trees:,} | {their_trees:,} |")
    print(f"| tree weight (sum of 1 - agreement) | {_weight(hashes, tree.a, tree.b):,.1f} | "
          f"{_weight(hashes, s, t):,.1f} |")  # fmt: skip
    print(f"| neighbours (s) | {knn_time:.1f} | {float(theirs['knn_time']):.1f} |")
    print(f"| whole TMAP (s) | {total:.1f} | {float(theirs['total_time']):.1f} |")
    draw(DATA / f"compare_chemsuite_{m}.png", tree.coords, tree.a, tree.b)
    xy = np.c_[theirs["x"], theirs["y"]].astype(float)
    draw(DATA / f"compare_tmap_silicon_{m}.png", xy, s, t)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("tmap")
    p.add_argument("--n", type=int, default=49_189)
    p.add_argument("--trees", type=int, default=16)
    p = sub.add_parser("report")
    p.add_argument("--n", type=int, default=49_189)
    p.add_argument("--sample", type=int, default=300)
    args = parser.parse_args()
    if args.step == "tmap":
        run_tmap(args.n, args.trees)
    else:
        report(args.n, args.sample)


if __name__ == "__main__":
    main()
