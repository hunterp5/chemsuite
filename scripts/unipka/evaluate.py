#!/usr/bin/env python3
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

"""chemsuite pKa weights on the held-out test sets, against MolGpKa (docs/UNIPKA_PLAN.md, U2).

    python scripts/unipka/evaluate.py <weights.npz> <tests folder> <cache.pkl>
        [--workers N] [--enumerated] [--sets novartis_acid,sampl6,...]

Each test row is a transition, "acid microstates>>base microstates", and its measured pKa.
Every microstate's free energy is computed once per fold and conformer (upstream's ten 3D
conformers and the 2D layout) and kept in *cache.pkl*, so the table shows what fewer
conformers and folds cost without running the network again.

- **given:** the pKa between the row's own microstates, as Uni-pKa's paper scores it;
- **enumerated** (``--enumerated``): from the molecule alone, as the app will: every microstate
  from chemsuite's enumerator, the macro pKa between the transition's two charges;
- **MolGpKa:** its site closest to the measured value, among sites of the transition's kind
  (an optimistic bound: the app cannot know which site is meant).
"""

from __future__ import annotations

import os

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")  # oversubscribed BLAS is pathologically slow here

import argparse  # noqa: E402
import csv  # noqa: E402
import pickle  # noqa: E402
from collections.abc import Callable, Iterable  # noqa: E402
from dataclasses import dataclass  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402

from chemsuite.ionisation import molgpka  # noqa: E402
from chemsuite.ionisation.unipka import conformers, microstates, model, thermo  # noqa: E402

SETS = ("novartis_acid", "novartis_base", "sampl6", "sampl7", "sampl8")
N_CONFORMERS = conformers.N_3D + 1  # the 2D layout last
CONFORMER_STEPS = (0, 1, 3, conformers.N_3D)  # 3D conformers kept (always with the 2D layout)


@dataclass(frozen=True)
class Row:
    test: str
    acid: tuple[str, ...]
    base: tuple[str, ...]
    pka: float

    @property
    def charge(self) -> int:
        """The acid side's total charge (the transition is charge → charge − 1)."""
        return Chem.GetFormalCharge(Chem.MolFromSmiles(self.acid[0]))


def read_rows(folder: Path, sets: Iterable[str]) -> list[Row]:
    out = []
    for name in sets:
        with (folder / f"{name}.tsv").open(encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                acid, base = r["SMILES"].split(">>")
                pair = tuple(acid.split(",")), tuple(base.split(","))
                out.append(Row(name, *pair, float(r["TARGET"])))
    return out


# -- free energies, in worker processes --------------------------------------------------------
_NETS: list[model.Network] = []


def _folds(path: Path) -> list[int]:
    with np.load(path) as archive:
        return sorted({int(k.split("/", 1)[0][1:]) for k in archive.files if "/" in k})


def _init(path: Path) -> None:
    RDLogger.DisableLog("rdApp.*")
    _NETS.extend(model.Network.load(path, k) for k in _folds(path))


def _energies(smiles: str) -> tuple[str, np.ndarray]:
    """(folds, conformers) free energies in kT; NaN when RDKit cannot read *smiles*."""
    if Chem.MolFromSmiles(smiles) is None:
        return smiles, np.full((len(_NETS), N_CONFORMERS), np.nan, dtype=np.float32)
    found = conformers.atoms(smiles)
    g = [model.free_energy(net, list(found.symbols), list(found.charges), found.coords)
         for net in _NETS]  # fmt: skip
    return smiles, np.array(g, dtype=np.float32)


def fill(cache: dict[str, np.ndarray], smiles: set[str], weights: Path, workers: int) -> None:
    todo = sorted(smiles - cache.keys())
    print(f"microstates: {len(smiles)}, to compute: {len(todo)}", flush=True)
    if not todo:
        return
    with Pool(workers, _init, (weights,)) as pool:
        for k, (s, g) in enumerate(pool.imap_unordered(_energies, todo, chunksize=4), 1):
            cache[s] = g
            if k % 50 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)}", flush=True)


# -- predictions --------------------------------------------------------------------------------
class Means:
    """The archive's per-fold means and their average (``thermo.common_mean``)."""

    folds: np.ndarray = np.zeros(0)
    mean: float = 0.0


def _g(smiles: str, g: np.ndarray, n_3d: int, folds: int) -> float:
    per_fold = g[:folds, [*range(n_3d), N_CONFORMERS - 1]].mean(axis=1)
    charge = Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))
    return thermo.common_mean(per_fold, charge, Means.folds[:folds], Means.mean)


def given(row: Row, cache: dict[str, np.ndarray], mean: float, n_3d: int, folds: int) -> float:
    ga = np.array([_g(s, cache[s], n_3d, folds) for s in row.acid])
    gb = np.array([_g(s, cache[s], n_3d, folds) for s in row.base])
    return thermo.macro_pka(ga, gb, mean)


def _molecule(row: Row) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(row.acid[0])
    return None if mol is None else molgpka.prepare(mol)


def enumerate_states(row: Row) -> dict[int, list[str]]:
    mol = _molecule(row)
    if mol is None:
        return {}
    return {q: [m.smiles for m in pool] for q, pool in microstates.ensemble(mol).items()}


def enumerated(row: Row, states: dict[int, list[str]], cache: dict[str, np.ndarray],
               mean: float, n_3d: int, folds: int) -> float:  # fmt: skip
    q = row.charge
    if not states.get(q) or not states.get(q - 1):
        return np.nan
    ga = np.array([_g(s, cache[s], n_3d, folds) for s in states[q]])
    gb = np.array([_g(s, cache[s], n_3d, folds) for s in states[q - 1]])
    return thermo.macro_pka(ga, gb, mean)


def molgpka_closest(row: Row) -> float:
    mol = _molecule(row)
    sites = [] if mol is None else molgpka.predict_sites(mol)
    kind = molgpka.BASE if row.charge > 0 else molgpka.ACID
    pool = [s.pka for s in sites if s.kind == kind] or [s.pka for s in sites]
    return min(pool, key=lambda p: abs(p - row.pka)) if pool else np.nan


# -- the table ----------------------------------------------------------------------------------
def _stats(want: np.ndarray, got: np.ndarray) -> tuple[float, float, float, int]:
    ok = ~np.isnan(got)
    err = got[ok] - want[ok]
    r2 = 1 - (err**2).sum() / ((want[ok] - want[ok].mean()) ** 2).sum()
    return float(np.sqrt((err**2).mean())), float(np.abs(err).mean()), float(r2), int((~ok).sum())


def table(rows: list[Row], methods: dict[str, Callable[[Row], float]]) -> str:
    sets = sorted({r.test for r in rows}, key=SETS.index)
    want = np.array([r.pka for r in rows])
    head = ["method", *[f"{s} ({sum(r.test == s for r in rows)})" for s in sets],
            "all RMSE", "all MAE", "all R²", "failed"]  # fmt: skip
    lines = ["| " + " | ".join(head) + " |", "|---" * len(head) + "|"]
    for name, predict in methods.items():
        got = np.array([predict(r) for r in rows], dtype=np.float64)
        cells = [name]
        for s in sets:
            mask = np.array([r.test == s for r in rows])
            cells.append(f"{_stats(want[mask], got[mask])[0]:.2f}")
        rmse, mae, r2, failed = _stats(want, got)
        cells += [f"{rmse:.2f}", f"{mae:.2f}", f"{r2:.2f}", str(failed)]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("weights", type=Path)
    parser.add_argument("tests", type=Path)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--enumerated", action="store_true")
    parser.add_argument("--sets", default=",".join(SETS))
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    rows = read_rows(args.tests, args.sets.split(","))
    with np.load(args.weights) as archive:
        mean = float(archive["mean"])
        Means.folds = np.asarray(archive["fold_means"], dtype=np.float64)
        Means.mean = mean
    n_folds = len(_folds(args.weights))
    store = {"weights": args.weights.name, "energies": {}, "states": {}}
    if args.cache.exists():
        store = pickle.loads(args.cache.read_bytes())
        if store["weights"] != args.weights.name:
            raise SystemExit(f"{args.cache} holds energies for {store['weights']}")
    cache, states = store["energies"], store["states"]
    wanted = {s for r in rows for s in (*r.acid, *r.base)}
    if args.enumerated:
        for r in rows:
            key = r.acid[0]
            if key not in states:
                states[key] = enumerate_states(r)
            wanted |= {s for pool in states[key].values() for s in pool}
    try:
        fill(cache, wanted, args.weights, args.workers)
    finally:
        args.cache.write_bytes(pickle.dumps(store))

    methods: dict[str, Callable[[Row], float]] = {"MolGpKa (closest site)": molgpka_closest}
    for folds in sorted({1, n_folds}):
        for n_3d in CONFORMER_STEPS:
            label = f"given, {n_3d} 3D + 2D, {folds} fold{'s' * (folds > 1)}"
            methods[label] = lambda r, n=n_3d, f=folds: given(r, cache, mean, n, f)
    if args.enumerated:
        for n_3d in (1, conformers.N_3D):
            methods[f"enumerated, {n_3d} 3D + 2D, {n_folds} folds"] = lambda r, n=n_3d: enumerated(
                r, states[r.acid[0]], cache, mean, n, n_folds
            )
    print(f"\n{args.weights.name}: mean pKa {mean:.4f}, {n_folds} folds\n")
    print(table(rows, methods))


if __name__ == "__main__":
    main()
