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

"""What a CYP inhibition model can honestly reach (Phase C, docs/SOM_PLAN.md).

    .venv/Scripts/python.exe scripts/cyp/benchmark.py            # the baseline, per isoform
    .venv/Scripts/python.exe scripts/cyp/benchmark.py --engine   # through Auto QSAR as well

Reads ``samples/cyp/instances.csv`` (see curate.py) and measures, per isoform, what a model
fitted on ECFP4 counts and the cheap 2D descriptors reaches on compounds it has not seen.

**Scaffold-grouped folds, never random.** A random split on a set with congeneric series puts
near-copies of a test compound in training and flatters every model — the same mistake the
sites-of-metabolism work found in ATTNSOM's published numbers, where it was worth 0.18 MCC. The
folds here keep a Bemis-Murcko scaffold whole.

**A shuffled-y null runs beside the real fit.** If a pipeline scores above noise on shuffled
labels, the number is wrong and nothing else in the run means anything. It is the cheapest
guard there is against a leak, so it is not optional.

**The 750-compound test file is a sanity check, not a blind one.** It was unblinded on
2026-11-04, so anyone may have trained on it. It is reported because an external set still
catches gross mistakes, and it is labelled as what it is wherever it is quoted.

The baseline is a plain ``HistGradientBoostingRegressor``, which is what the engine's own zoo
would reach for and which needs nothing that is still being written. ``--engine`` additionally
runs chemsuite's Auto QSAR engine (``chemsuite/qsar``) over the same pool and folds, which both
measures what the full search adds and exercises that engine on a real, messy dataset.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase

REPO = Path(__file__).resolve().parents[2]
INSTANCES = REPO / "samples" / "cyp" / "instances.csv"

ISOFORMS = ("1A2", "2C9", "2D6", "3A4")
#: Folds per run and how many seeds the whole thing is repeated over, as the SoM work reports.
FOLDS, SEEDS = 5, 3
#: The fewest rows an isoform needs before it is worth reporting at all.
MIN_ROWS = 100


@dataclass(frozen=True)
class Data:
    """The curated table, as arrays."""

    smiles: tuple[str, ...]
    mols: tuple[Chem.Mol | None, ...]
    split: np.ndarray  # "train" or "test"
    values: dict[str, np.ndarray]  # isoform -> pIC50, NaN where not measured

    def rows_for(self, isoform: str, split: str) -> np.ndarray:
        """Rows of *split* that have a value for *isoform*."""
        has = np.isfinite(self.values[isoform])
        return np.flatnonzero(has & (self.split == split))


def read(path: Path) -> Data:
    """The curated table."""
    smiles, split = [], []
    values: dict[str, list[float]] = {isoform: [] for isoform in ISOFORMS}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            smiles.append(row["smiles"])
            split.append(row["split"])
            for isoform in ISOFORMS:
                text = (row[f"pIC50_{isoform}"] or "").strip()
                values[isoform].append(float(text) if text else float("nan"))
    with rdBase.BlockLogs():
        mols = tuple(Chem.MolFromSmiles(s) for s in smiles)
    return Data(
        smiles=tuple(smiles),
        mols=mols,
        split=np.array(split),
        values={k: np.asarray(v, dtype=np.float64) for k, v in values.items()},
    )


def pool_for(data: Data) -> tuple[np.ndarray, tuple[str, ...]]:
    """ECFP4 counts and the cheap 2D descriptors, the features Phase C named.

    Both come from the app's own code, so a feature here means the same number the table's
    Calculate Descriptors and similarity tools would give. Predicted pKa is deliberately left
    out: it costs seconds per molecule against microseconds for the rest, and
    ``qsar/recipe.py`` excludes it from its default pool for that reason.
    """
    from chemsuite.analysis import fingerprints as fp
    from chemsuite.descriptors.compute import compute_descriptors
    from chemsuite.qsar.recipe import default_descriptors

    cells = list(data.smiles)
    keys = default_descriptors()
    columns = compute_descriptors(cells, keys)
    described = np.column_stack(
        [np.asarray(columns[key].values(), dtype=np.float64) for key in keys]
    )
    counts, _valid = fp.compute_counts(cells, fp.variant("morgan2", counts=True))
    bits = tuple(f"f:ecfp4:{i}" for i in range(counts.shape[1]))
    names = tuple(f"d:{key}" for key in keys) + bits
    return np.hstack([described, counts.astype(np.float64)]), names


def scaffold_folds(data: Data, rows: np.ndarray, folds: int, seed: int) -> list[np.ndarray]:
    """*rows* dealt into *folds* groups, a Bemis-Murcko scaffold never split across two.

    Largest group first into the emptiest fold, which keeps the folds even while keeping every
    scaffold whole; acyclic rows have no scaffold to share and are dealt one at a time.
    """
    from chemsuite.qsar.splits import scaffold_groups

    groups = scaffold_groups([data.mols[i] for i in rows])
    members: dict[str, list[int]] = {}
    loose: list[int] = []
    for position, group in enumerate(groups):
        if group:
            members.setdefault(str(group), []).append(int(rows[position]))
        else:
            loose.append(int(rows[position]))

    rng = np.random.default_rng(seed)
    parcels = [np.array(v) for v in members.values()] + [np.array([i]) for i in loose]
    rng.shuffle(parcels)
    parcels.sort(key=len, reverse=True)
    out: list[list[int]] = [[] for _ in range(folds)]
    for parcel in parcels:
        out[min(range(folds), key=lambda f: len(out[f]))] += list(parcel)
    return [np.sort(np.array(f, dtype=np.int64)) for f in out]


def scores_of(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    """R², RMSE, MAE and Spearman, on the rows that have both."""
    from scipy.stats import spearmanr

    keep = np.isfinite(observed) & np.isfinite(predicted)
    observed, predicted = observed[keep], predicted[keep]
    if observed.size < 3:
        return {"r2": float("nan"), "rmse": float("nan"), "mae": float("nan"), "rho": float("nan")}
    residual = observed - predicted
    spread = observed - observed.mean()
    return {
        "r2": float(1.0 - residual @ residual / (spread @ spread)),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mae": float(np.mean(np.abs(residual))),
        "rho": float(spearmanr(observed, predicted).statistic),
    }


def out_of_fold(x: np.ndarray, y: np.ndarray, folds: list[np.ndarray], seed: int) -> np.ndarray:
    """Predict every row from a model that never saw its fold."""
    from sklearn.ensemble import HistGradientBoostingRegressor

    predicted = np.full(y.shape, np.nan)
    for held in folds:
        trained = np.setdiff1d(np.concatenate(folds), held)
        model = HistGradientBoostingRegressor(
            max_iter=300, learning_rate=0.05, early_stopping=False, random_state=seed
        )
        model.fit(x[trained], y[trained])
        predicted[held] = model.predict(x[held])
    return predicted


def measure(data: Data, x: np.ndarray, isoform: str, seeds: int) -> dict[str, object]:
    """One isoform, over *seeds* scaffold-split repeats, with its shuffled-label null."""
    rows = data.rows_for(isoform, "train")
    y = data.values[isoform]
    real: list[dict[str, float]] = []
    null: list[dict[str, float]] = []
    for seed in range(seeds):
        folds = scaffold_folds(data, rows, FOLDS, seed)
        every = np.concatenate(folds)
        real.append(scores_of(y[every], out_of_fold(x, y, folds, seed)[every]))

        shuffled = y.copy()
        rng = np.random.default_rng(1000 + seed)
        shuffled[rows] = rng.permutation(y[rows])
        null.append(scores_of(shuffled[every], out_of_fold(x, shuffled, folds, seed)[every]))
    return {"n": int(rows.size), "real": real, "null": null}


def external(data: Data, x: np.ndarray, isoform: str, seed: int = 0) -> dict[str, float]:
    """Fit on every training row and score the unblinded test file (a sanity check, not blind)."""
    from sklearn.ensemble import HistGradientBoostingRegressor

    trained = data.rows_for(isoform, "train")
    held = data.rows_for(isoform, "test")
    y = data.values[isoform]
    model = HistGradientBoostingRegressor(
        max_iter=300, learning_rate=0.05, early_stopping=False, random_state=seed
    )
    model.fit(x[trained], y[trained])
    return scores_of(y[held], model.predict(x[held]))


def averaged(runs: list[dict[str, float]], key: str) -> tuple[float, float]:
    """Mean and spread of *key* over the repeats."""
    values = [run[key] for run in runs if np.isfinite(run[key])]
    if not values:
        return float("nan"), float("nan")
    return float(np.mean(values)), float(np.std(values))


def report(results: dict[str, dict[str, object]], outside: dict[str, dict[str, float]]) -> None:
    """Markdown tables, to paste into docs/SOM_PLAN.md."""
    print("\n**Scaffold-split, out of fold** (mean ± sd over seeds):\n")
    print("| isoform | n | R² | RMSE | MAE | Spearman |")
    print("|---|---|---|---|---|---|")
    for isoform, found in results.items():
        runs = found["real"]  # type: ignore[index]
        r2, r2sd = averaged(runs, "r2")
        rmse, _ = averaged(runs, "rmse")
        mae, _ = averaged(runs, "mae")
        rho, _ = averaged(runs, "rho")
        print(
            f"| {isoform} | {found['n']} | {r2:.3f} ± {r2sd:.3f} | "
            f"{rmse:.2f} | {mae:.2f} | {rho:.3f} |"
        )

    print("\n**Shuffled labels** (the same pipeline; anything above noise means a leak):\n")
    print("| isoform | R² | Spearman |")
    print("|---|---|---|")
    for isoform, found in results.items():
        runs = found["null"]  # type: ignore[index]
        r2, _ = averaged(runs, "r2")
        rho, _ = averaged(runs, "rho")
        print(f"| {isoform} | {r2:.3f} | {rho:.3f} |")

    print("\n**The unblinded 750-compound file** (public since 2026-11-04: a check, not a test):\n")
    print("| isoform | R² | RMSE | MAE | Spearman |")
    print("|---|---|---|---|---|")
    for isoform, found in outside.items():
        print(
            f"| {isoform} | {found['r2']:.3f} | {found['rmse']:.2f} | "
            f"{found['mae']:.2f} | {found['rho']:.3f} |"
        )


def through_engine(data: Data, x: np.ndarray, names: tuple[str, ...], isoform: str) -> None:
    """The same pool and target through chemsuite's Auto QSAR engine, for comparison."""
    from chemsuite.qsar import search
    from chemsuite.qsar.recipe import Recipe
    from chemsuite.qsar.splits import scaffold_groups

    rows = data.rows_for(isoform, "train")
    recipe = Recipe(
        task="regression",
        target=f"pIC50 {isoform}",
        structures="SMILES",
        descriptors=(),
        fingerprints=(),
        split="scaffold",
        budget="quick",
    )
    groups = np.full(len(data.smiles), "", dtype=object)
    groups[rows] = scaffold_groups([data.mols[i] for i in rows])
    started = time.perf_counter()
    outcome = search.run(x, data.values[isoform], rows, recipe, groups=groups)
    seconds = time.perf_counter() - started
    best = outcome.leaderboard[0] if outcome.leaderboard else None
    print(
        f"  {isoform}: {outcome.scores} in {seconds:.0f} s"
        + (f", best {best.candidate}" if best is not None else "")
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, default=INSTANCES)
    parser.add_argument("--seeds", type=int, default=SEEDS)
    parser.add_argument("--engine", action="store_true", help="also run Auto QSAR over the pool")
    arguments = parser.parse_args(argv[1:])

    if not arguments.instances.exists():
        print(f"{arguments.instances} is not there; run scripts/cyp/curate.py", file=sys.stderr)
        return 1

    data = read(arguments.instances)
    print(f"{len(data.smiles)} compounds; building the feature pool…")
    started = time.perf_counter()
    x, names = pool_for(data)
    print(f"  {x.shape[1]} features in {time.perf_counter() - started:.0f} s")

    results: dict[str, dict[str, object]] = {}
    outside: dict[str, dict[str, float]] = {}
    for isoform in ISOFORMS:
        rows = data.rows_for(isoform, "train")
        if rows.size < MIN_ROWS:
            print(f"  {isoform}: only {rows.size} rows, skipped")
            continue
        started = time.perf_counter()
        results[isoform] = measure(data, x, isoform, arguments.seeds)
        outside[isoform] = external(data, x, isoform)
        print(f"  {isoform}: {rows.size} rows, {time.perf_counter() - started:.0f} s")

    report(results, outside)

    if arguments.engine:
        print("\nThrough Auto QSAR (quick budget):")
        for isoform in results:
            through_engine(data, x, names, isoform)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
