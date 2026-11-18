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
#
"""Measure and train chemsuite's sites-of-metabolism model. docs/SOM_PLAN.md, M1 and M2.

Runs in the app's own venv -- scikit-learn is already a dependency, a fit takes seconds,
and the pivot away from the graph network took training out of WSL and off the GPU.

**The features come from the app**, ``chemsuite.metabolism.features``, not from a copy
kept here: training and prediction then cannot drift apart, which is the usual way a
model that measured well ships badly.

**The input** is the curated instance table of the independent ``som-benchmark``
(Will de Snoo, MIT), whose curation carries the upstream 1-indexed SOM labels through
standardization into canonical SMILES atom order and excludes any record that fails a
round-trip check (see its ``src/som/data.py``). Build it there, once::

    python scripts/01_download.py && python scripts/02_curate.py

then copy ``data/processed/instances.csv`` and ``molecules.csv`` to ``samples/som/``.
The Zaretzki data itself never ships, and neither does this table; only the fitted trees
do.

**The labels are ``PRIMARY_SOM`` only** (decided with the user): the harder and more
defensible question, and the definition every number in docs/SOM_PLAN.md uses.

**The splits group by molecule**, never by molecule-isoform instance. The published
numbers in this field split instances, so 86% of test instances have their own molecule
in training under another isoform (``leakage_check.py``), which is worth about 0.18 MCC.
Only the grouped numbers are quoted.

Typical use::

    python scripts/som/train_gbm.py                       # measure, 10 seeds
    python scripts/som/train_gbm.py --export chemsuite/metabolism/data/som_gbm.npz
"""

from __future__ import annotations

import argparse
import csv
import statistics
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, matthews_corrcoef

from chemsuite.metabolism import ISOFORMS
from chemsuite.metabolism.features import (
    DOMAIN_RADII,
    FEATURE_NAMES,
    atom_block,
    atom_environments,
    isoform_one_hot,
)

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_INSTANCES = REPO / "samples" / "som" / "instances.csv"

#: ``som-benchmark``'s grid, and its protocol: the four configurations are scored on the
#: validation split of the first seed by average precision, and the choice is reused.
GRID = (
    {"learning_rate": 0.05, "max_leaf_nodes": 31},
    {"learning_rate": 0.05, "max_leaf_nodes": 63},
    {"learning_rate": 0.10, "max_leaf_nodes": 31},
    {"learning_rate": 0.10, "max_leaf_nodes": 63},
)
#: Fixed across the grid, from ``som-benchmark``'s ``_fit_gbm``.
FIXED = {
    "max_iter": 300,
    "early_stopping": False,
    "l2_regularization": 1.0,
    "min_samples_leaf": 20,
}
TEST_FRACTION = 0.2
VALIDATION_FRACTION = 0.1


# ---------------------------------------------------------------------------------- data


@dataclass(frozen=True)
class Instance:
    """One (molecule, isoform) measurement: its primary sites, in SMILES atom order."""

    isoform: str
    smiles: str
    sites: tuple[int, ...]


@dataclass
class Table:
    """The flat atom-level design matrix, with the bookkeeping to group it again."""

    x: np.ndarray
    y: np.ndarray
    instance: np.ndarray  # which Instance each row belongs to

    def rows_of(self, instances: np.ndarray) -> np.ndarray:
        return np.isin(self.instance, instances)


def read_instances(path: Path) -> list[Instance]:
    """The curated table, primary sites only."""
    out = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["isoform"] not in ISOFORMS:
                raise ValueError(f"unknown isoform in {path}: {row['isoform']!r}")
            out.append(
                Instance(
                    isoform=row["isoform"],
                    smiles=row["smiles"],
                    sites=tuple(int(index) for index in row["som_primary"].split()),
                )
            )
    return out


def build_table(instances: list[Instance]) -> Table:
    """Featurize every atom of every instance, each molecule featurized once."""
    blocks: dict[str, np.ndarray] = {}
    rows, labels, instance_of = [], [], []
    for position, instance in enumerate(instances):
        if instance.smiles not in blocks:
            mol = Chem.MolFromSmiles(instance.smiles)
            if mol is None:
                raise ValueError(f"curated SMILES failed to parse: {instance.smiles}")
            blocks[instance.smiles] = atom_block(mol)
        block = blocks[instance.smiles]
        atoms = block.shape[0]
        if instance.sites and max(instance.sites) >= atoms:
            raise ValueError(f"site index past the last atom: {instance}")

        one_hot = np.tile(isoform_one_hot(instance.isoform), (atoms, 1))
        rows.append(np.hstack([block, one_hot]))
        label = np.zeros(atoms, dtype=np.int8)
        label[list(instance.sites)] = 1
        labels.append(label)
        instance_of.append(np.full(atoms, position, dtype=np.int32))

    table = Table(
        x=np.vstack(rows),
        y=np.concatenate(labels),
        instance=np.concatenate(instance_of),
    )
    if table.x.shape[1] != len(FEATURE_NAMES):
        raise ValueError("the design matrix does not match FEATURE_NAMES")
    return table


def split(instances: list[Instance], seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Partition molecules, then map the partition onto instance indices.

    Test is drawn first and validation out of what is left, so tuning never sees a test
    molecule under any isoform.
    """
    rng = np.random.default_rng(seed)
    molecules = sorted({instance.smiles for instance in instances})
    rng.shuffle(molecules)
    test_size = int(round(TEST_FRACTION * len(molecules)))
    validation_size = int(round(VALIDATION_FRACTION * len(molecules)))
    test = set(molecules[:test_size])
    validation = set(molecules[test_size : test_size + validation_size])

    train_i, validation_i, test_i = [], [], []
    for position, instance in enumerate(instances):
        if instance.smiles in test:
            test_i.append(position)
        elif instance.smiles in validation:
            validation_i.append(position)
        else:
            train_i.append(position)
    return (
        np.asarray(train_i, dtype=np.int64),
        np.asarray(validation_i, dtype=np.int64),
        np.asarray(test_i, dtype=np.int64),
    )


# ------------------------------------------------------------------------------- metrics


def top_k_hits(scores: np.ndarray, table: Table, instances: np.ndarray, k: int) -> list[bool]:
    """For each instance, whether any true site is among its *k* highest-scoring atoms.

    Ties are broken by lowest atom index, exactly as the app breaks them, so a molecule
    with symmetry-equivalent atoms is scored the way it will be displayed.
    """
    hits = []
    for identifier in instances:
        rows = np.flatnonzero(table.instance == identifier)
        if not rows.size:
            continue
        atom_scores = scores[rows]
        ranked = np.lexsort((np.arange(atom_scores.size), -atom_scores))[:k]
        hits.append(bool(table.y[rows][ranked].any()))
    return hits


def evaluate(scores: np.ndarray, table: Table, instances: np.ndarray) -> dict[str, float]:
    """MCC over atoms at a 0.5 cut, and top-1/2/3 over instances."""
    rows = table.rows_of(instances)
    truth, predicted = table.y[rows], (scores >= 0.5).astype(np.int8)
    result = {"mcc": float(matthews_corrcoef(truth, predicted))}
    for k in (1, 2, 3):
        hits = top_k_hits(scores_by_row(scores, table, rows), table, instances, k)
        result[f"top{k}"] = float(np.mean(hits)) if hits else float("nan")
    return result


def scores_by_row(scores: np.ndarray, table: Table, rows: np.ndarray) -> np.ndarray:
    """*scores*, which cover only *rows*, scattered back over every row of the table."""
    full = np.full(table.y.size, -np.inf)
    full[rows] = scores
    return full


def fit(table: Table, train: np.ndarray, parameters: dict, seed: int):
    """One fit on the instances *train*."""
    rows = table.rows_of(train)
    model = HistGradientBoostingClassifier(**FIXED, random_state=seed, **parameters)
    model.fit(table.x[rows], table.y[rows])
    return model


def tune(table: Table, train: np.ndarray, validation: np.ndarray) -> dict:
    """The grid configuration with the best validation average precision."""
    validation_rows = table.rows_of(validation)
    best, best_score = GRID[0], -np.inf
    for parameters in GRID:
        model = fit(table, train, parameters, seed=0)
        probability = model.predict_proba(table.x[validation_rows])[:, 1]
        score = average_precision_score(table.y[validation_rows], probability)
        print(f"  {parameters}  validation AP {score:.4f}")
        if score > best_score:
            best, best_score = parameters, score
    return dict(best)


# ----------------------------------------------------------------------------------- run


def measure(table: Table, instances: list[Instance], seeds: int) -> dict:
    """Grouped-split accuracy over *seeds* seeds, overall and per isoform."""
    train, validation, _ = split(instances, 0)
    print("tuning on the validation split of seed 0:")
    parameters = tune(table, train, validation)
    print(f"chosen: {parameters}\n")

    overall: list[dict[str, float]] = []
    per_isoform: dict[str, list[dict[str, float]]] = {name: [] for name in ISOFORMS}
    for seed in range(seeds):
        train, _, test = split(instances, seed)
        model = fit(table, train, parameters, seed)
        rows = table.rows_of(test)
        scores = model.predict_proba(table.x[rows])[:, 1]
        overall.append(evaluate(scores, table, test))
        for name in ISOFORMS:
            subset = np.asarray([i for i in test if instances[i].isoform == name], dtype=np.int64)
            if subset.size:
                subset_rows = table.rows_of(subset)
                keep = subset_rows[rows]  # the test rows that are this isoform's
                per_isoform[name].append(evaluate(scores[keep], table, subset))
        print(
            f"seed {seed}: MCC {overall[-1]['mcc']:.3f}  "
            f"top-1 {overall[-1]['top1']:.3f}  top-2 {overall[-1]['top2']:.3f}  "
            f"top-3 {overall[-1]['top3']:.3f}"
        )
    return {"parameters": parameters, "overall": overall, "per_isoform": per_isoform}


def summarize(runs: list[dict[str, float]]) -> dict[str, tuple[float, float]]:
    """Mean and standard deviation of each metric over the seeds."""
    out = {}
    for metric in ("mcc", "top1", "top2", "top3"):
        values = [run[metric] for run in runs]
        deviation = statistics.stdev(values) if len(values) > 1 else 0.0
        out[metric] = (statistics.fmean(values), deviation)
    return out


def report(result: dict) -> None:
    summary = summarize(result["overall"])
    print(f"\nmolecule-grouped, PRIMARY_SOM only, over {len(result['overall'])} seeds:")
    print("| metric | mean | sd |\n|---|---|---|")
    for metric, (mean, deviation) in summary.items():
        print(f"| {metric} | {mean:.3f} | {deviation:.3f} |")

    print("\n| isoform | MCC | top-1 | top-2 | top-3 |\n|---|---|---|---|---|")
    for name, runs in result["per_isoform"].items():
        if not runs:
            continue
        values = summarize(runs)
        print(
            f"| {name} | {values['mcc'][0]:.3f} | {values['top1'][0]:.3f} | "
            f"{values['top2'][0]:.3f} | {values['top3'][0]:.3f} |"
        )


def domain_index(instances: list[Instance]) -> dict[str, np.ndarray]:
    """The applicability-domain index: every training atom's Morgan environments.

    One sorted array of unique identifiers per radius. An atom of a new molecule whose
    environments are all in the index is one the model has seen the like of; a
    radius-2 environment nobody has seen is the useful warning.
    """
    seen: dict[int, set[int]] = {radius: set() for radius in DOMAIN_RADII}
    for smiles in {instance.smiles for instance in instances}:
        mol = Chem.MolFromSmiles(smiles)
        for pairs in atom_environments(mol):
            for radius, identifier in pairs:
                seen[radius].add(identifier)
    return {
        f"domain_r{radius}": np.asarray(sorted(values), dtype=np.int64)
        for radius, values in seen.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, default=DEFAULT_INSTANCES)
    parser.add_argument("--seeds", type=int, default=10)
    parser.add_argument(
        "--export",
        type=Path,
        help="fit on every instance with the tuned settings and write the .npz",
    )
    arguments = parser.parse_args()

    instances = read_instances(arguments.instances)
    molecules = len({instance.smiles for instance in instances})
    print(f"{len(instances)} instances, {molecules} molecules, {arguments.instances}")
    table = build_table(instances)
    print(
        f"{table.x.shape[0]} atoms, {table.x.shape[1]} columns, "
        f"{int(table.y.sum())} labelled sites "
        f"({table.y.mean():.4f} of atoms)\n"
    )

    result = measure(table, instances, arguments.seeds)
    report(result)

    if arguments.export:
        from export_gbm import write

        print("\nfitting on every instance with the tuned settings")
        every = np.arange(len(instances), dtype=np.int64)
        model = fit(table, every, result["parameters"], seed=0)
        extra = {
            "isoforms": np.asarray(ISOFORMS),
            "feature_names": np.asarray(FEATURE_NAMES),
            **domain_index(instances),
        }
        written = write(arguments.export, model, check_on=table.x, extra=extra)
        print(
            f"wrote {arguments.export}: {int(written['trees'])} trees, "
            f"{int(written['nodes'])} nodes, {written['megabytes']:.2f} MB; "
            f"largest difference from scikit-learn "
            f"{written['probability_difference']:.2g}"
        )


if __name__ == "__main__":
    main()
