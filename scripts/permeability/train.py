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

"""Fit the permeability models, one arm at a time, and write the ones to ship (P2).

    .venv/Scripts/python.exe scripts/permeability/train.py --arm plain
    .venv/Scripts/python.exe scripts/permeability/train.py --arm pka --endpoint caco2
    .venv/Scripts/python.exe scripts/permeability/train.py --ship plain

Reads the CSVs ``curate.py`` wrote and fits through **chemsuite's own Auto QSAR engine** —
``qsar.features.build`` for the pool and ``qsar.fit.fit`` for the search — so the shipped model
is the same kind of object the Data ▸ QSAR tool produces and the app applies it with no
scikit-learn at run time. Nothing here reimplements modelling.

Every number printed is the engine's own **out-of-fold** score (``Fitted.scores``): what the
model would have said about each training compound had it never seen it, under the recipe's
scaffold split. That is the only honest figure available without a second holdout, and it is
the figure docs/PERMEABILITY_PLAN.md quotes.

**The four arms**, all on the same split and seed, so the only difference is the feature pool:

* ``plain`` — every descriptor that needs no conformer and no model, plus Morgan count
  fingerprints. The baseline, and what the evidence says should win.
* ``pka`` — ``plain`` plus the pKa-derived descriptors (``cLogD 7.4`` and friends). These cost
  a pKa model run per molecule, seconds rather than milliseconds, so they ship only if they
  beat ``plain`` by :data:`PKA_MARGIN`. **Slow: hours for the larger endpoints.**
* ``descriptors`` — ``plain`` without the fingerprint block, to see whether 2,048 sparse
  columns earn their width on sets this size.
* ``shape`` — ``plain`` plus the 3D shape descriptors, over one ETKDG conformer per molecule.
  This is the cheap proxy for the question a Uni-Mol head would answer expensively, and its
  result is the gate on that work (PERMEABILITY_PLAN, last section).

Each run writes ``<endpoint>.<arm>.json`` beside the data so the arms can be run separately,
in any order, and compared later without refitting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from chemsuite.descriptors.catalog import CATALOG, PKA_KEYS, SHAPE_GROUP
from chemsuite.qsar import features
from chemsuite.qsar.fit import Fitted, fit, scaffold_groups
from chemsuite.qsar.recipe import BUDGETS, Recipe, default_descriptors, default_fingerprints
from chemsuite.qsar.splits import as_groups

REPO = Path(__file__).resolve().parents[2]
CURATED = REPO / "samples" / "permeability"
MODELS = REPO / "chemsuite" / "permeability" / "data"

#: ``pka`` ships only if it beats ``plain`` by this much MAE (log units). Below it, seconds per
#: molecule buys nothing a chemist would notice against an assay that scatters by 0.3-0.5.
PKA_MARGIN = 0.05

#: Which curated file each endpoint reads, its task, and what to call the shipped model.
ENDPOINTS = {
    "mdck_er": ("regression", "MDR1-MDCK efflux ratio"),
    "caco2": ("regression", "Caco-2 permeability"),
    "pampa": ("classification", "PAMPA permeability"),
}


@dataclass(frozen=True)
class Arm:
    """One feature pool to try."""

    key: str
    #: extra catalog keys beyond the plain set
    extra: tuple[str, ...] = ()
    fingerprints: bool = True
    #: whether the rows need a 3D conformer before the pool is built
    needs_3d: bool = False


def shape_keys() -> tuple[str, ...]:
    """The 3D shape descriptors, which need a conformer."""
    return tuple(key for key, d in CATALOG.items() if d.group == SHAPE_GROUP and not d.is_text)


def arms() -> dict[str, Arm]:
    return {
        "plain": Arm("plain"),
        "pka": Arm("pka", extra=tuple(sorted(PKA_KEYS))),
        "descriptors": Arm("descriptors", fingerprints=False),
        "shape": Arm("shape", extra=shape_keys(), needs_3d=True),
    }


def recipe_for(endpoint: str, arm: Arm, seed: int, budget: str) -> Recipe:
    """The recipe for one arm. Everything but the feature pool is held fixed."""
    task, _title = ENDPOINTS[endpoint]
    return Recipe(
        task=task,
        target=endpoint,
        budget=budget,
        # the name of the column a *caller* reads structures from. A recipe will not validate
        # without one, and for a shipped model it is documentation only: the tool hands the
        # pool the structure column the user chose, not this name.
        structures="Structure",
        descriptors=default_descriptors() + arm.extra,
        fingerprints=default_fingerprints() if arm.fingerprints else (),
        split="scaffold",
        seed=seed,
    )


def read_curated(path: Path, endpoint: str) -> tuple[list[str], np.ndarray]:
    """The SMILES and the values from a curated CSV, in file order."""
    smiles: list[str] = []
    values: list[float] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            smiles.append(row["smiles"])
            values.append(float(row[endpoint]))
    return smiles, np.asarray(values, dtype=np.float64)


def with_conformers(smiles: list[str]) -> list[bytes | str]:
    """One ETKDG conformer per molecule, as RDKit blobs; the SMILES text where embedding fails.

    The shape descriptors are the only thing that needs this, and a molecule that will not
    embed simply leaves them NaN — which the pool already handles as a missing feature.
    """
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem

    cells: list[bytes | str] = []
    with rdBase.BlockLogs():
        for text in smiles:
            mol = Chem.MolFromSmiles(text)
            if mol is None:
                cells.append(text)
                continue
            mol = Chem.AddHs(mol)
            parameters = AllChem.ETKDGv3()
            parameters.randomSeed = 0xF00D
            if AllChem.EmbedMolecule(mol, parameters) != 0:
                cells.append(text)
                continue
            cells.append(Chem.RemoveHs(mol).ToBinary())
    return cells


def classes_for(task: str, values: np.ndarray) -> tuple[np.ndarray | None, np.ndarray]:
    """The class labels and the target the engine wants: indices for a classification."""
    if task != "classification":
        return None, values
    labels = np.array([f"{value:g}" for value in sorted(set(values.tolist()))], dtype=np.str_)
    index = {name: i for i, name in enumerate(labels.tolist())}
    return labels, np.array([index[f"{value:g}"] for value in values], dtype=np.float64)


def run_arm(endpoint: str, arm: Arm, seed: int, workers: int, budget: str) -> tuple[Fitted, dict]:
    """Fit one arm of one endpoint and come back with the fit and what to record."""
    task, title = ENDPOINTS[endpoint]
    smiles, values = read_curated(CURATED / f"{endpoint}.csv", endpoint)
    cells: list[bytes | str] = list(smiles)
    if arm.needs_3d:
        print(f"  embedding {len(smiles)} conformers…", flush=True)
        cells = with_conformers(smiles)

    recipe = recipe_for(endpoint, arm, seed, budget)
    labels, y = classes_for(task, values)
    started = time.monotonic()
    print(f"  building the pool ({len(recipe.descriptors)} descriptors)…", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pool = features.build(cells, recipe, executor=executor)
    usable = pool.usable & np.isfinite(y)
    print(
        f"  pool {pool.rows} x {pool.columns}, {int(usable.sum())} usable; searching…",
        flush=True,
    )
    fitted = fit(
        pool.x,
        pool.names,
        np.nan_to_num(y, nan=0.0),
        np.flatnonzero(usable),
        recipe,
        keys=features.keys_for(cells),
        classes=labels,
        groups=as_groups(scaffold_groups(features.mols_for(cells))),
        label=f"{title} ({arm.key})",
    )
    record = {
        "endpoint": endpoint,
        "arm": arm.key,
        "task": task,
        "seed": seed,
        "budget": budget,
        "rows": int(usable.sum()),
        "pool_columns": pool.columns,
        "features_kept": fitted.model.features,
        "candidate": fitted.model.candidate,
        "scores": fitted.scores,
        "minutes": round((time.monotonic() - started) / 60, 1),
    }
    return fitted, record


def report(record: dict) -> str:
    """One line per arm, with the metrics that matter for its task."""
    scores = record["scores"]
    if record["task"] == "regression":
        shown = " ".join(
            f"{name} {scores[name]:.3f}" for name in ("mae", "rmse", "r2") if name in scores
        )
    else:
        shown = " ".join(
            f"{name} {scores[name]:.3f}"
            for name in ("roc_auc", "balanced_accuracy", "mcc")
            if name in scores
        )
    return (
        f"  {record['arm']:12s} {record['candidate']:22s} {shown}"
        f"  [{record['features_kept']}/{record['pool_columns']} features,"
        f" {record['rows']} rows, {record['minutes']} min]"
    )


def main(argv: list[str]) -> int:
    known = arms()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", action="append", default=[], choices=sorted(known))
    parser.add_argument("--endpoint", action="append", default=[], choices=sorted(ENDPOINTS))
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--workers", type=int, default=0, help="0 for the default pool size")
    parser.add_argument(
        "--budget",
        default="",
        choices=["", *sorted(BUDGETS)],
        help="how much search to spend; the default is quick to compare arms and standard to "
        "ship, because every arm of a comparison has to be given the same budget",
    )
    parser.add_argument(
        "--ship",
        default="",
        choices=sorted(known),
        help="write this arm's models into the package's data folder",
    )
    parser.add_argument("--into", type=Path, default=MODELS)
    arguments = parser.parse_args(argv[1:])

    asked = arguments.arm or ([arguments.ship] if arguments.ship else [])
    chosen_arms = [known[key] for key in asked]
    if not chosen_arms:
        print("choose at least one --arm, or --ship an arm", file=sys.stderr)
        return 2
    endpoints = arguments.endpoint or list(ENDPOINTS)
    workers = arguments.workers or None
    budget = arguments.budget or ("standard" if arguments.ship else "quick")

    for endpoint in endpoints:
        print(f"\n{ENDPOINTS[endpoint][1]} [{endpoint}], {budget} budget")
        for arm in chosen_arms:
            fitted, record = run_arm(endpoint, arm, arguments.seed, workers, budget)
            (CURATED / f"{endpoint}.{arm.key}.{budget}.json").write_text(
                json.dumps(record, indent=2), encoding="utf-8"
            )
            print(report(record))
            if arguments.ship == arm.key:
                arguments.into.mkdir(parents=True, exist_ok=True)
                path = fitted.model.save(arguments.into / f"{endpoint}")
                print(f"  shipped {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
                print(pin_line(endpoint, path))
    return 0


def pin_line(endpoint: str, path: Path) -> str:
    """The line to paste into ``PINS`` in ``chemsuite/permeability/models.py``.

    Printed rather than written: the pins live in Python source a reviewer reads, as the pKa
    weights' do, so updating one is a deliberate edit and never a side effect of a fit.
    """
    data = path.read_bytes()
    return (
        f'    "{endpoint}": Pin({len(data)}, '
        f'"{hashlib.sha256(data).hexdigest()}"),  # paste into models.PINS'
    )


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
