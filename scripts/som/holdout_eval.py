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
"""Does averaging the folds buy anything? Measured on a truly external set.

docs/SOM_PLAN.md, M1 step 3. This decides how the weights ship (M3): one fold is 34 MB, so
float16 it is about 17 MB and could be vendored in-tree beside MolGpKa's files, while a ten-fold
ensemble is ~341 MB and has to be a pinned download with a consent prompt. The pKa work found
five folds bought only 0.01-0.02 RMSE over one, so this is worth measuring before building any
download machinery.

A k-fold run cannot answer it: for any fold's test molecules, every *other* fold trained on
them, so the ensemble cannot be scored without leakage. AZ-ExactSOM is the way out -- a separate
120-compound set from human hepatocyte assays that the Zaretzki-trained folds never saw. Any
molecule it shares with Zaretzki is dropped, and the overlap is reported.

AZ carries no isoform labels (its sites are aggregate hepatocyte outcomes), so a Zaretzki model
is asked for an isoform-agnostic answer two ways: the 3A4 head alone (3A4 dominates hepatic
clearance) and the per-atom maximum over all nine heads ("metabolized by something").

Run in WSL, in the training env. Nothing here ships.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", required=True, type=Path, help="the ATTNSOM clone")
    parser.add_argument("--ckpt-dir", required=True, type=Path, help="holds fold*_best.pt")
    parser.add_argument("--zaretzki-dir", default=None, help="to detect overlap; default clone")
    parser.add_argument("--az-csv", default=None, help="default <clone>/cyp_dataset/az/...")
    return parser.parse_args()


def metrics_for(trues: list[np.ndarray], probs: list[np.ndarray]) -> dict[str, float]:
    """Top-k over ranked atoms, plus MCC at 0.5. Molecules with no labelled site are skipped."""
    out: dict[str, float] = {}
    for k in (1, 2, 3):
        hits = 0
        counted = 0
        for true, prob in zip(trues, probs, strict=True):
            if true.sum() == 0:
                continue
            counted += 1
            top = np.argsort(-prob)[:k]
            if true[top].sum() > 0:
                hits += 1
        out[f"top{k}"] = hits / max(counted, 1)
    flat_true = np.concatenate([t for t in trues if t.sum() > 0])
    flat_prob = np.concatenate([p for t, p in zip(trues, probs, strict=True) if t.sum() > 0])
    pred = (flat_prob >= 0.5).astype(int)
    tp = float(((pred == 1) & (flat_true == 1)).sum())
    tn = float(((pred == 0) & (flat_true == 0)).sum())
    fp = float(((pred == 1) & (flat_true == 0)).sum())
    fn = float(((pred == 0) & (flat_true == 1)).sum())
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out["mcc"] = (tp * tn - fp * fn) / denom if denom > 0 else 0.0
    out["n_molecules"] = float(sum(1 for t in trues if t.sum() > 0))
    return out


def external_set(modules, cyps: list[str], az_csv: str, zaretzki_dir: Path) -> list[dict]:
    """AZ-ExactSOM compounds with exact sites that no Zaretzki fold can have trained on."""
    chem = modules["Chem"]
    trained_on: set[str] = set()
    for iso in cyps:
        for record in modules["load_cyp_sdf"](str(zaretzki_dir / f"{iso}.sdf")):
            trained_on.add(chem.MolToSmiles(record["mol"]))

    az = modules["load_az_csv"](az_csv)
    kept = []
    overlapping = 0
    for record in az:
        if chem.MolToSmiles(record["mol"]) in trained_on:
            overlapping += 1
            continue
        if record["som_idxs"]:
            kept.append(record)
    print(
        f"AZ compounds: {len(az)}; dropped for Zaretzki overlap: {overlapping}; "
        f"with exact sites and unseen: {len(kept)}"
    )
    return kept


def fold_predictions(modules, cyps, kept, ckpt_dir: Path) -> list[list[np.ndarray]]:
    """Per-atom probabilities for every fold: probs[fold][molecule] -> (isoforms, atoms)."""
    chem = modules["Chem"]
    config = json.loads((ckpt_dir / "config.json").read_text())
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    paths = sorted(
        ckpt_dir.glob("fold*_best.pt"),
        key=lambda p: int("".join(c for c in p.name if c.isdigit())),
    )
    print(f"{len(paths)} checkpoints from {ckpt_dir}")

    per_fold: list[list[np.ndarray]] = []
    for path in paths:
        model = modules["ATTNSOM"](**config).to(device)
        model.load_state_dict(torch.load(path, map_location=device, weights_only=True))
        model.eval()

        fold_probs: list[np.ndarray] = []
        for record in kept:
            mol = record["mol"]
            smiles = chem.MolToSmiles(mol)
            graphs = []
            for cyp_idx in range(len(cyps)):
                data, _ = modules["mol_to_graph"](mol, smiles, som_indices=[])
                data.cyp_idx = torch.full((data.num_nodes,), cyp_idx, dtype=torch.long)
                graphs.append(data)
            # All nine isoforms in one batch: the SAGPool branch, as evaluated upstream.
            with torch.no_grad():
                logits, _, _ = model(modules["Batch"].from_data_list(graphs).to(device))
            fold_probs.append(
                torch.sigmoid(logits).view(len(cyps), mol.GetNumAtoms()).cpu().numpy()
            )
        per_fold.append(fold_probs)
    return per_fold


def labels_of(kept: list[dict]) -> list[np.ndarray]:
    trues = []
    for record in kept:
        label = np.zeros(record["mol"].GetNumAtoms())
        for i in record["som_idxs"]:
            if 0 <= i < len(label):
                label[i] = 1.0
        trues.append(label)
    return trues


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(args.upstream))

    from dataset import load_az_csv, load_cyp_sdf
    from dataset_utils import mol_to_graph
    from model import ATTNSOM
    from rdkit import Chem, RDLogger
    from torch_geometric.data import Batch

    RDLogger.DisableLog("rdApp.*")
    modules = {
        "load_az_csv": load_az_csv,
        "load_cyp_sdf": load_cyp_sdf,
        "mol_to_graph": mol_to_graph,
        "ATTNSOM": ATTNSOM,
        "Chem": Chem,
        "Batch": Batch,
    }

    cyps = ["1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"]
    clone = args.upstream
    az_csv = args.az_csv or str(clone / "cyp_dataset" / "az" / "az_120_compounds.csv")
    zaretzki_dir = Path(args.zaretzki_dir or (clone / "cyp_dataset"))

    kept = external_set(modules, cyps, az_csv, zaretzki_dir)
    per_fold = fold_predictions(modules, cyps, kept, args.ckpt_dir)
    trues = labels_of(kept)

    for view, pick in (
        ("3A4 head", lambda a: a[cyps.index("3A4")]),
        ("max over heads", lambda a: a.max(axis=0)),
    ):
        singles = [metrics_for(trues, [pick(a) for a in fold]) for fold in per_fold]
        mean_probs = [
            pick(np.mean([per_fold[f][m] for f in range(len(per_fold))], axis=0))
            for m in range(len(kept))
        ]
        ensemble = metrics_for(trues, mean_probs)

        print()
        print(f"=== {view} ===  ({int(ensemble['n_molecules'])} molecules scored)")
        print(f"{'metric':>8} {'single fold (mean ± std)':>26} {'ensemble':>10} {'gain':>8}")
        for key in ("top1", "top2", "top3", "mcc"):
            values = np.array([s[key] for s in singles])
            gain = ensemble[key] - values.mean()
            print(
                f"{key:>8} {values.mean():>14.3f} ± {values.std():<9.3f} "
                f"{ensemble[key]:>10.3f} {gain:>+8.3f}"
            )

    # How many folds are actually worth shipping? Average random subsets of each size, so the
    # curve shows where the gain flattens against the download it costs (34.15 MB per fold).
    rng = np.random.default_rng(0)
    n_folds = len(per_fold)
    print()
    print("=== ensemble size sweep (max over heads) ===")
    print(f"{'folds':>5} {'float32':>9} {'float16':>9} {'top1':>14} {'top2':>14} {'mcc':>14}")
    for size in (1, 2, 3, 5, 7, n_folds):
        if size > n_folds:
            continue
        draws = 1 if size == n_folds else 12
        rows = []
        for _ in range(draws):
            chosen = rng.choice(n_folds, size=size, replace=False)
            probs = [
                np.mean([per_fold[f][m] for f in chosen], axis=0).max(axis=0)
                for m in range(len(kept))
            ]
            rows.append(metrics_for(trues, probs))
        summary = {
            key: (np.mean([r[key] for r in rows]), np.std([r[key] for r in rows]))
            for key in ("top1", "top2", "mcc")
        }
        mb32 = size * 34.15
        print(
            f"{size:>5} {mb32:>8.0f}M {mb32 / 2:>8.0f}M "
            + " ".join(
                f"{summary[key][0]:>8.3f} ± {summary[key][1]:<3.3f}"
                for key in ("top1", "top2", "mcc")
            )
        )

    print()
    print("A gain inside the single-fold spread means one fold is enough: ship one fold in-tree")
    print("and drop M3's download, consent and cache path (docs/SOM_PLAN.md M3).")


if __name__ == "__main__":
    main()
