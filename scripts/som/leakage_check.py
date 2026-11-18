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
"""Measure how much of upstream's cross-validation leaks the same molecule across isoforms.

docs/SOM_PLAN.md, M1. Upstream splits **molecule-isoform instances**:

    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    for trainval_idx, test_idx in skf.split(np.arange(len(graphs)), y_strat):

A molecule appears once per isoform it was measured against (2,003 instances over 679
molecules), so the same structure can be trained on for 1A2 and tested on for 3A4.
``apply_no_leakage_to_dataloaders`` masks only the auxiliary attention targets and says so:
"The atom-level labels ``y`` are untouched."

The paper's own premise is that SoM annotations are highly similar across isoform pairs, so
this inflates the reported top-k and MCC as estimates of performance on a *new* molecule.
This script quantifies it, exactly reproducing upstream's split. Nothing here ships.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", required=True, type=Path, help="the ATTNSOM clone")
    parser.add_argument("--dataset-dir", default=None, help="defaults to <upstream>/cyp_dataset")
    parser.add_argument("--n-splits", default=10, type=int)
    parser.add_argument("--seed", default=42, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(args.upstream))

    from dataset import load_multi_cyp
    from dataset_utils import make_strata_labels
    from sklearn.model_selection import StratifiedKFold

    cyps = ["1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"]
    base = args.dataset_dir or str(args.upstream / "cyp_dataset")
    graphs = load_multi_cyp(base, cyps)
    smiles = [g.smiles for g in graphs]

    per_molecule = Counter(smiles)
    print(f"instances {len(graphs)}, unique molecules {len(per_molecule)}")
    spread = Counter(per_molecule.values())
    print("isoforms per molecule:", dict(sorted(spread.items())))
    singles = sum(1 for s in smiles if per_molecule[s] == 1)
    print(f"instances whose molecule appears once only: {singles} ({singles / len(graphs):.1%})")
    print()

    strata = make_strata_labels(graphs)
    splitter = StratifiedKFold(n_splits=args.n_splits, shuffle=True, random_state=args.seed)

    leaked_total = 0
    tested_total = 0
    print(f"{'fold':>4} {'test':>5} {'leaked':>7} {'share':>7}")
    for fold, (trainval_idx, test_idx) in enumerate(
        splitter.split(np.arange(len(graphs)), strata), start=1
    ):
        train_smiles = {smiles[i] for i in trainval_idx}
        leaked = sum(1 for i in test_idx if smiles[i] in train_smiles)
        leaked_total += leaked
        tested_total += len(test_idx)
        print(f"{fold:>4} {len(test_idx):>5} {leaked:>7} {leaked / len(test_idx):>6.1%}")

    print()
    print(
        f"overall: {leaked_total} of {tested_total} test instances "
        f"({leaked_total / tested_total:.1%}) share their molecule with the training split "
        f"under a different isoform"
    )
    print()
    print("A molecule-grouped split is the honest comparison; see docs/SOM_PLAN.md M1.")


if __name__ == "__main__":
    main()
