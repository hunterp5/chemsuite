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

"""Convert our trained Uni-pKa folds to the numpy archive chemsuite reads (UNIPKA_PLAN.md, U3).

Runs in the training env (torch, and Uni-pKa's source in ~/bench/unipka/src), never in the app:

    python convert_weights.py <save dir with fold_*/checkpoint_best.pt> <train.lmdb> <out.npz>
        [--nfolds 5] [--cv-seed 42]

Uni-pKa centres pKa on the mean of each fold's training split, which upstream's own split code
recomputes here from *train.lmdb* (the same nfolds and seed as ``train.sh``). The archive keeps
them as ``fold_means`` and their average as ``mean``; the predictor moves every fold onto
``mean`` (``thermo.common_mean``). Each fold's inference parameters are stored as fp16 under
``f<k>/<name>``, without the ``unimol.`` prefix, ``k`` being the fold's number.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch

KEEP = ("embed_tokens.", "embed_charges.", "gbf.", "gbf_proj.", "encoder.emb_layer_norm.",
        "encoder.layers.", "encoder.final_layer_norm.", "classification_heads.")  # fmt: skip


def fold_means(lmdb_path: Path, nfolds: int, seed: int) -> list[float]:
    """Each fold's training-split mean pKa, by upstream's cross-validation split."""
    sys.path.insert(0, str(Path.home() / "bench/unipka/src"))
    from unicore.data import LMDBDataset
    from unimol.data.lmdb_dataset import FoldLMDBDataset

    data = LMDBDataset(str(lmdb_path))
    target = np.array([data[i]["target"] for i in range(len(data))], dtype=np.float64)
    keys, start, end = FoldLMDBDataset(data, seed, 0, nfolds=nfolds).get_fold_info()
    ordered = target[keys]
    return [float(np.delete(ordered, np.s_[start[k] : end[k]]).mean()) for k in range(nfolds)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("folder", type=Path)
    parser.add_argument("lmdb", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--nfolds", type=int, default=5)
    parser.add_argument("--cv-seed", type=int, default=42)
    args = parser.parse_args()
    every = fold_means(args.lmdb, args.nfolds, args.cv_seed)
    arrays: dict[str, np.ndarray] = {}
    means = []
    for path in sorted(args.folder.glob("fold_*/checkpoint_best.pt")):
        k = int(path.parent.name.removeprefix("fold_"))
        state = torch.load(path, map_location="cpu", weights_only=False)["model"]
        for name, value in state.items():
            short = name.removeprefix("unimol.")
            if short.startswith(KEEP):
                arrays[f"f{k}/{short}"] = value.float().numpy().astype(np.float16)
        means.append(every[k])
        print(path, "→", f"f{k}", f"mean {every[k]:.4f}")
    arrays["fold_means"] = np.array(means, dtype=np.float64)
    arrays["mean"] = np.array(np.mean(means), dtype=np.float64)
    np.savez_compressed(args.out, **arrays)
    print(args.out, f"{args.out.stat().st_size / 1e6:.1f} MB, mean {np.mean(means):.4f}")


if __name__ == "__main__":
    main()
