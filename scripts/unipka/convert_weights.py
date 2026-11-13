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

Runs in the training env (torch), never in the app:

    python convert_weights.py <save dir with fold_*/checkpoint_best.pt> <out.npz> <mean pKa>

*mean pKa* is the training set's mean target, which Uni-pKa subtracts while training; the
archive keeps it as ``mean``. Each fold's inference parameters are stored as fp16 under
``f<k>/<name>``, without the ``unimol.`` prefix.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import torch

KEEP = ("embed_tokens.", "embed_charges.", "gbf.", "gbf_proj.", "encoder.emb_layer_norm.",
        "encoder.layers.", "encoder.final_layer_norm.", "classification_heads.")  # fmt: skip


def main() -> None:
    folder, out, mean = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
    arrays = {"mean": np.array(mean, dtype=np.float64)}
    for k, path in enumerate(sorted(folder.glob("fold_*/checkpoint_best.pt"))):
        state = torch.load(path, map_location="cpu", weights_only=False)["model"]
        for name, value in state.items():
            short = name.removeprefix("unimol.")
            if short.startswith(KEEP):
                arrays[f"f{k}/{short}"] = value.float().numpy().astype(np.float16)
        print(path, "→", f"f{k}")
    np.savez_compressed(out, **arrays)
    print(out, f"{out.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
