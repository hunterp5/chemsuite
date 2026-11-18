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
"""Convert a trained site-of-metabolism checkpoint to the .npz the app ships.

docs/SOM_PLAN.md, M1/M3. The app has no torch: it runs the forward pass in numpy from this
archive, which is vendored in `chemsuite/metabolism/data/` (one fold, about 17 MB, decided
2026-11-18) and loaded with importlib.resources, as `ionisation/molgpka.py` loads its own.

Keys keep their PyTorch parameter names, so the numpy port reads like the paper and a diff
against upstream's `model.py` stays possible. Large matrices are stored float16 to halve the
file; small ones (norm gains, Chebyshev coefficients, the GIN epsilon, group scales) stay
float32, because they cost nothing and precision there is not worth risking.

``--check`` re-loads the archive into the torch model and compares logits against the
float32 original on real molecules, which is what justifies float16 at all.

Run in WSL, in the training env.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

# Anything smaller than this stays float32: the total is a few hundred kilobytes.
FLOAT32_BELOW = 65_536


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ckpt", required=True, type=Path, help="a fold*_best.pt")
    parser.add_argument("--config", required=True, type=Path, help="the run's config.json")
    parser.add_argument("--out", required=True, type=Path, help="the .npz to write")
    parser.add_argument(
        "--check",
        default=False,
        action="store_true",
        help="compare float16 logits against float32 on real molecules",
    )
    parser.add_argument("--upstream", default=None, type=Path, help="the clone, for --check")
    parser.add_argument("--sdf", default=None, type=Path, help="molecules for --check")
    return parser.parse_args()


def to_archive(state: dict[str, torch.Tensor]) -> tuple[dict[str, np.ndarray], int, int]:
    arrays: dict[str, np.ndarray] = {}
    halved = 0
    kept = 0
    for key, tensor in state.items():
        value = tensor.detach().cpu().numpy()
        if value.size >= FLOAT32_BELOW:
            arrays[key] = value.astype(np.float16)
            halved += value.size
        else:
            arrays[key] = value.astype(np.float32)
            kept += value.size
    return arrays, halved, kept


def check(args: argparse.Namespace, state: dict[str, torch.Tensor], config: dict) -> None:
    """Round-trip the archive through the torch model and compare logits."""
    import sys

    sys.path.insert(0, str(args.upstream))
    from dataset_utils import mol_to_graph
    from model import ATTNSOM
    from rdkit import Chem, RDLogger
    from torch_geometric.data import Batch

    RDLogger.DisableLog("rdApp.*")
    cyps = config["cyp_names"]

    exact = ATTNSOM(**config)
    exact.load_state_dict(state)
    exact.eval()

    with np.load(args.out) as archive:
        restored = {k: torch.from_numpy(archive[k].astype(np.float32)) for k in archive.files}
    rounded = ATTNSOM(**config)
    rounded.load_state_dict(restored)
    rounded.eval()

    mols = []
    for mol in Chem.SDMolSupplier(str(args.sdf)):
        if mol is not None and mol.GetNumAtoms() > 3:
            mols.append(mol)
        if len(mols) >= 25:
            break

    worst_logit = 0.0
    worst_prob = 0.0
    flipped = 0
    rows = 0
    # A flip only matters if the two atoms were meaningfully apart to begin with; a 1e-4 shift
    # reordering a near-tie is a tie-break, not a loss of accuracy.
    flip_margins: list[float] = []
    for mol in mols:
        smiles = Chem.MolToSmiles(mol)
        graphs = []
        for cyp_idx in range(len(cyps)):
            data, _ = mol_to_graph(mol, smiles, som_indices=[])
            data.cyp_idx = torch.full((data.num_nodes,), cyp_idx, dtype=torch.long)
            graphs.append(data)
        batch = Batch.from_data_list(graphs)
        with torch.no_grad():
            a, _, _ = exact(batch)
            b, _, _ = rounded(batch)
        a = a.view(len(cyps), mol.GetNumAtoms())
        b = b.view(len(cyps), mol.GetNumAtoms())
        worst_logit = max(worst_logit, (a - b).abs().max().item())
        worst_prob = max(worst_prob, (a.sigmoid() - b.sigmoid()).abs().max().item())
        for row_a, row_b in zip(a, b, strict=True):
            rows += 1
            if row_a.argmax() != row_b.argmax():
                flipped += 1
                ranked = row_a.sigmoid().sort(descending=True).values
                flip_margins.append((ranked[0] - ranked[1]).item())

    print()
    print(f"--check on {len(mols)} molecules x {len(cyps)} isoforms:")
    print(f"  largest logit difference       : {worst_logit:.6f}")
    print(f"  largest probability difference : {worst_prob:.6f}")
    print(f"  top-1 site changed             : {flipped} of {rows}")
    if flip_margins:
        gaps = ", ".join(f"{m:.2e}" for m in flip_margins)
        print(f"  their float32 top-1/top-2 gap : {gaps}")
    widest = max(flip_margins, default=0.0)
    if worst_prob < 0.005 and widest < 0.01:
        print("  => float16 is harmless: the only reordering is of atoms that were tied.")
        print("     Ship float16 (17 MB).")
    else:
        print("  => float16 moves a real prediction. Store the large matrices float32 (34 MB).")


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text())
    state = torch.load(args.ckpt, map_location="cpu", weights_only=True)

    arrays, halved, kept = to_archive(state)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Not compressed: npz deflate buys little on weights and costs load time on every start.
    np.savez(args.out, **arrays)

    total = halved + kept
    size = args.out.stat().st_size
    print(f"{args.ckpt.name} -> {args.out}")
    print(f"  tensors {len(arrays)}, parameters {total:,}")
    print(f"  float16 {halved:,} ({halved / total:.1%}), float32 {kept:,}")
    print(f"  file {size / 1e6:.2f} MB (float32 throughout would be {total * 4 / 1e6:.2f} MB)")

    if args.check:
        if not args.upstream or not args.sdf:
            raise SystemExit("--check needs --upstream and --sdf")
        check(args, state, config)


if __name__ == "__main__":
    main()
