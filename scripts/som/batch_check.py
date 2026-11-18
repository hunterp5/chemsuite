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
"""Measure whether a site-of-metabolism prediction depends on how the batch is built.

docs/SOM_PLAN.md, M1. Upstream pools the FiLM graph vector two ways (model.py):

    if batch.max() == 0 and batch.size(0) > 1:   # exactly one graph in the batch
        g = gap(h, batch)                        # plain mean over all atoms
    else:
        g = gap(*sagpool(h, ...))                # mean over the kept 80%

Training and evaluation always batch, so the published numbers come from the SAGPool branch,
and upstream's own ``build_batch`` batches nine isoform copies of one molecule, which also
takes it. The mean branch fires only when a single graph is passed on its own.

This script quantifies what that costs, so the port can decide knowingly:

  A. nine-isoform batch (upstream inference) vs one graph at a time -> does the answer change?
  B. nine-isoform batch vs a large mixed batch -> is the SAGPool branch itself batch-independent?

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
    parser.add_argument("--ckpt", required=True, type=Path, help="one fold*_best.pt")
    parser.add_argument("--config", required=True, type=Path, help="the run's config.json")
    parser.add_argument("--sdf", required=True, type=Path, help="an isoform SDF to sample from")
    parser.add_argument("--n-mols", default=12, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # The clone's modules import their siblings by bare name.
    sys.path.insert(0, str(args.upstream))
    from dataset_utils import mol_to_graph
    from model import ATTNSOM
    from rdkit import Chem, RDLogger
    from torch_geometric.data import Batch

    RDLogger.DisableLog("rdApp.*")

    config = json.loads(args.config.read_text())
    cyps = config["cyp_names"]
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = ATTNSOM(**config).to(device)
    model.load_state_dict(torch.load(args.ckpt, map_location=device))
    model.eval()

    supplier = Chem.SDMolSupplier(str(args.sdf))
    mols = []
    for mol in supplier:
        if mol is not None and mol.GetNumAtoms() > 3:
            mols.append(mol)
        if len(mols) >= args.n_mols:
            break
    print(f"{len(mols)} molecules from {args.sdf.name}, {len(cyps)} isoforms, device {device}")

    def graphs_for(mol: Chem.Mol) -> list:
        smiles = Chem.MolToSmiles(mol)
        out = []
        for cyp_idx in range(len(cyps)):
            data, _ = mol_to_graph(mol, smiles, som_indices=[])
            data.cyp_idx = torch.full((data.num_nodes,), cyp_idx, dtype=torch.long)
            out.append(data)
        return out

    @torch.no_grad()
    def logits_of(graph_list: list) -> torch.Tensor:
        batch = Batch.from_data_list(graph_list).to(device)
        logits, _, _ = model(batch)
        return logits.cpu()

    worst_gap = 0.0
    flipped_top1 = 0
    compared = 0

    for mol in mols:
        per_isoform = graphs_for(mol)

        # A: the whole nine-isoform batch, as upstream's inference does.
        batched = logits_of(per_isoform)
        n_atoms = per_isoform[0].num_nodes
        batched = batched.view(len(cyps), n_atoms)

        # A': each isoform's graph alone, which trips the mean-pool branch.
        alone = torch.stack([logits_of([g]).view(n_atoms) for g in per_isoform])

        gap = (batched - alone).abs().max().item()
        worst_gap = max(worst_gap, gap)
        for row_b, row_a in zip(batched, alone, strict=True):
            compared += 1
            if int(row_b.argmax()) != int(row_a.argmax()):
                flipped_top1 += 1

    print()
    print("A. nine-isoform batch vs one graph alone (SAGPool branch vs mean branch)")
    print(f"   largest absolute logit difference : {worst_gap:.4f}")
    print(f"   top-1 site changed               : {flipped_top1} of {compared} molecule-isoforms")

    # B: is the SAGPool branch itself independent of how much else is in the batch?
    every_graph: list = []
    for mol in mols:
        every_graph.extend(graphs_for(mol))
    big = logits_of(every_graph)

    offset = 0
    worst_b = 0.0
    for mol in mols:
        per_isoform = graphs_for(mol)
        n_atoms = per_isoform[0].num_nodes
        small = logits_of(per_isoform)
        span = len(cyps) * n_atoms
        worst_b = max(worst_b, (big[offset : offset + span] - small).abs().max().item())
        offset += span

    print()
    print("B. nine-isoform batch vs one big mixed batch (both on the SAGPool branch)")
    print(f"   largest absolute logit difference : {worst_b:.6f}")
    print(f"   graphs in the big batch           : {len(every_graph)}")

    print()
    verdict_a = "IDENTICAL" if worst_gap < 1e-4 else "DIFFERENT"
    verdict_b = "IDENTICAL" if worst_b < 1e-4 else "DIFFERENT"
    print(f"verdict: A {verdict_a}, B {verdict_b}")
    if verdict_b == "IDENTICAL" and verdict_a == "DIFFERENT":
        print("=> SAGPool is deterministic given >=2 graphs; the port must always batch the")
        print("   isoforms together and must implement SAGPooling. The one-graph mean branch")
        print("   is a separate, unreachable-by-design path.")
    print(f"np {np.__version__}, torch {torch.__version__}")


if __name__ == "__main__":
    main()
