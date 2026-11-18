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
"""Record what upstream's torch model predicts, so the numpy port can be held to it.

docs/SOM_PLAN.md, M1. Writes tests/metabolism/som_golden.json:

  * ``molecules``  - per-atom probabilities for all nine isoforms, for a spread of molecules,
    computed through the **batched** path (nine isoform graphs in one batch), which is the
    configuration the published accuracy was measured in and the only one the app will use.
  * ``layers``     - a digest of the tensor after every stage (atom encoder, each GraphCliff
    layer, SAGPool, FiLM, attention) for a few molecules: shape, the first sixteen values, and
    mean/std/min/max/norm. Enough to bisect which stage a mismatched port diverges at, without
    committing megabytes of 768-wide activations.
  * ``layers_full`` - every stage's complete tensor for one deliberately tiny molecule, so the
    port can also be checked elementwise once the digests point at a stage.
  * ``features``   - the raw 33-dim atom and 13-dim bond matrices, to catch featurizer drift
    (the Crippen and TPSA contributions come from private RDKit APIs).

The molecule list deliberately includes a symmetric molecule, because SAGPooling's top-k has
to break ties between equally scored atoms and the port must break them the same way.

Run in WSL, in the training env, after convert_weights.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

# A spread of CYP substrates plus awkward cases: symmetry (for tie-breaks), charges, P and S,
# halogens, fused rings, a single heavy atom, and a molecule with no bonds at all.
SMILES: list[str] = [
    "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",  # caffeine
    "COc1ccc2cc1OC[C@@H]1CCN(C)CC[C@H]21",  # dextromethorphan core
    "CC(=O)Nc1ccc(O)cc1",  # paracetamol
    "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",  # ibuprofen
    "Clc1ccccc1C1=NCC(=O)Nc2ccc(Cl)cc21",  # a benzodiazepine
    "CCOP(=S)(OCC)Oc1ccc(cc1)[N+](=O)[O-]",  # parathion: P, S, nitro
    "c1ccc2c(c1)[nH]c1ccccc12",  # carbazole, symmetric
    "c1ccccc1",  # benzene: fully symmetric, every atom tied
    "OC(=O)c1ccccc1O",  # salicylic acid
    "FC(F)(F)c1ccccc1N",  # CF3 aniline
    "C[S+](C)[O-]",  # charged sulfur
    "CCCCCCCC",  # flexible chain, symmetric ends
    "c1ccc(cc1)c1ccccc1",  # biphenyl, symmetric
    "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21",  # chlorpromazine
    "COc1cc2c(cc1OC)C(=O)C(CC2)Cc1ccc(OC)c(OC)c1",  # many equivalent methoxys
    "O",  # water: one heavy atom, no bonds
    "CO",  # smallest two-atom case
    "c1ccncc1",  # pyridine
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CCN(CC)C(=O)c1ccccc1",  # a simple amide
]

# Digests are kept for these; full activations only for FULL_MOLECULE, which is two heavy atoms
# and one bond, so the archive stays small while still exercising message passing.
DIGEST_MOLECULES = ("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "c1ccccc1", "CO")
FULL_MOLECULE = "CO"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--ckpt", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--decimals", default=6, type=int)
    return parser.parse_args()


def rounded(array: np.ndarray, decimals: int) -> list:
    return np.round(np.asarray(array, dtype=np.float64), decimals).tolist()


def digest(array: np.ndarray, decimals: int) -> dict:
    """Shape, the leading values, and summary statistics: enough to localize a divergence."""
    flat = np.asarray(array, dtype=np.float64).ravel()
    return {
        "shape": list(np.shape(array)),
        "head": rounded(flat[:16], decimals),
        "mean": round(float(flat.mean()), decimals),
        "std": round(float(flat.std()), decimals),
        "min": round(float(flat.min()), decimals),
        "max": round(float(flat.max()), decimals),
        "norm": round(float(np.linalg.norm(flat)), decimals),
    }


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(args.upstream))

    from dataset_utils import mol_to_graph
    from model import ATTNSOM
    from rdkit import Chem, RDLogger
    from torch_geometric.data import Batch

    RDLogger.DisableLog("rdApp.*")

    config = json.loads(args.config.read_text())
    cyps = config["cyp_names"]
    model = ATTNSOM(**config)
    model.load_state_dict(torch.load(args.ckpt, map_location="cpu", weights_only=True))
    model.eval()

    # Capture every stage by hooking the modules, rather than re-deriving the forward pass.
    captured: dict[str, np.ndarray] = {}

    def grab(name: str):
        def hook(_module, _inputs, output):
            value = output[0] if isinstance(output, tuple) else output
            if isinstance(value, torch.Tensor):
                captured[name] = value.detach().cpu().numpy()

        return hook

    handles = [model.atom_encoder.register_forward_hook(grab("atom_encoder"))]
    for index, layer in enumerate(model.encoder.layers):
        handles.append(layer.register_forward_hook(grab(f"graphcliff_{index}")))
    handles.append(model.sagpool.register_forward_hook(grab("sagpool")))
    handles.append(model.graph_film.register_forward_hook(grab("film")))
    handles.append(model.attn_head.register_forward_hook(grab("attention")))

    out: dict[str, object] = {
        "note": (
            "Generated by scripts/som/golden.py from upstream ATTNSOM (MIT) through the batched "
            "nine-isoform path. Probabilities are sigmoid(logits). docs/SOM_PLAN.md M1."
        ),
        "isoforms": cyps,
        "config": config,
        "checkpoint": args.ckpt.name,
        "rdkit": Chem.rdBase.rdkitVersion,
        "torch": torch.__version__,
        "molecules": [],
        "layers": {},
        "layers_full": {},
        "features": {},
    }

    for smiles in SMILES:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise SystemExit(f"could not parse {smiles!r}")
        canonical = Chem.MolToSmiles(mol)

        graphs = []
        for cyp_idx in range(len(cyps)):
            data, _ = mol_to_graph(mol, canonical, som_indices=[])
            data.cyp_idx = torch.full((data.num_nodes,), cyp_idx, dtype=torch.long)
            graphs.append(data)

        captured.clear()
        with torch.no_grad():
            logits, _, _ = model(Batch.from_data_list(graphs).to("cpu"))

        n_atoms = mol.GetNumAtoms()
        probabilities = torch.sigmoid(logits).view(len(cyps), n_atoms).numpy()

        out["molecules"].append(
            {
                "smiles": smiles,
                "canonical": canonical,
                "n_atoms": n_atoms,
                "symbols": [a.GetSymbol() for a in mol.GetAtoms()],
                # [isoform][atom], isoforms in the order of "isoforms" above.
                "probabilities": rounded(probabilities, args.decimals),
            }
        )

        if smiles in DIGEST_MOLECULES:
            single, _ = mol_to_graph(mol, canonical, som_indices=[])
            out["features"][smiles] = {
                "x": rounded(single.x.numpy(), args.decimals),
                "edge_index": single.edge_index.numpy().tolist(),
                "edge_attr": rounded(single.edge_attr.numpy(), args.decimals),
            }
            out["layers"][smiles] = {
                name: digest(value, args.decimals) for name, value in sorted(captured.items())
            }
            if smiles == FULL_MOLECULE:
                out["layers_full"][smiles] = {
                    name: rounded(value, args.decimals) for name, value in sorted(captured.items())
                }

    for handle in handles:
        handle.remove()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1))
    size = args.out.stat().st_size
    print(f"{len(SMILES)} molecules; digests for {len(out['layers'])}, full for {FULL_MOLECULE}")
    print(f"wrote {args.out} ({size / 1e6:.2f} MB)")
    print("stages captured:", ", ".join(sorted(out["layers"][DIGEST_MOLECULES[0]])))


if __name__ == "__main__":
    main()
