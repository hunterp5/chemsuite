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
"""Write ``tests/metabolism/som_golden.json``: the shipped model's scores, pinned.

Run after ``train_gbm.py --export`` replaces the archive, and read the diff: a change in
the features or the trees shows up here as changed numbers, and a change in the model
shows up as all of them changing. The molecules are marketed drugs with known CYP routes
plus a few small cases that exercise rings, halogens and charges.

    python scripts/som/golden.py
"""

from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem

from chemsuite.metabolism import ISOFORMS
from chemsuite.metabolism.sites import predict_sites

GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "metabolism" / "som_golden.json"

MOLECULES = {
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "testosterone": "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O",
    "diclofenac": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    "propranolol": "CC(C)NCC(O)COc1cccc2ccccc12",
    "dextromethorphan": "COc1ccc2c(c1)C13CCCCC1N(C)CCC3C2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "omeprazole": "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
    "warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "ethanol": "CCO",
    "toluene": "Cc1ccccc1",
    "nicotine": "CN1CCCC1c1cccnc1",
    "aniline": "Nc1ccccc1",
    "benzamidine": "NC(=N)c1ccccc1",
}


def main() -> None:
    molecules = {}
    for name, smiles in MOLECULES.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"{name}: {smiles} does not parse")
        ranked = predict_sites(mol)
        molecules[smiles] = {
            isoform: [round(site.score, 7) for site in sorted(sites, key=lambda site: site.atom)]
            for isoform, sites in ranked.items()
        }

    GOLDEN.write_text(
        json.dumps(
            {"isoforms": list(ISOFORMS), "names": MOLECULES, "molecules": molecules},
            indent=1,
        )
        + "\n",
        encoding="utf-8",
    )
    atoms = sum(len(next(iter(by_isoform.values()))) for by_isoform in molecules.values())
    print(
        f"wrote {GOLDEN}: {len(molecules)} molecules, {atoms} atoms, "
        f"{GOLDEN.stat().st_size / 1024:.0f} KB"
    )


if __name__ == "__main__":
    main()
