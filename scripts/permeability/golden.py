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

"""Write ``tests/permeability/permeability_golden.json``: the shipped models' numbers, pinned.

Run after ``train.py --ship`` replaces an archive, and read the diff: a change in the
descriptor catalog, the feature order or the prepare step shows up here as changed numbers,
and a change in a model shows up as all of that endpoint's numbers changing.

    .venv/Scripts/python.exe scripts/permeability/golden.py

The molecules are marketed drugs whose permeability and efflux behavior is textbook, plus a
few small cases that exercise charges, halogens and rings.
"""

from __future__ import annotations

import json
from pathlib import Path

from chemsuite.permeability import ENDPOINTS
from chemsuite.permeability import models as shipped
from chemsuite.permeability.predict import predict_permeability

GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "permeability" / "permeability_golden.json"

MOLECULES = {
    # the classic Caco-2 pair: lipophilic and permeable, polar and not
    "propranolol": "CC(C)NCC(O)COc1cccc2ccccc12",
    "atenolol": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
    # a P-glycoprotein substrate, which should show a high efflux ratio
    "verapamil": "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
    # everyday small drugs
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "diclofenac": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    "warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "metformin": "CN(C)C(=N)N=C(N)N",
    # small cases that are not drugs at all
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
    "glycine": "NCC(=O)O",
    "benzamidine": "NC(=N)c1ccccc1",
    # a bare trisaccharide: nothing like the training chemistry, and a probe that the
    # applicability domain notices. It is deliberately not labelled as a drug — an earlier
    # version of this file called it digoxin, which it is not: digoxin is C41H64O14 with a
    # steroid core, and this is only the sugar chain, C18H32O13.
    "sugar_trimer": "CC1OC(OC2C(O)C(OC3OC(C)C(O)C(O)C3O)OC(C)C2O)C(O)C(O)C1O",
}

#: decimals kept, which is far tighter than the assay and the point: this pins the arithmetic
PLACES = 7


def main() -> None:
    missing = [e.id for e in ENDPOINTS if not shipped.installed(e)]
    if missing:
        raise SystemExit(
            f"not every model is built: {', '.join(missing)}.\n"
            "Run scripts/permeability/train.py --ship <arm> first."
        )

    smiles = list(MOLECULES.values())
    results = predict_permeability(smiles, ENDPOINTS)
    endpoints: dict[str, dict] = {}
    for endpoint in ENDPOINTS:
        said = results.per_endpoint[endpoint.id].as_prediction()
        entry: dict[str, object] = {
            "task": endpoint.task,
            "predicted": [round(float(value), PLACES) for value in said.predicted],
            "distance": [round(float(value), PLACES) for value in said.distance],
            "inside": [bool(flag) for flag in said.inside],
        }
        if endpoint.task == "regression":
            entry["low"] = [round(float(value), PLACES) for value in said.low]
            entry["high"] = [round(float(value), PLACES) for value in said.high]
        else:
            entry["scores"] = [
                [round(float(value), PLACES) for value in row] for row in said.scores
            ]
        endpoints[endpoint.id] = entry

    GOLDEN.write_text(
        json.dumps({"names": MOLECULES, "endpoints": endpoints}, indent=1) + "\n",
        encoding="utf-8",
    )
    print(
        f"wrote {GOLDEN}: {len(MOLECULES)} molecules x {len(ENDPOINTS)} endpoints, "
        f"{GOLDEN.stat().st_size / 1024:.0f} KB"
    )


if __name__ == "__main__":
    main()
