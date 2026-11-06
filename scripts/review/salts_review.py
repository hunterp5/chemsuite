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

"""Sketcher review: the salts-and-solvents table against PubChem and Ketcher.

    python scripts/review/salts_review.py [KETCHER_SALTS_SDF]

Each entry of ``chemsuite.sketch.salts.TABLE`` is looked up on PubChem by its ``pubchem``
name (PUG REST, network needed) and compared by InChIKey: the full key, or the connectivity
block where PubChem's record has no stereo. With KETCHER_SALTS_SDF (Ketcher's
``packages/ketcher-react/src/templates/salts-and-solvents.sdf``, commit 2b031cc2) each Ketcher
record is also matched to an entry by InChIKey, and records with no matching entry are listed
with their structure. Results print as Markdown for docs/sketcher/salts-review.md.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

from rdkit import Chem, RDLogger

from chemsuite.sketch.salts import TABLE

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES,InChIKey/JSON"


class LookupFailed(RuntimeError):
    pass


def pubchem(name: str) -> tuple[str, str] | None:
    """(InChIKey, SMILES) of PubChem's first record for *name*; None when it has none.

    Fetched with curl, which uses the system's certificate store (Python's own may lack it).
    """
    url = PUG.format(urllib.parse.quote(name))
    for _attempt in range(3):
        done = subprocess.run(
            ["curl", "-s", "-m", "20", "-w", "\n%{http_code}", url],
            capture_output=True,
            text=True,
            check=False,
        )
        body, _, status = done.stdout.rpartition("\n")
        if status == "404":
            return None
        if status == "200":
            props = json.loads(body)["PropertyTable"]["Properties"][0]
            return props["InChIKey"], props.get("SMILES", "")
        time.sleep(1.0)
    raise LookupFailed(name)


def key(smiles: str) -> str:
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    print("| Entry | chemsuite SMILES | PubChem name | Verdict |")
    print("|---|---|---|---|")
    tally = {"full": 0, "connectivity": 0, "differ": 0, "missing": 0}
    for entry in TABLE:
        ours = key(entry.smiles)
        try:
            found = pubchem(entry.pubchem)
        except LookupFailed:
            sys.exit(f"PubChem could not be reached for {entry.pubchem!r}; nothing was compared.")
        time.sleep(0.25)  # PubChem asks for at most 5 requests a second
        if found is None:
            verdict, kind = "not found on PubChem", "missing"
        elif found[0] == ours:
            verdict, kind = "agrees (InChIKey)", "full"
        elif found[0].split("-")[0] == ours.split("-")[0] and found[0].split("-")[1].startswith(
            "UHFFFAOYS"
        ):
            verdict, kind = "agrees (connectivity; PubChem record has no stereo)", "connectivity"
        else:
            verdict, kind = f"**differs**: PubChem `{found[1]}`", "differ"
        tally[kind] += 1
        print(f"| {entry.name} | `{entry.smiles}` | {entry.pubchem} | {verdict} |")
    print(
        f"\n{len(TABLE)} entries: {tally['full']} agree in full, {tally['connectivity']} by "
        f"connectivity, {tally['differ']} differ, {tally['missing']} not found."
    )
    if len(sys.argv) > 1:
        keys = {key(e.smiles): e.name for e in TABLE}
        print("\n| Ketcher record | Ketcher structure | chemsuite entry |")
        print("|---|---|---|")
        for mol in Chem.SDMolSupplier(str(Path(sys.argv[1]))):
            name = mol.GetProp("name").replace("\xa0", " ")
            match = keys.get(Chem.MolToInchiKey(mol), "**none**")
            print(f"| {name} | `{Chem.MolToSmiles(Chem.RemoveHs(mol))}` | {match} |")


if __name__ == "__main__":
    main()
