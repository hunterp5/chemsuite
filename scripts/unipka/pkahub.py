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

"""pKaHub's measured pKa values as chemsuite's training set (docs/UNIPKA_PLAN.md, U1).

pKaHub (Sipos-Szabó, Bajusz, Balogh and Keserű, J. Chem. Inf. Model. 2026,
doi:10.1021/acs.jcim.6c00107; https://github.com/keserulab/pkahub, MIT) compiles published
pKa values and assigns each to a macrostate transition with Epik, listing the microstates on
both sides. Only the sources in :data:`SHIPPED` are used (decided 2026-11-12): not the IUPAC
digitised set (CC BY-NC), DataWarrior (unknown provenance) or the AttenGpKa set (which may
repeat them), and never the test sets.

    python scripts/unipka/pkahub.py <pka_all_datasets.tsv> <out.tsv> --tests <folder>

The output is Uni-pKa's format: "acid microstates>>base microstates" and the pKa.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[2]
SHIPPED = frozenset(["Baltruschat ChEMBL", "Hunt", "Organic Oxygen Acids and Nitrogen Bases",
                     "Settimo", "Manchester", "Caine", "Jensen", "OCHEM",
                     "QSARToolbox"])  # fmt: skip
TESTS = frozenset(["Novartis", "AvLiLuMoVe", "SAMPL6", "SAMPL7", "SAMPL8", "euroSAMPL1"])
ELEMENTS = frozenset(["C", "N", "O", "S", "H", "Cl", "F", "Br", "I", "Si", "P", "B"])
MAX_HEAVY = 60
MAX_ASSIGNMENT_ERROR = 2.0  # |measured − Epik| beyond this: the transition may be misassigned
TEMPERATURE = (20.0, 30.0)  # °C, when given
SPREAD = 1.0  # repeated values of one transition spread wider than this are dropped


def _skeleton(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(rdMolStandardize.ChargeParent(mol)).split("-")[0]


def held_out(tests: Path, rows: list[dict[str, str]]) -> set[str]:
    """Skeletons of pKaHub's test sets, the Uni-pKa test files and the neutron set's ligands."""
    out = {_skeleton(r["smiles"]) for r in rows if r["data_source"] in TESTS}
    for path in sorted(tests.glob("*.tsv")):
        for line in path.read_text("utf-8").splitlines()[1:]:
            for smiles in re.split(r"[\t,>]+", line):
                if smiles and Chem.MolFromSmiles(smiles) is not None:
                    out.add(_skeleton(smiles))
    for cif in sorted((ROOT / "samples" / "prepare" / "ccd").glob("*.cif")):
        found = re.search(r"InChIKey\s+InChI\s+\S+\s+\"?([A-Z]{14})-", cif.read_text("utf-8"))
        if found:
            out.add(found.group(1))
    out.discard(None)
    return out


def _usable(smiles: str) -> bool:
    """Every microstate parses, has only the model's elements, at most :data:`MAX_HEAVY` heavy
    atoms, and no atom charged beyond ±1 (the charge dictionary's range)."""
    for s in smiles.split(","):
        mol = Chem.MolFromSmiles(s)
        if mol is None or mol.GetNumHeavyAtoms() > MAX_HEAVY:
            return False
        if any(a.GetSymbol() not in ELEMENTS or abs(a.GetFormalCharge()) > 1
               for a in Chem.AddHs(mol).GetAtoms()):  # fmt: skip
            return False
    return True


def build(src: Path, dst: Path, tests: Path) -> None:
    RDLogger.DisableLog("rdApp.*")
    rows = list(csv.DictReader(src.open(encoding="utf-8"), delimiter="\t"))
    hold = held_out(tests, rows)
    print(f"held-out skeletons: {len(hold)}")
    counts: Counter[str] = Counter(rows=len(rows))
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    skeletons: dict[str, str | None] = {}
    for r in rows:
        if r["data_source"] not in SHIPPED:
            counts["source not shipped"] += 1
            continue
        if r["temperature"] and not TEMPERATURE[0] <= float(r["temperature"]) <= TEMPERATURE[1]:
            counts["temperature"] += 1
            continue
        if r["assignment_error"] and float(r["assignment_error"]) > MAX_ASSIGNMENT_ERROR:
            counts["assignment error"] += 1
            continue
        pre, post = r["microspecies_pre"].strip('"'), r["microspecies_post"].strip('"')
        if not pre or not post:
            counts["no microspecies"] += 1
            continue
        if r["molid"] not in skeletons:
            skeletons[r["molid"]] = _skeleton(r["smiles"])
        if skeletons[r["molid"]] in hold:
            counts["held out"] += 1
            continue
        key = (pre, post)
        values[key].append(float(r["pka_value"]))
        sources[key].add(r["data_source"])
    out = []
    for (pre, post), found in values.items():
        if max(found) - min(found) > SPREAD:
            counts["values disagree"] += 1
            continue
        if not _usable(pre) or not _usable(post):
            counts["structure not usable"] += 1
            continue
        out.append(
            (f"{pre}>>{post}", float(np.median(found)), ";".join(sorted(sources[pre, post])))
        )
    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", encoding="utf-8", newline="\n") as f:
        f.write("\tSMILES\tTARGET\tref.\n")
        for k, (pair, value, ref) in enumerate(out):
            f.write(f"{k}\t{pair}\t{value:.3f}\t{ref}\n")
    for name, n in counts.items():
        print(f"{name}: {n}")
    targets = np.array([v for _p, v, _r in out])
    print(f"transitions written: {len(out)}, mean pKa {targets.mean():.4f}")


if __name__ == "__main__":
    if len(sys.argv) < 5 or "--tests" not in sys.argv:
        sys.exit(__doc__)
    build(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[sys.argv.index("--tests") + 1]))
