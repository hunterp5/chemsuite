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

"""Sketcher review: the sketcher's R/S and E/Z labels on the CIP Validation Suite.

    python scripts/review/cip_review.py COMPOUNDS_2D_SDF

COMPOUNDS_2D_SDF is ``compounds_2d.sdf`` from github.com/CIPValidationSuite/ValidationSuite
(Hanson, Mayfield, Vainio, Yerin, Musacchio and Redkin, J. Chem. Inf. Model. 2018, 58, 1755:
300 structures with the labels the authors agreed). Each record is opened in the sketcher as a
user would paste it (``convert.from_mol``: drawn coordinates and wedges kept), labelled by
``stereo.labels`` and compared with the suite's tetrahedral (R/S/r/s) and double-bond (E/Z)
labels. Records whose stereo units include axial, helical or extended (allene) units, which
the sketcher does not label, are counted separately. Results print as Markdown for
docs/sketcher/stereo-review.md.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdCIPLabeler

from chemsuite.sketch import convert, stereo

_LABEL = re.compile(r"(\d+)([RSrsEZMPmp])")


def expected(mol: Chem.Mol) -> dict[int, str]:
    """Suite labels as 0-based atom index -> label."""
    return {int(n) - 1: code for n, code in _LABEL.findall(mol.GetProp("CIP_LABELS"))}


def ours(mol: Chem.Mol) -> dict[int, str]:
    """The sketcher's labels, as 0-based atom index -> label (both ends of an E/Z bond)."""
    sketch = convert.from_mol(mol)
    order = convert.atom_order(sketch)
    index = {aid: i for i, aid in enumerate(order)}
    found = stereo.labels(sketch)
    out = {index[aid]: code for aid, code in found.atoms.items()}
    for bid, code in found.bonds.items():
        bond = sketch.bonds[bid]
        out[index[bond.begin]] = code
        out[index[bond.end]] = code
    return out


def direct(mol: Chem.Mol) -> dict[int, str]:
    """RDKit's new labeller on the suite's own molfile, without the sketcher."""
    work = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(work)
    out = {a.GetIdx(): a.GetProp("_CIPCode") for a in work.GetAtoms() if a.HasProp("_CIPCode")}
    for bond in work.GetBonds():
        if bond.HasProp("_CIPCode"):
            out[bond.GetBeginAtomIdx()] = bond.GetProp("_CIPCode")
            out[bond.GetEndAtomIdx()] = bond.GetProp("_CIPCode")
    return out


def _fmt(labels: dict[int, str] | str) -> str:
    if isinstance(labels, str):
        return labels
    return " ".join(f"{i + 1}{c}" for i, c in sorted(labels.items())) or "(none)"


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    path = Path(sys.argv[1])
    agree, differ, skipped, unreadable = [], [], [], 0
    for mol in Chem.SDMolSupplier(str(path), removeHs=False):
        if mol is None:
            unreadable += 1
            continue
        if not mol.HasProp("CIP_LABELS"):
            skipped.append((mol.GetProp("STRUCTURE_ID"), "no labels"))
            continue
        units = set(mol.GetProp("STEREOGENIC_UNITS").replace(",", " ").split())
        if units - {"TH", "CT"}:
            skipped.append((mol.GetProp("STRUCTURE_ID"), " ".join(sorted(units))))
            continue
        want = expected(mol)
        try:
            got = ours(mol)
        except convert.SketchError:
            differ.append(
                (mol.GetProp("STRUCTURE_ID"), want, "unreadable by the sketcher", direct(mol))
            )
            continue
        if got == want:
            agree.append(mol.GetProp("STRUCTURE_ID"))
        else:
            differ.append((mol.GetProp("STRUCTURE_ID"), want, got, direct(mol)))
    tested = len(agree) + len(differ)
    print(f"{tested} records with only tetrahedral and E/Z units: {len(agree)} agree, "
          f"{len(differ)} differ. {len(skipped)} skipped, {unreadable} unreadable.\n")  # fmt: skip
    same_as_rdkit = sum(got == rdkit for *_, got, rdkit in differ)
    print(f"Of the {len(differ)} that differ, {same_as_rdkit} give exactly what RDKit's labeller "
          "gives on the suite's own molfile.\n")  # fmt: skip
    if differ:
        print("| Record | Suite (atom + label) | Sketcher | RDKit on the suite molfile |")
        print("|---|---|---|---|")
        for rid, want, got, rdkit in differ:
            print(f"| {rid} | {_fmt(want)} | {_fmt(got)} | {_fmt(rdkit)} |")
    kinds: dict[str, int] = {}
    for _rid, why in skipped:
        kinds[why] = kinds.get(why, 0) + 1
    print("\nSkipped by stereo unit: " + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))


if __name__ == "__main__":
    main()
