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

"""Sketcher review: chemsuite's abbreviation table against Ketcher's and RDKit's.

    python scripts/review/abbreviation_review.py KETCHER_FG_SDF

KETCHER_FG_SDF is Ketcher's ``packages/ketcher-react/src/templates/fg.sdf`` (commit 2b031cc2).
RDKit's list is ``rdkit.Chem.rdAbbreviations.GetDefaultAbbreviations()``. Groups are compared
by the canonical SMILES of the group bonded to ``[*]``. Results print as Markdown and are
discussed in docs/sketcher/abbreviation-review.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdAbbreviations

from chemsuite.sketch.abbrev import TABLE


def _canon(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def ketcher_records(path: Path) -> dict[str, str]:
    """name -> canonical group SMILES on [*], or a reason the record is unusable."""
    out: dict[str, str] = {}
    for block in path.read_text(encoding="utf-8").split("$$$$\n"):
        block = block.lstrip("\n")  # Ketcher puts a blank line after each $$$$
        if "M  END" not in block:
            continue
        name = block.splitlines()[0].strip()
        molblock = block[: block.index("M  END") + len("M  END")]
        mol = Chem.MolFromMolBlock(molblock)
        if mol is None:
            out[name] = "unreadable"
            continue
        sap = [line for line in molblock.splitlines() if line.startswith("M  SAP")]
        if not sap:
            out[name] = "no attachment point"
            continue
        rw = Chem.RWMol(mol)
        star = rw.AddAtom(Chem.Atom(0))
        rw.AddBond(int(sap[0].split()[4]) - 1, star, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except (ValueError, RuntimeError) as err:
            out[name] = f"invalid ({type(err).__name__})"
            continue
        out[name] = Chem.MolToSmiles(rw)
    return out


def rdkit_records() -> dict[str, str]:
    return {a.label: Chem.MolToSmiles(a.mol) for a in rdAbbreviations.GetDefaultAbbreviations()}


def compare(ours: dict[str, str], theirs: dict[str, str], source: str) -> list[str]:
    rows = [f"| Label | chemsuite | {source} | Verdict |", "|---|---|---|---|"]
    tally = {"agree": 0, "differ": 0, "unusable": 0, "absent": 0}
    lower = {k.lower(): v for k, v in theirs.items()}
    for entry in TABLE:
        mine = ours[entry.label]
        found = next(
            (lower[k.lower()] for k in (entry.label, *entry.aliases) if k.lower() in lower), None
        )
        if found is None:
            tally["absent"] += 1
            continue
        if not found.startswith("*") and "*" not in found:
            verdict, tally["unusable"] = f"{source} record {found}", tally["unusable"] + 1
        elif _canon(found) == mine:
            tally["agree"] += 1
            continue
        else:
            verdict, tally["differ"] = "**differ**", tally["differ"] + 1
        rows.append(f"| {entry.label} | `{mine}` | `{found}` | {verdict} |")
    head = ", ".join(f"{k} {v}" for k, v in tally.items())
    return [f"**{source}: {head}.** Rows below are the ones that do not agree.", "", *rows]


def main(argv: list[str]) -> int:
    RDLogger.DisableLog("rdApp.*")
    ours = {e.label: _canon(e.smiles) for e in TABLE}
    ketcher = ketcher_records(Path(argv[1]))
    out = sys.stdout
    out.write(f"## chemsuite table: {len(TABLE)} abbreviations\n\n")
    out.write("## Against Ketcher fg.sdf\n\n")
    out.write("\n".join(compare(ours, ketcher, "Ketcher")) + "\n\n")
    bad = {k: v for k, v in ketcher.items() if "*" not in v}
    out.write(f"Ketcher records RDKit cannot use ({len(bad)}): ")
    out.write(", ".join(f"{k} ({v})" for k, v in sorted(bad.items())) + "\n\n")
    out.write("## Against RDKit's default abbreviations\n\n")
    out.write("\n".join(compare(ours, rdkit_records(), "RDKit")) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
