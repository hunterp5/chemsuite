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
"""Write ``chemsuite/macromolecule/look/data/amber_charges.json`` from PDB2PQR's ``AMBER.DAT``.

The charges are Amber's RESP charges for the standard residues (Cornell et al., J. Am. Chem.
Soc. 1995, 117, 5179), which ff99SB and ff14SB keep; PDB2PQR distributes them under BSD-3
(its ``LICENSE.md`` is copied beside the table). Only the numbers are used: the file is read
as text and pdb2pqr is never imported.

    python scripts/build_amber_charges.py path/to/AMBER.DAT path/to/LICENSE.md

Download both from https://github.com/Electrostatics/pdb2pqr (``pdb2pqr/dat/AMBER.DAT`` and
``LICENSE.md``) into a folder of their own first.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "chemsuite" / "macromolecule" / "look" / "data"
SOURCE = "https://github.com/Electrostatics/pdb2pqr/blob/master/pdb2pqr/dat/AMBER.DAT"


def read(path: Path) -> dict[str, dict[str, float]]:
    """``{residue: {atom: charge}}`` from AMBER.DAT's tab-separated lines."""
    out: dict[str, dict[str, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 3 or line.startswith(("#", "---")):
            continue
        residue, atom, charge = fields[0], fields[1], float(fields[2])
        out.setdefault(residue, {})[atom] = round(charge, 6)
    return out


def main(dat: str, licence: str) -> None:
    table = read(Path(dat))
    OUT.mkdir(parents=True, exist_ok=True)
    document = {
        "about": "Amber ff94/ff99 RESP partial charges per residue and atom (Cornell et al., "
        "J. Am. Chem. Soc. 1995, 117, 5179; kept by ff99SB and ff14SB), from PDB2PQR's "
        "AMBER.DAT (BSD-3-Clause, see LICENSE-pdb2pqr.md beside this file). Residue names "
        "follow Amber: HID/HIE/HIP, ASH, GLH, LYN, CYM/CYX, N- and C-terminal forms prefixed "
        "N and C, caps ACE/NME/NHE.",
        "source": SOURCE,
        "charges": dict(sorted(table.items())),
    }
    text = json.dumps(document, indent=1, sort_keys=False) + "\n"
    (OUT / "amber_charges.json").write_text(text, encoding="utf-8", newline="\n")
    shutil.copyfile(licence, OUT / "LICENSE-pdb2pqr.md")
    print(f"{len(table)} residues written to {OUT / 'amber_charges.json'}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
