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
"""Download the protein test files into ``tests/data/protein`` (checked in).

    python scripts/build_protein_fixtures.py

PDB entries are in the public domain (wwPDB, CC0 1.0). Each covers something the protein
package must handle:

- 1CRN, crambin: small, one chain, helices and a sheet (PDB and mmCIF)
- 3HTB, T4 lysozyme L99A: a ligand (JZ4), waters, additives and alternate locations
- 1ATP, protein kinase A: ATP and two Mn ions, phosphorylated residues, a peptide chain
- 1PPB, thrombin: insertion codes and a covalently bound inhibitor (PPACK, 0G6)
- 1D3Z, ubiquitin (NMR): cut to its first two models
"""

from __future__ import annotations

import gzip
import sys
import urllib.request
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "data" / "protein"
URL = "https://files.rcsb.org/download/{}"
ENTRIES = ["1CRN.pdb", "1CRN.cif", "3HTB.cif", "1ATP.cif", "1PPB.cif"]


def fetch(name: str) -> bytes:
    with urllib.request.urlopen(URL.format(name), timeout=60) as response:
        return response.read()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ENTRIES:
        (OUT / f"{name.lower()}.gz").write_bytes(gzip.compress(fetch(name), mtime=0))
    nmr = gemmi.read_structure_string(fetch("1D3Z.cif").decode(), format=gemmi.CoorFormat.Mmcif)
    while len(nmr) > 2:
        del nmr[len(nmr) - 1]
    text = nmr.make_mmcif_document().as_string()
    (OUT / "1d3z_2models.cif.gz").write_bytes(gzip.compress(text.encode(), mtime=0))
    for path in sorted(OUT.iterdir()):
        print(f"{path.name}: {path.stat().st_size // 1024} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
