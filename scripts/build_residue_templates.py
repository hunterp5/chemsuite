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
"""Build ``chemsuite/protein/data/residues.json`` from the wwPDB Chemical Component Dictionary.

    python scripts/build_residue_templates.py

The file is checked in, so users never need this. For each residue below it keeps the bonds
between atom names (with their order and aromaticity) and the component's SMILES; the protein
package bonds residues from it instead of guessing from distances. The CCD is in the public
domain (wwPDB, CC0 1.0).
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "chemsuite" / "protein" / "data" / "residues.json"
URL = "https://files.rcsb.org/ligands/download/{}.cif"

AMINO = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "SEC",
    "PYL",
    "MSE",
    "SEP",
    "TPO",
    "PTR",
    "HYP",
    "MLY",
    "M3L",
    "CSO",
    "CSD",
    "OCS",
    "CME",
    "KCX",
    "LLP",
    "ALY",
    "CGU",
    "FME",
]
CAPS = ["ACE", "NH2", "NME"]
NUCLEIC = ["A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI"]
#: cofactors and ligands common enough in drug-discovery entries to bond by template
COFACTORS = [
    "HEM",
    "HEC",
    "ATP",
    "ADP",
    "AMP",
    "ANP",
    "GTP",
    "GDP",
    "GNP",
    "NAD",
    "NAP",
    "NDP",
    "FAD",
    "FMN",
    "SAM",
    "SAH",
    "COA",
    "PLP",
    "ACO",
    "TPP",
]
ORDERS = {"SING": 1, "DOUB": 2, "TRIP": 3}


def template(code: str) -> dict:
    with urllib.request.urlopen(URL.format(code), timeout=30) as response:
        block = gemmi.cif.read_string(response.read().decode()).sole_block()
    bonds = [
        [a, b, ORDERS.get(order.upper(), 1), int(aromatic.upper() == "Y")]
        for a, b, order, aromatic in block.find(
            "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"]
        )
    ]
    canonical = {
        gemmi.cif.as_string(program): gemmi.cif.as_string(text)
        for kind, program, text in block.find(
            "_pdbx_chem_comp_descriptor.", ["type", "program", "descriptor"]
        )
        if kind == "SMILES_CANONICAL"
    }
    smiles = canonical.get("OpenEye OEToolkits") or canonical.get("CACTVS", "")
    names = [[gemmi.cif.as_string(a), gemmi.cif.as_string(b), o, r] for a, b, o, r in bonds]
    return {"bonds": names, "smiles": smiles}


def main() -> int:
    codes = AMINO + CAPS + NUCLEIC + COFACTORS
    templates = {}
    for code in codes:
        templates[code] = template(code)
        print(code, len(templates[code]["bonds"]), "bonds")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(templates, separators=(",", ":"), sort_keys=True) + "\n")
    print(f"wrote {len(templates)} templates to {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
