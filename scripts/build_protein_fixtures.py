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

    python scripts/build_protein_fixtures.py            # all of them
    python scripts/build_protein_fixtures.py 5gnk       # only those whose name starts so

PDB entries are in the public domain (wwPDB, CC0 1.0). Each covers something the protein
package must handle:

- 1CRN, crambin: small, one chain, helices and a sheet (PDB and mmCIF)
- 3HTB, T4 lysozyme L99A: a ligand (JZ4), waters, additives and alternate locations
- 1ATP, protein kinase A: ATP and two Mn ions, phosphorylated residues, a peptide chain
- 1PPB, thrombin: insertion codes and a covalently bound inhibitor (PPACK, 0G6)
- 1D3Z, ubiquitin (NMR): cut to its first two models
- 5GNK, EGFR T790M with LXX-6-34 (80U): cut to the residues within 12 Å of 80U. A covalent
  warhead on Cys797, a halogen bond to Leu788 O, and a drug-like ligand that is not among the
  bundled residue templates; the entry's own ``_chem_comp_atom``/``_chem_comp_bond`` rows are
  kept for the components left, as RCSB files carry them
- 0G6 (PPACK) as the Chemical Component Dictionary serves it, with its leaving atom flagged:
  a stand-in for the online dictionary
"""

from __future__ import annotations

import gzip
import sys
import urllib.request
from pathlib import Path

import gemmi

from chemsuite.core.download import ssl_context

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "data" / "protein"
URL = "https://files.rcsb.org/download/{}"
CCD_URL = "https://files.rcsb.org/ligands/download/{}.cif"
#: components fetched from the CCD itself
COMPONENTS = ["0G6"]
ENTRIES = ["1CRN.pdb", "1CRN.cif", "3HTB.cif", "1ATP.cif", "1PPB.cif"]
#: (entry, ligand, radius in Å): cut to the residues near the ligand
POCKETS = [("5GNK", "80U", 12.0)]
#: the dictionary categories copied for the components left in a cut entry
COMPONENT_TABLES = ("_chem_comp_atom.", "_chem_comp_bond.")


def fetch(name: str, url: str = URL) -> bytes:
    with urllib.request.urlopen(url.format(name), timeout=60, context=ssl_context()) as response:
        return response.read()


def pocket(text: str, ligand: str, radius: float) -> str:
    """The entry cut to the residues within ``radius`` of ``ligand``, as mmCIF: the links
    between the residues left, and their rows of the component dictionary tables."""
    block = gemmi.cif.read_string(text).sole_block()
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    model = st[0]
    search = gemmi.NeighborSearch(model, st.cell, 5).populate()
    keep = set()
    for atom in (a for chain in model for res in chain if res.name == ligand for a in res):
        for mark in search.find_atoms(atom.pos, "\0", radius=radius):
            cra = mark.to_cra(model)
            keep.add((cra.chain.name, str(cra.residue.seqid), cra.residue.name))
    for chain in model:
        for i in reversed(range(len(chain))):
            if (chain.name, str(chain[i].seqid), chain[i].name) not in keep:
                del chain[i]
    st.remove_empty_chains()

    def kept(p: gemmi.AtomAddress) -> bool:
        return (p.chain_name, str(p.res_id.seqid), p.res_id.name) in keep

    st.connections = gemmi.ConnectionList(
        [c for c in st.connections if kept(c.partner1) and kept(c.partner2)]
    )
    out = st.make_mmcif_document()
    left = {res.name for chain in model for res in chain}
    for category in COMPONENT_TABLES:
        table = block.find_mmcif_category(category)
        tags = [tag.removeprefix(category) for tag in table.tags]
        loop = out.sole_block().init_mmcif_loop(category, tags)
        for row in table:
            if row[0] in left:
                loop.add_row([row.str(i) if row.has(i) else "?" for i in range(len(tags))])
    return out.as_string()


def main(only: list[str]) -> int:
    def wanted(name: str) -> bool:
        return not only or any(name.lower().startswith(o.lower()) for o in only)

    OUT.mkdir(parents=True, exist_ok=True)
    for name in filter(wanted, ENTRIES):
        (OUT / f"{name.lower()}.gz").write_bytes(gzip.compress(fetch(name), mtime=0))
    if wanted("1d3z"):
        text = fetch("1D3Z.cif").decode()
        nmr = gemmi.read_structure_string(text, format=gemmi.CoorFormat.Mmcif)
        while len(nmr) > 2:
            del nmr[len(nmr) - 1]
        text = nmr.make_mmcif_document().as_string()
        (OUT / "1d3z_2models.cif.gz").write_bytes(gzip.compress(text.encode(), mtime=0))
    for entry, ligand, radius in POCKETS:
        if wanted(entry):
            text = pocket(fetch(f"{entry}.cif").decode(), ligand, radius)
            path = OUT / f"{entry.lower()}_pocket.cif.gz"
            path.write_bytes(gzip.compress(text.encode(), mtime=0))
    for name in COMPONENTS:
        if wanted(f"{name}_ccd"):
            data = fetch(name, CCD_URL)
            (OUT / f"{name.lower()}_ccd.cif.gz").write_bytes(gzip.compress(data, mtime=0))
    for path in sorted(OUT.iterdir()):
        print(f"{path.name}: {path.stat().st_size // 1024} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
