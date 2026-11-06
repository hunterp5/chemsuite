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

"""Sketcher review: Ketcher's valence model compared with RDKit's, case by case.

``ketcher_hydrogens`` is ``chemsuite.sketch.ctfile.default_hydrogens``, a port of
``Atom.calcValence`` from Ketcher (github.com/epam/ketcher,
packages/ketcher-core/src/domain/entities/atom.ts, commit 2b031cc2, Apache-2.0) in its
default valence mode. For every element in
``ELEMENTS``, charge -2..+2, radical electrons 0..2 and bond-order sum 0..6 it asks both
models for the implicit hydrogen count and prints the disagreements as Markdown:

    python scripts/review/valence_review.py > /tmp/valence-cases.md

The findings are discussed in docs/sketcher/valence-review.md.
"""

from __future__ import annotations

import sys
from collections import defaultdict

from rdkit import Chem, RDLogger

from chemsuite.sketch.ctfile import ELEMENTS, INVALID
from chemsuite.sketch.ctfile import default_hydrogens as ketcher_hydrogens


def rdkit_hydrogens(label: str, charge: int, conn: int, rad: int) -> int | None:
    """RDKit's implicit H count for the same atom (bonded to *conn* dummy atoms)."""
    mol = Chem.RWMol()
    centre = Chem.Atom(label)
    centre.SetFormalCharge(charge)
    centre.SetNumRadicalElectrons(rad)
    index = mol.AddAtom(centre)
    for _ in range(conn):
        mol.AddBond(index, mol.AddAtom(Chem.Atom(0)), Chem.BondType.SINGLE)
    try:
        mol.UpdatePropertyCache(strict=True)
    except (ValueError, RuntimeError):  # valence error
        return INVALID
    return mol.GetAtomWithIdx(index).GetNumImplicitHs()


#: Real species as a sketcher would hold them (no hydrogen counts), the atom to check (index
#: in the SMILES), its true hydrogen count, and a name. Hydrogens in the SMILES are dropped
#: before testing: both models must work them out.
REFERENCE = [
    ("CC(C)(C)C", 1, 0, "neopentane, central C"),
    ("CC", 0, 3, "ethane, methyl C"),
    ("C[N+]#[C-]", 2, 0, "methyl isocyanide, terminal C"),
    ("c1ccccc1[C+](c1ccccc1)c1ccccc1", 6, 0, "trityl cation"),
    ("CN(C)C", 1, 0, "trimethylamine"),
    ("CN", 1, 2, "methylamine N"),
    ("C[N+](C)(C)C", 1, 0, "tetramethylammonium"),
    ("C[N+](C)C", 1, 1, "trimethylammonium"),
    ("[N+]", 0, 4, "ammonium"),
    ("CC(C)[N-]C(C)C", 3, 0, "diisopropylamide"),
    ("C[N-]", 1, 1, "methylamide anion"),
    ("O=N(=O)c1ccccc1", 1, 0, "nitrobenzene drawn pentavalent (CTfile convention)"),
    ("O=[N+]([O-])c1ccccc1", 1, 0, "nitrobenzene, charge-separated"),
    ("COC", 1, 0, "dimethyl ether"),
    ("CO", 1, 1, "methanol O"),
    ("C[O-]", 1, 0, "methoxide"),
    ("C[O+](C)C", 1, 0, "trimethyloxonium"),
    ("[O+]", 0, 3, "hydronium"),
    ("CB(C)C", 1, 0, "trimethylborane"),
    ("F[B-](F)(F)F", 1, 0, "tetrafluoroborate"),
    ("[B-]", 0, 4, "borohydride"),
    ("[Al-]", 0, 4, "tetrahydroaluminate"),
    ("C[Si](C)(C)C", 1, 0, "tetramethylsilane"),
    ("CC[Si](CC)CC", 2, 1, "triethylsilane Si"),
    ("F[Si-2](F)(F)(F)(F)F", 1, 0, "hexafluorosilicate"),
    ("CCCC[Sn](CCCC)CCCC", 4, 1, "tributyltin hydride Sn"),
    ("c1ccccc1P(c1ccccc1)c1ccccc1", 6, 0, "triphenylphosphine"),
    ("CP", 1, 2, "methylphosphine P"),
    ("OP(=O)(O)O", 1, 0, "phosphoric acid"),
    ("ClP(Cl)(Cl)(Cl)Cl", 1, 0, "phosphorus pentachloride"),
    ("F[P-](F)(F)(F)(F)F", 1, 0, "hexafluorophosphate"),
    ("F[As-](F)(F)(F)(F)F", 1, 0, "hexafluoroarsenate"),
    ("F[Sb-](F)(F)(F)(F)F", 1, 0, "hexafluoroantimonate"),
    ("CSC", 1, 0, "dimethyl sulfide"),
    ("CS", 1, 1, "methanethiol S"),
    ("CS(C)=O", 1, 0, "DMSO drawn S=O"),
    ("CS(C)(=O)=O", 1, 0, "dimethyl sulfone"),
    ("FS(F)(F)(F)(F)F", 1, 0, "sulfur hexafluoride"),
    ("C[S+](C)C", 1, 0, "trimethylsulfonium"),
    ("C[S-]", 1, 0, "methanethiolate"),
    ("C[Se]C", 1, 0, "dimethyl selenide"),
    ("CF", 1, 0, "fluoromethane F"),
    ("[F-]", 0, 0, "fluoride"),
    ("CCl", 1, 0, "chloromethane Cl"),
    ("Cl", 0, 1, "hydrogen chloride"),
    ("[Cl-]", 0, 0, "chloride"),
    ("FCl(F)F", 1, 0, "chlorine trifluoride"),
    ("OCl(=O)=O", 1, 0, "chloric acid drawn with Cl=O"),
    ("OCl(=O)(=O)=O", 1, 0, "perchloric acid drawn with Cl=O"),
    ("c1ccccc1[Br+]c1ccccc1", 6, 0, "diphenylbromonium"),
    ("Ic1ccccc1", 0, 0, "iodobenzene I"),
    ("CC(=O)OI(OC(C)=O)c1ccccc1", 4, 0, "(diacetoxyiodo)benzene"),
    ("FI(F)(F)(F)F", 1, 0, "iodine pentafluoride"),
    ("FI(F)(F)(F)(F)(F)F", 1, 0, "iodine heptafluoride"),
    ("c1ccccc1[I+]c1ccccc1", 6, 0, "diphenyliodonium"),
    ("I[I-]I", 1, 0, "triiodide, central I"),
    ("CCCC[Li]", 4, 0, "n-butyllithium Li"),
    ("[Na+]", 0, 0, "sodium ion"),
    ("C[Mg]Br", 1, 0, "methylmagnesium bromide Mg"),
    ("CC[Zn]CC", 2, 0, "diethylzinc"),
    ("N[Pt](N)(Cl)Cl", 1, 0, "cisplatin drawn with four single bonds"),
]


def _drawn(smiles: str) -> Chem.Mol:
    """*smiles* as a sketch: heavy atoms, charges and bonds; no hydrogen information."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    return mol


def _mark(value: int | None, want: int) -> str:
    return show(value) + ("" if value == want else " ✗")


def reference_table() -> tuple[str, int, int]:
    """Markdown rows scoring both models on REFERENCE, and how many each gets right.

    RDKit is asked the way chemsuite would ask it: the sketch is written as a V2000 molfile
    without hydrogen counts and read back with sanitisation, which normalises hypervalent
    drawings such as nitro and perchlorate.
    """
    rows = []
    k_ok = r_ok = 0
    for smiles, index, want, species in REFERENCE:
        drawn = _drawn(smiles)
        atom = drawn.GetAtomWithIdx(index)
        conn = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
        k = ketcher_hydrogens(atom.GetSymbol(), atom.GetFormalCharge(), conn, 0)
        read = Chem.MolFromMolBlock(Chem.MolToMolBlock(drawn, kekulize=False))
        r = INVALID if read is None else read.GetAtomWithIdx(index).GetTotalNumHs()
        k_ok += k == want
        r_ok += r == want

        rows.append(f"| {species} | `{smiles}` | {want} | {_mark(k, want)} | {_mark(r, want)} |")
    return "\n".join(rows), k_ok, r_ok


def show(value: int | None) -> str:
    return "error" if value is INVALID else str(value)


def main() -> int:
    RDLogger.DisableLog("rdApp.*")
    out = sys.stdout
    cases = disagree = 0
    by_element: dict[str, list[str]] = defaultdict(list)
    for label in ELEMENTS:
        for charge in (-2, -1, 0, 1, 2):
            for rad in (0, 1, 2):
                for conn in range(7):
                    cases += 1
                    k = ketcher_hydrogens(label, charge, conn, rad)
                    r = rdkit_hydrogens(label, charge, conn, rad)
                    if k != r:
                        disagree += 1
                        by_element[label].append(
                            f"| {label} | {charge:+d} | {rad} | {conn} | {show(k)} | {show(r)} |"
                        )
    rows, k_ok, r_ok = reference_table()
    out.write(
        f"## Reference species: Ketcher {k_ok}/{len(REFERENCE)}, RDKit {r_ok}/{len(REFERENCE)}\n\n"
    )
    out.write("| Species | Drawn as | Expected H | Ketcher | RDKit (molfile) |\n")
    out.write("|---|---|---|---|---|\n" + rows + "\n\n")
    out.write(f"## All combinations: {cases} cases, {disagree} disagreements\n\n")
    for label, rows in by_element.items():
        out.write(f"### {label}: {len(rows)}\n\n")
        out.write("| Element | Charge | Radical e⁻ | Bond-order sum | Ketcher H | RDKit H |\n")
        out.write("|---|---|---|---|---|---|\n")
        out.write("\n".join(rows) + "\n\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
