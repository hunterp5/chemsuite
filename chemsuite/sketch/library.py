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
"""The sketcher's template library (Shift+T): ring systems a medicinal chemist draws daily.

Ketcher ships a ``library.sdf`` of 276 templates; this library is written afresh rather than
copied (SKETCHER_PLAN §1): each entry is a SMILES, laid out by the sketcher's own Clean
(``layout.py``), so every template has the same bond length and look as the rest of a drawing.
The tests check that every template, once placed, is the molecule its SMILES names (by
InChIKey), that hand-worked formulas match, and that no layout puts atoms on top of each other.

**Attachment.** The SMILES's first atom is where the template attaches when dropped on an atom
(the 4-position of pyridine, C5 of pyrimidine, C5 of indole…), and the first ring bond at that
atom is the one laid onto a clicked bond when fusing.

Templates are placed by the ring tool's rules (``templates.py``): on empty canvas, on an atom
(spiro or substituent), or fused onto a bond.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from rdkit import Chem

from chemsuite.sketch.templates import Template

HETEROAROMATIC = "Aromatic heterocycles"
SATURATED = "Saturated heterocycles"
FUSED = "Fused bicycles"
BRIDGED = "Bridged and spiro"


@dataclass(frozen=True)
class Entry:
    name: str
    smiles: str
    group: str
    formula: str = ""  # worked out by hand (checked by the tests where given)


def _e(name: str, smiles: str, group: str, formula: str = "") -> Entry:
    return Entry(name, smiles, group, formula)


# fmt: off
TABLE: tuple[Entry, ...] = (
    _e("Pyridine", "c1ccncc1", HETEROAROMATIC, "C5H5N"),
    _e("Pyrimidine", "c1cncnc1", HETEROAROMATIC, "C4H4N2"),
    _e("Pyrazine", "c1cnccn1", HETEROAROMATIC, "C4H4N2"),
    _e("Pyridazine", "c1ccnnc1", HETEROAROMATIC, "C4H4N2"),
    _e("1,3,5-Triazine", "c1ncncn1", HETEROAROMATIC, "C3H3N3"),
    _e("Pyrrole", "c1cc[nH]c1", HETEROAROMATIC, "C4H5N"),
    _e("Furan", "c1ccoc1", HETEROAROMATIC, "C4H4O"),
    _e("Thiophene", "c1ccsc1", HETEROAROMATIC, "C4H4S"),
    _e("Imidazole", "c1cnc[nH]1", HETEROAROMATIC, "C3H4N2"),
    _e("Pyrazole", "c1cn[nH]c1", HETEROAROMATIC, "C3H4N2"),
    _e("Oxazole", "c1cocn1", HETEROAROMATIC, "C3H3NO"),
    _e("Isoxazole", "c1cnoc1", HETEROAROMATIC, "C3H3NO"),
    _e("Thiazole", "c1cscn1", HETEROAROMATIC, "C3H3NS"),
    _e("Isothiazole", "c1cnsc1", HETEROAROMATIC),
    _e("1,2,3-Triazole", "c1c[nH]nn1", HETEROAROMATIC, "C2H3N3"),
    _e("1,2,4-Triazole", "c1nc[nH]n1", HETEROAROMATIC, "C2H3N3"),
    _e("Tetrazole", "c1nnn[nH]1", HETEROAROMATIC, "CH2N4"),
    _e("1,2,4-Oxadiazole", "c1ncon1", HETEROAROMATIC, "C2H2N2O"),
    _e("1,3,4-Oxadiazole", "c1nnco1", HETEROAROMATIC, "C2H2N2O"),
    _e("1,3,4-Thiadiazole", "c1nncs1", HETEROAROMATIC),
    _e("2-Pyridone", "C1=CC(=O)NC=C1", HETEROAROMATIC, "C5H5NO"),

    _e("Azetidine", "C1CNC1", SATURATED, "C3H7N"),
    _e("Oxetane", "C1COC1", SATURATED, "C3H6O"),
    _e("Pyrrolidine", "C1CCNC1", SATURATED, "C4H9N"),
    _e("Tetrahydrofuran", "C1CCOC1", SATURATED, "C4H8O"),
    _e("Piperidine", "C1CCNCC1", SATURATED, "C5H11N"),
    _e("Piperazine", "C1CNCCN1", SATURATED, "C4H10N2"),
    _e("Morpholine", "C1COCCN1", SATURATED, "C4H9NO"),
    _e("Thiomorpholine", "C1CSCCN1", SATURATED),
    _e("Tetrahydropyran", "C1CCOCC1", SATURATED, "C5H10O"),
    _e("1,4-Dioxane", "C1COCCO1", SATURATED, "C4H8O2"),
    _e("Azepane", "C1CCCNCC1", SATURATED, "C6H13N"),
    _e("Piperidin-4-one", "C1CC(=O)CCN1", SATURATED, "C5H9NO"),
    _e("2-Pyrrolidinone", "C1CC(=O)NC1", SATURATED, "C4H7NO"),

    _e("Naphthalene", "c1ccc2ccccc2c1", FUSED, "C10H8"),
    _e("Indole", "c1ccc2[nH]ccc2c1", FUSED, "C8H7N"),
    _e("7-Azaindole", "c1cnc2[nH]ccc2c1", FUSED, "C7H6N2"),
    _e("Indazole", "c1ccc2[nH]ncc2c1", FUSED, "C7H6N2"),
    _e("Benzimidazole", "c1ccc2[nH]cnc2c1", FUSED, "C7H6N2"),
    _e("Benzofuran", "c1ccc2occc2c1", FUSED, "C8H6O"),
    _e("Benzothiophene", "c1ccc2sccc2c1", FUSED, "C8H6S"),
    _e("Benzoxazole", "c1ccc2ocnc2c1", FUSED, "C7H5NO"),
    _e("Benzothiazole", "c1ccc2scnc2c1", FUSED, "C7H5NS"),
    _e("Quinoline", "c1ccc2ncccc2c1", FUSED, "C9H7N"),
    _e("Isoquinoline", "c1ccc2cnccc2c1", FUSED, "C9H7N"),
    _e("Quinazoline", "c1ccc2ncncc2c1", FUSED, "C8H6N2"),
    _e("Quinoxaline", "c1ccc2nccnc2c1", FUSED, "C8H6N2"),
    _e("Purine", "c1ncc2[nH]cnc2n1", FUSED, "C5H4N4"),
    _e("Indane", "c1ccc2c(c1)CCC2", FUSED, "C9H10"),
    _e("Tetralin", "c1ccc2c(c1)CCCC2", FUSED, "C10H12"),
    _e("Indoline", "c1ccc2c(c1)CCN2", FUSED, "C8H9N"),
    _e("1,2,3,4-Tetrahydroisoquinoline", "c1ccc2c(c1)CCNC2", FUSED, "C9H11N"),
    _e("Chromane", "c1ccc2c(c1)CCCO2", FUSED, "C9H10O"),
    _e("1,3-Benzodioxole", "c1ccc2c(c1)OCO2", FUSED, "C7H6O2"),

    _e("Bicyclo[1.1.1]pentane", "C1C2CC1C2", BRIDGED, "C5H8"),
    _e("Bicyclo[2.2.2]octane", "C1CC2CCC1CC2", BRIDGED, "C8H14"),
    _e("Norbornane", "C1CC2CCC1C2", BRIDGED, "C7H12"),
    _e("Adamantane", "C1C2CC3CC1CC(C2)C3", BRIDGED, "C10H16"),
    _e("Quinuclidine", "C1CN2CCC1CC2", BRIDGED, "C7H13N"),
    _e("8-Azabicyclo[3.2.1]octane", "C1CC2CCC(C1)N2", BRIDGED, "C7H13N"),
    _e("3-Azabicyclo[3.1.0]hexane", "C1C2CC2CN1", BRIDGED, "C5H9N"),
    _e("Spiro[3.3]heptane", "C1CC2(C1)CCC2", BRIDGED, "C7H12"),
    _e("2-Azaspiro[3.3]heptane", "C1CC2(C1)CNC2", BRIDGED, "C6H11N"),
    _e("2-Oxa-6-azaspiro[3.3]heptane", "C1NCC12COC2", BRIDGED, "C5H9NO"),
)
# fmt: on

GROUPS: tuple[str, ...] = tuple(dict.fromkeys(e.group for e in TABLE))
BY_NAME = {e.name: e for e in TABLE}


@cache
def template(name: str) -> Template:
    """The library entry *name* as a placeable ``Template`` (laid out once, then cached)."""
    from chemsuite.sketch import convert, layout
    from chemsuite.sketch.edits import History

    entry = BY_NAME[name]
    mol = Chem.MolFromSmiles(entry.smiles)
    aromatic = any(a.GetIsAromatic() for a in mol.GetAtoms())
    Chem.Kekulize(mol, clearAromaticFlags=True)  # drawn Kekulé; re-assigned on fusing
    sketch = convert.from_mol(mol)
    history = History(sketch)
    change = history.change("Clean")
    layout.clean(change)
    history.commit(change)
    sketch = history.sketch
    ids = sorted(sketch.atoms)  # from_mol keeps the SMILES atom order
    index = {atom_id: i for i, atom_id in enumerate(ids)}
    atoms = tuple((sketch.atoms[i].x, sketch.atoms[i].y, sketch.atoms[i].symbol) for i in ids)
    bonds = tuple(
        (index[b.begin], index[b.end], b.order) for _bond_id, b in sorted(sketch.bonds.items())
    )
    ring_info = mol.GetRingInfo()
    attach_bond = next(
        (
            k
            for k, (i, j, _order) in enumerate(bonds)
            if 0 in (i, j) and ring_info.NumBondRings(mol.GetBondBetweenAtoms(i, j).GetIdx())
        ),
        0,
    )
    return Template(name, atoms, bonds, attach_atom=0, attach_bond=attach_bond, aromatic=aromatic)


def search(text: str) -> list[Entry]:
    """Entries whose name contains *text* (any case)."""
    needle = text.strip().lower()
    return [e for e in TABLE if needle in e.name.lower()]
