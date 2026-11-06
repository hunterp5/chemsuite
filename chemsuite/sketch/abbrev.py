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

"""Abbreviations (Ph, Boc, CO₂Me …) and typed atom labels (SKETCHER_PLAN §4.5).

The table is chemsuite's own. Neither candidate source could be shipped as is
(docs/sketcher/abbreviation-review.md):
- Ketcher's ``fg.sdf`` has eight records RDKit cannot read, a wrong Mes and a valence-less SO2.
- RDKit's default abbreviations define iBu as sec-butyl.

Each entry states the substituent's IUPAC name and the formula that name implies, written by
hand. tests/sketch/test_abbrev.py checks every structure against its formula and attachment
atom, and scripts/review/abbreviation_review.py compares the table with both sources.

An abbreviation is stored as real atoms plus a ``Group`` record, so the molecule, its valence
and its SMILES never depend on whether the label is shown contracted.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, replace
from statistics import median

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor

from chemsuite.sketch import geometry
from chemsuite.sketch.edits import Change
from chemsuite.sketch.model import ELEMENT_FIELDS, Atom, Group, Order, query_fields

Point = tuple[float, float]


@dataclass(frozen=True)
class Abbreviation:
    """``smiles`` starts with ``*``, the point of attachment, bonded to the attachment atom.

    ``formula`` is the substituent's formula (without the atom it attaches to), from its name.
    ``left`` is the label as written when the bond comes in from the right (MeO₂C-).
    """

    label: str
    smiles: str
    formula: str
    name: str
    aliases: tuple[str, ...] = ()
    left: str = ""

    @property
    def left_label(self) -> str:
        return self.left or self.label


A = Abbreviation
TABLE: tuple[Abbreviation, ...] = (
    # alkyl and aryl
    A("Me", "*C", "CH3", "methyl"),
    A("Et", "*CC", "C2H5", "ethyl"),
    A("Pr", "*CCC", "C3H7", "propyl", ("nPr", "n-Pr")),
    A("iPr", "*C(C)C", "C3H7", "propan-2-yl", ("i-Pr",)),
    A("cPr", "*C1CC1", "C3H5", "cyclopropyl", ("c-Pr",)),
    A("Bu", "*CCCC", "C4H9", "butyl", ("nBu", "n-Bu")),
    A("iBu", "*CC(C)C", "C4H9", "2-methylpropyl", ("i-Bu",)),
    A("sBu", "*C(C)CC", "C4H9", "butan-2-yl", ("s-Bu", "secBu")),
    A("tBu", "*C(C)(C)C", "C4H9", "tert-butyl", ("t-Bu", "tertBu")),
    A("Cy", "*C1CCCCC1", "C6H11", "cyclohexyl", ("cHex",)),
    A("Ph", "*c1ccccc1", "C6H5", "phenyl"),
    A("Bn", "*Cc1ccccc1", "C7H7", "benzyl"),
    A("Tol", "*c1ccc(C)cc1", "C7H7", "4-methylphenyl (p-tolyl)", ("pTol", "p-Tol")),
    A("Mes", "*c1c(C)cc(C)cc1C", "C9H11", "2,4,6-trimethylphenyl (mesityl)"),
    A("Tr", "*C(c1ccccc1)(c1ccccc1)c1ccccc1", "C19H15", "triphenylmethyl (trityl)", ("CPh3",)),
    A("PMB", "*Cc1ccc(OC)cc1", "C8H9O", "(4-methoxyphenyl)methyl"),
    # acyl and carbonyl groups
    A("Ac", "*C(C)=O", "C2H3O", "acetyl"),
    A("Bz", "*C(=O)c1ccccc1", "C7H5O", "benzoyl"),
    A("Piv", "*C(=O)C(C)(C)C", "C5H9O", "2,2-dimethylpropanoyl (pivaloyl)"),
    A("Boc", "*C(=O)OC(C)(C)C", "C5H9O2", "tert-butoxycarbonyl"),
    A("Cbz", "*C(=O)OCc1ccccc1", "C8H7O2", "benzyloxycarbonyl", ("Z",)),
    A("Fmoc", "*C(=O)OCC1c2ccccc2-c2ccccc21", "C15H11O2", "(9H-fluoren-9-ylmethoxy)carbonyl"),
    A("CHO", "*C=O", "CHO", "formyl", left="OHC"),
    A("CO2H", "*C(=O)O", "CHO2", "carboxy", ("COOH",), left="HO2C"),
    A("CO2Me", "*C(=O)OC", "C2H3O2", "methoxycarbonyl", ("COOMe",), left="MeO2C"),
    A("CO2Et", "*C(=O)OCC", "C3H5O2", "ethoxycarbonyl", ("COOEt",), left="EtO2C"),
    A("CO2tBu", "*C(=O)OC(C)(C)C", "C5H9O2", "tert-butoxycarbonyl", ("COOtBu",), left="tBuO2C"),
    A("CONH2", "*C(N)=O", "CH2NO", "carbamoyl", left="H2NOC"),
    A("CN", "*C#N", "CN", "cyano", left="NC"),
    # halogenated
    A("CF3", "*C(F)(F)F", "CF3", "trifluoromethyl", left="F3C"),
    A("CCl3", "*C(Cl)(Cl)Cl", "CCl3", "trichloromethyl", left="Cl3C"),
    # oxygen
    A("OMe", "*OC", "CH3O", "methoxy", left="MeO"),
    A("OEt", "*OCC", "C2H5O", "ethoxy", left="EtO"),
    A("OAc", "*OC(C)=O", "C2H3O2", "acetyloxy", left="AcO"),
    A("OBn", "*OCc1ccccc1", "C7H7O", "benzyloxy", left="BnO"),
    A("OCF3", "*OC(F)(F)F", "CF3O", "trifluoromethoxy", left="F3CO"),
    A("OTf", "*OS(=O)(=O)C(F)(F)F", "CF3O3S", "trifluoromethanesulfonyloxy", left="TfO"),
    A("OTs", "*OS(=O)(=O)c1ccc(C)cc1", "C7H7O3S", "4-methylbenzenesulfonyloxy", left="TsO"),
    A("OMs", "*OS(C)(=O)=O", "CH3O3S", "methanesulfonyloxy", left="MsO"),
    A("MOM", "*COC", "C2H5O", "methoxymethyl"),
    A("THP", "*C1CCCCO1", "C5H9O", "oxan-2-yl (tetrahydropyranyl)"),
    A("SEM", "*COCC[Si](C)(C)C", "C6H15OSi", "[2-(trimethylsilyl)ethoxy]methyl"),
    # nitrogen
    A("NO2", "*[N+](=O)[O-]", "NO2", "nitro", left="O2N"),
    A("NHAc", "*NC(C)=O", "C2H4NO", "acetamido", left="AcHN"),
    A("NHBoc", "*NC(=O)OC(C)(C)C", "C5H10NO2", "(tert-butoxycarbonyl)amino", left="BocHN"),
    A("NMe2", "*N(C)C", "C2H6N", "dimethylamino", left="Me2N"),
    A("NCO", "*N=C=O", "CNO", "isocyanato", left="OCN"),
    A("NCS", "*N=C=S", "CNS", "isothiocyanato", left="SCN"),
    # sulfur
    A("SMe", "*SC", "CH3S", "methylsulfanyl", left="MeS"),
    A("SCN", "*SC#N", "CNS", "thiocyanato", left="NCS"),
    A("Ms", "*S(C)(=O)=O", "CH3O2S", "methanesulfonyl", ("SO2Me",)),
    A("Ts", "*S(=O)(=O)c1ccc(C)cc1", "C7H7O2S", "4-methylbenzenesulfonyl (tosyl)", ("Tos",)),
    A("Tf", "*S(=O)(=O)C(F)(F)F", "CF3O2S", "trifluoromethanesulfonyl"),
    A("SO2NH2", "*S(N)(=O)=O", "H2NO2S", "sulfamoyl", left="H2NO2S"),
    A("SO2Cl", "*S(=O)(=O)Cl", "ClO2S", "chlorosulfonyl", left="ClO2S"),
    A("SO3H", "*S(=O)(=O)O", "HO3S", "sulfo", left="HO3S"),
    # silicon, boron, phosphorus
    A("TMS", "*[Si](C)(C)C", "C3H9Si", "trimethylsilyl"),
    A("TES", "*[Si](CC)(CC)CC", "C6H15Si", "triethylsilyl"),
    A("TBS", "*[Si](C)(C)C(C)(C)C", "C6H15Si", "tert-butyldimethylsilyl", ("TBDMS",)),
    A("TIPS", "*[Si](C(C)C)(C(C)C)C(C)C", "C9H21Si", "tri(propan-2-yl)silyl"),
    A("TBDPS", "*[Si](c1ccccc1)(c1ccccc1)C(C)(C)C", "C16H19Si", "tert-butyldiphenylsilyl"),
    A("Bpin", "*B1OC(C)(C)C(C)(C)O1", "C6H12BO2", "4,4,5,5-tetramethyl-1,3,2-dioxaborolan-2-yl"),
    A("PO3H2", "*P(=O)(O)O", "H2O3P", "phosphono", left="H2O3P"),
)  # fmt: skip
del A

_BY_KEY = {}
for _abbr in TABLE:
    for _key in (_abbr.label, *_abbr.aliases):
        _BY_KEY.setdefault(_key.lower(), _abbr)


_EXACT = {key: entry for entry in TABLE for key in (entry.label, *entry.aliases)}


def find(text: str) -> Abbreviation | None:
    """The abbreviation a label names (case-insensitive, aliases included)."""
    return _BY_KEY.get(text.strip().lower())


_ELEMENT_LABEL = re.compile(r"^(\d+)?([A-Z][a-z]?)(?:H(\d*))?([+-−]\d*|\d*[+-−])?$")


def parse_label(text: str):
    """What a typed label means: ``("atom", symbol, isotope, charge)``, ``("abbrev", entry)``,
    ``("query", fields)`` or None.

    Atom labels may carry an isotope, the hydrogens the atom would have anyway, and a charge:
    ``N``, ``Cl``, ``13C``, ``OH``, ``NH2``, ``N+``, ``O-``, ``D``.

    Order: an abbreviation written exactly (so ``Ac``, ``Pr`` and ``Ts`` are acetyl, propyl and
    tosyl, as in ChemDraw, not actinium, praseodymium and tennessine), then R-groups and query
    atoms (``R1``, ``A``, ``Q``, ``X``, ``[C,N,O]``; see ``model.query_fields``), then element
    symbols (so ``co`` is not taken for cobalt but ``Co`` is), then abbreviations in any case.
    """
    text = text.strip()
    exact = _EXACT.get(text)
    if exact is not None:
        return ("abbrev", exact)
    query = query_fields(text)
    if query is not None:
        return ("query", query)
    if text in ("D", "T"):
        return ("atom", "H", 2 if text == "D" else 3, 0)
    match = _ELEMENT_LABEL.match(text)
    if match:
        isotope, symbol, _h, charge = match.groups()
        if _is_element(symbol):
            return ("atom", symbol, int(isotope or 0), _charge(charge))
    entry = find(text)
    return ("abbrev", entry) if entry else None


def _is_element(symbol: str) -> bool:
    try:
        with rdBase.BlockLogs():
            return Chem.GetPeriodicTable().GetAtomicNumber(symbol) > 0
    except (RuntimeError, ValueError):
        return False


def _charge(text: str | None) -> int:
    if not text:
        return 0
    sign = 1 if "+" in text else -1
    digits = text.strip("+-−")
    return sign * int(digits or 1)


def display(label: str) -> list[tuple[str, str]]:
    """Label runs with digits after a letter or bracket subscripted (CO₂Me), as in a formula."""
    runs: list[tuple[str, str]] = []
    for part in re.findall(r"\d+|\D+", label):
        if part.isdigit() and runs:
            runs.append((part, "sub"))
        elif part.isdigit():
            runs.append((part, "sup"))  # a leading number is an isotope
        else:
            runs.append((part, "normal"))
    return runs


def _layout(entry: Abbreviation) -> tuple[Chem.Mol, list[Point]]:
    """The group as a Kekulé molecule with 2D positions in bond lengths; atom 0 is the ``*``."""
    mol = Chem.MolFromSmiles(entry.smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    # RDKit's own depictor, whatever the global preference (chemsuite.core.depict prefers
    # CoordGen, which crowds bulky groups such as trityl and TIPS onto each other).
    prefer = rdDepictor.GetPreferCoordGen()
    rdDepictor.SetPreferCoordGen(False)
    try:
        rdDepictor.Compute2DCoords(mol)
    finally:
        rdDepictor.SetPreferCoordGen(prefer)
    conf = mol.GetConformer()
    lengths = [
        (
            conf.GetAtomPosition(b.GetBeginAtomIdx()) - conf.GetAtomPosition(b.GetEndAtomIdx())
        ).Length()
        for b in mol.GetBonds()
    ]
    scale = median(lengths)
    points = [(conf.GetAtomPosition(i).x / scale, conf.GetAtomPosition(i).y / scale)
              for i in range(mol.GetNumAtoms())]  # fmt: skip
    return mol, points


_ORDERS = {Chem.BondType.DOUBLE: Order.DOUBLE, Chem.BondType.TRIPLE: Order.TRIPLE}


def place(change: Change, atom: int, entry: Abbreviation) -> int:
    """Turn ``atom`` into the abbreviation: it becomes the group's attachment atom and the rest
    of the group is laid out away from its neighbour. Returns the group id.

    Raises ValueError if the atom has more than one bond (an abbreviation has one attachment).
    """
    sketch = change.sketch
    neighbours = sketch.neighbours(atom)
    if len(neighbours) > 1:
        raise ValueError(f"{entry.label} attaches by one bond; this atom has {len(neighbours)}.")
    mol, points = _layout(entry)
    star = mol.GetAtomWithIdx(0)
    (attach_rd,) = [n.GetIdx() for n in star.GetNeighbors()]
    here = sketch.atoms[atom].pos
    want = geometry.angle_of(sketch.atoms[neighbours[0]].pos, here) if neighbours else 0.0
    have = geometry.angle_of(points[0], points[attach_rd])
    turn = want - have
    c, s = math.cos(turn), math.sin(turn)
    ax, ay = points[attach_rd]
    ids: dict[int, int] = {}
    for rd in mol.GetAtoms():
        i = rd.GetIdx()
        if i == 0:
            continue
        x, y = points[i][0] - ax, points[i][1] - ay
        pos = (here[0] + x * c - y * s, here[1] + x * s + y * c)
        value = Atom(*pos, symbol=rd.GetSymbol(), charge=rd.GetFormalCharge())
        if i == attach_rd:
            old = sketch.atoms[atom]
            change.set_atom(atom, replace(value, x=old.x, y=old.y))
            ids[i] = atom
        else:
            ids[i] = change.add_atom(value)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if 0 in (a, b):
            continue
        change.add_bond(ids[a], ids[b], _ORDERS.get(bond.GetBondType(), Order.SINGLE))
    return change.add_group(Group(entry.label, frozenset(ids.values()), atom))


def set_label(change: Change, atom: int, text: str) -> str:
    """Apply a typed label to an atom. Returns "" or why the label could not be used."""
    parsed = parse_label(text)
    if parsed is None:
        return f"“{text}” is not an element, query atom or known abbreviation."
    if parsed[0] == "query":
        old = change.sketch.atoms[atom]
        new = replace(old, **parsed[1])
        if new != old:
            change.set_atom(atom, new)
        return ""
    if parsed[0] == "abbrev":
        try:
            place(change, atom, parsed[1])
        except ValueError as err:
            return str(err)
        return ""
    _, symbol, isotope, charge = parsed
    old = change.sketch.atoms[atom]
    change.set_atom(
        atom, replace(old, symbol=symbol, isotope=isotope, charge=charge, **ELEMENT_FIELDS)
    )
    return ""


def set_expanded(change: Change, gid: int, expanded: bool) -> None:
    group = change.sketch.groups[gid]
    change.set_group(gid, replace(group, expanded=expanded))


def ungroup(change: Change, gid: int) -> None:
    """Forget the abbreviation; its atoms stay as drawn."""
    change.set_group(gid, None)
