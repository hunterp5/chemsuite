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

"""Salt forms: the fragments a strip removes, and making a salt of a parent structure.

Stripping keeps the largest organic fragment; :func:`removed_fragments` says what went, so a
counter-ion is recorded rather than lost.

:func:`form_salt` goes the other way. A salt is charge balanced, so the proton moves: an
anionic counter-ion (chloride, mesylate) takes the parent's most basic site to a cation, and a
cationic one (sodium, choline) takes its most acidic site to an anion. The site is chosen by
the SMARTS order in :data:`BASIC_SITES` and :data:`ACIDIC_SITES`, strongest first, which
follows the usual ordering of these groups in water; it is a rule of thumb, not a pKa
calculation (Predict pKa does that). A parent that already carries the needed charge keeps it
and only gains the counter-ion. A parent with no suitable site gives None, so the caller can
say why.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from rdkit import Chem, rdBase

#: Where a proton is taken from the parent, strongest base first.
BASIC_SITES: tuple[tuple[str, str], ...] = (
    ("guanidine", "[NX3;H2,H1,H0;!$(N[#6]=[O,S]);!$(N[S,P]=O)][CX3](=[NX2])[NX3]"),
    ("amidine", "[CX3](=[NX2;!$(N[#8]):1])[NX3;!$(N[#6]=[O,S])]"),
    (
        "aliphatic amine",
        "[NX3;H2,H1,H0;!$(N[#6]=[O,S,N]);!$(N[S,P]=O);!$(N[#8]);!$(N[#7]);!$(N=*);!$(Nc);!$(N#*)]",
    ),  # fmt: skip
    ("aromatic nitrogen", "[nX2;!$(n[#8]);!$(n[#7])]"),
    ("aniline", "[NX3;H2,H1;!$(N[#6]=[O,S,N]);!$(N[S,P]=O);!$(N[#8]);$(Nc)]"),
)
#: Where a proton is given up to the counter-ion, strongest acid first.
ACIDIC_SITES: tuple[tuple[str, str], ...] = (
    ("sulfonic acid", "[SX4](=[OX1])(=[OX1])[OX2H1]"),
    ("carboxylic acid", "[CX3](=[OX1])[OX2H1]"),
    ("tetrazole", "[nX3H1,NX3H1]1[nX2,NX2][nX2,NX2][nX2,NX2][cX3,CX3]1"),
    ("acyl sulfonamide", "[SX4](=[OX1])(=[OX1])[NX3H1][CX3]=[OX1]"),
    ("phosphonic acid", "[PX4](=[OX1])([OX2H1])"),
    ("phenol", "[OX2H1]c"),
)


@dataclass(frozen=True)
class Counterion:
    """A counter-ion that can be added to a parent structure.

    *smiles* is the ion as it stands in the salt. *charge* is its charge, so a parent needs
    ``abs(charge)`` protons moved for each copy. *acid* says which way the proton goes: an
    acid's conjugate base (chloride) protonates the parent, a base's conjugate acid (sodium)
    deprotonates it.
    """

    id: str
    name: str
    smiles: str
    charge: int
    acid: bool
    #: the everyday abbreviation, shown beside the name ("HCl", "TFA")
    short: str = ""

    @property
    def protons(self) -> int:
        """How many sites of the parent one copy of this ion changes."""
        return abs(self.charge)

    @property
    def label(self) -> str:
        """The name as a menu shows it."""
        return f"{self.name} ({self.short})" if self.short else self.name


COUNTERIONS: tuple[Counterion, ...] = (
    Counterion("hydrochloride", "Hydrochloride", "[Cl-]", -1, True, "HCl"),
    Counterion("hydrobromide", "Hydrobromide", "[Br-]", -1, True, "HBr"),
    Counterion("sulfate", "Sulfate", "[O-]S(=O)(=O)[O-]", -2, True),
    Counterion("phosphate", "Phosphate", "[O-]P(=O)([O-])[O-]", -3, True),
    Counterion("mesylate", "Mesylate", "CS(=O)(=O)[O-]", -1, True),
    Counterion("tosylate", "Tosylate", "Cc1ccc(S(=O)(=O)[O-])cc1", -1, True),
    Counterion("besylate", "Besylate", "[O-]S(=O)(=O)c1ccccc1", -1, True),
    Counterion("trifluoroacetate", "Trifluoroacetate", "[O-]C(=O)C(F)(F)F", -1, True, "TFA"),
    Counterion("acetate", "Acetate", "CC(=O)[O-]", -1, True),
    Counterion("citrate", "Citrate", "[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]", -3, True),
    Counterion("fumarate", "Fumarate", "[O-]C(=O)/C=C/C(=O)[O-]", -2, True),
    Counterion("maleate", "Maleate", "[O-]C(=O)/C=C\\C(=O)[O-]", -2, True),
    Counterion("succinate", "Succinate", "[O-]C(=O)CCC(=O)[O-]", -2, True),
    Counterion("tartrate", "Tartrate", "[O-]C(=O)[C@H](O)[C@@H](O)C(=O)[O-]", -2, True),
    Counterion("oxalate", "Oxalate", "[O-]C(=O)C(=O)[O-]", -2, True),
    Counterion("sodium", "Sodium", "[Na+]", 1, False),
    Counterion("potassium", "Potassium", "[K+]", 1, False),
    Counterion("lithium", "Lithium", "[Li+]", 1, False),
    Counterion("calcium", "Calcium", "[Ca+2]", 2, False),
    Counterion("magnesium", "Magnesium", "[Mg+2]", 2, False),
    Counterion("zinc", "Zinc", "[Zn+2]", 2, False),
    Counterion("ammonium", "Ammonium", "[NH4+]", 1, False),
    Counterion("choline", "Choline", "C[N+](C)(C)CCO", 1, False),
    Counterion("meglumine", "Meglumine", "C[NH2+]C[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO", 1, False),
    Counterion("tromethamine", "Tromethamine", "OCC([NH3+])(CO)CO", 1, False, "Tris"),
)

_BY_ID = {ion.id: ion for ion in COUNTERIONS}


def by_id(ion_id: str) -> Counterion | None:
    return _BY_ID.get(ion_id)


@cache
def _pattern(smarts: str) -> Chem.Mol | None:
    with rdBase.BlockLogs():
        return Chem.MolFromSmarts(smarts)


def _site(mol: Chem.Mol, sites: tuple[tuple[str, str], ...], taken: set[int]) -> int | None:
    """The first atom matching *sites* in order, skipping atoms already used."""
    for _name, smarts in sites:
        query = _pattern(smarts)
        if query is None:
            continue
        for match in mol.GetSubstructMatches(query):
            for index in match:
                atom = mol.GetAtomWithIdx(index)
                if index in taken or atom.GetAtomicNum() not in (7, 8, 15, 16):
                    continue
                if atom.GetAtomicNum() in (8, 15, 16) and atom.GetTotalNumHs() == 0:
                    continue  # an acidic site needs its proton
                return index
    return None


def _protonate_site(mol: Chem.Mol, index: int) -> None:
    atom = mol.GetAtomWithIdx(index)
    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
    atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
    atom.SetNoImplicit(True)


def _deprotonate_site(mol: Chem.Mol, index: int) -> None:
    atom = mol.GetAtomWithIdx(index)
    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
    atom.SetNumExplicitHs(max(atom.GetTotalNumHs() - 1, 0))
    atom.SetNoImplicit(True)


def _charged_parent(mol: Chem.Mol, ion: Counterion, wanted: int) -> tuple[Chem.Mol, int] | None:
    """*mol* with up to *wanted* protons moved, and how many moved. None: no site at all."""
    work = Chem.RWMol(mol)
    sites = BASIC_SITES if ion.acid else ACIDIC_SITES
    taken: set[int] = set()
    moved = 0
    while moved < wanted:
        index = _site(work, sites, taken)
        if index is None:
            break
        taken.add(index)
        if ion.acid:
            _protonate_site(work, index)
        else:
            _deprotonate_site(work, index)
        moved += 1
    if not moved:
        return None
    return work.GetMol(), moved


def form_salt(mol: Chem.Mol, ion: Counterion, count: int = 1) -> Chem.Mol | None:
    """*mol* as a salt of *count* copies of *ion*, or None when it has no suitable site.

    The salt is charge balanced: each proton the ion needs is taken from, or given to, the
    strongest site :data:`BASIC_SITES` or :data:`ACIDIC_SITES` finds, and charge the parent
    already carries counts towards the balance. A multiply charged ion that the parent cannot
    balance on its own repeats the parent, so ethylamine and sulfate give the 2:1 salt.
    """
    if count < 1:
        return None
    counterion = Chem.MolFromSmiles(ion.smiles)
    if counterion is None:
        return None
    wanted = ion.protons * count
    have = Chem.GetFormalCharge(mol)
    needed = wanted - (have if ion.acid else -have)
    if needed <= 0:  # the parent already balances the ion
        parent, copies, moved = Chem.Mol(mol), 1, 0
    else:
        charged = _charged_parent(mol, ion, needed)
        if charged is None:
            return None
        parent, moved = charged
        copies = -(-needed // moved)  # repeat the parent until the charges balance
        if copies * moved != needed:  # a part charge would be left over
            return None
    out = parent
    for _ in range(copies - 1):
        out = Chem.CombineMols(out, parent)
    for _ in range(count):
        out = Chem.CombineMols(out, counterion)
    try:
        Chem.SanitizeMol(out)
    except (ValueError, RuntimeError):
        return None
    return out


def removed_fragments(before: Chem.Mol, after: Chem.Mol) -> str:
    """The fragments of *before* that *after* no longer has, as dot-joined SMILES.

    One copy is removed for each fragment kept, so a symmetrical disalt still records the
    copy that went.
    """
    with rdBase.BlockLogs():
        was = [
            Chem.MolToSmiles(f) for f in Chem.GetMolFrags(before, asMols=True, sanitizeFrags=False)
        ]
        kept = [
            Chem.MolToSmiles(f) for f in Chem.GetMolFrags(after, asMols=True, sanitizeFrags=False)
        ]
    for smiles in kept:
        if smiles in was:
            was.remove(smiles)
    return ".".join(was)


__all__ = [
    "ACIDIC_SITES",
    "BASIC_SITES",
    "COUNTERIONS",
    "Counterion",
    "by_id",
    "form_salt",
    "removed_fragments",
]
