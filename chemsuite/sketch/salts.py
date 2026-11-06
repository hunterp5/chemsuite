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

"""Salts and solvents the sketcher can drop beside a drawing (SKETCHER_PLAN P2).

The list follows the 135 records of Ketcher's ``salts-and-solvents.sdf`` (commit 2b031cc2),
but every structure here is written out afresh and checked against PubChem by InChIKey
(``scripts/review/salts_review.py``, docs/sketcher/salts-review.md). Two Ketcher records are
drawn wrongly and are corrected here: its DBU has an N–N bond and a C=C bond (a different
compound, under the wrong name 1,5-diaza…), and its guanidine hydrochloride draws the cation
protonated on an amino nitrogen, H₃N⁺–C(=NH)NH₂, instead of as guanidinium, (H₂N)₂C=NH₂⁺
(InChI cannot tell the two apart, so the review flags only DBU). Its "perfluorocyclic ether"
is named for what it is, perfluoro(2-butyltetrahydrofuran).

``pubchem`` is the name (or CAS number) the review looks up; ``label`` is the short name
shown in menus where chemists use one (DMF, THF).
"""

from __future__ import annotations

from dataclasses import dataclass

from chemsuite.sketch import convert
from chemsuite.sketch.edits import Change, insert
from chemsuite.sketch.model import Sketch


@dataclass(frozen=True)
class Entry:
    name: str
    smiles: str
    group: str
    label: str = ""
    pubchem: str = ""

    @property
    def title(self) -> str:
        return f"{self.name} ({self.label})" if self.label else self.name


def _e(name, smiles, group, label="", pubchem=""):
    return Entry(name, smiles, group, label, pubchem or name)


ACIDS, ALCOHOLS, SMALL = "Acids and anhydrides", "Alcohols and diols", "Water and gases"
AROMATICS, BASES, APROTIC = "Aromatic hydrocarbons", "Amines and bases", "Polar aprotic"
ESTERS, ETHERS, HALO = "Esters and carbonates", "Ethers and glycols", "Halogenated"
ALKANES, KETONES, SALTS = "Alkanes", "Ketones", "Salts and buffers"

# fmt: off
TABLE: tuple[Entry, ...] = (
    _e("acetic acid", "CC(=O)O", ACIDS, "AcOH"),
    _e("acetic anhydride", "CC(=O)OC(C)=O", ACIDS, "Ac2O"),
    _e("formic acid", "O=CO", ACIDS),
    _e("methanesulfonic acid", "CS(=O)(=O)O", ACIDS, "MsOH"),
    _e("propionic acid", "CCC(=O)O", ACIDS),
    _e("chloroacetic acid", "OC(=O)CCl", ACIDS),
    _e("trichloroacetic acid", "OC(=O)C(Cl)(Cl)Cl", ACIDS, "TCA"),
    _e("trifluoroacetic acid", "OC(=O)C(F)(F)F", ACIDS, "TFA"),
    _e("1,2-propanediol", "CC(O)CO", ALCOHOLS),
    _e("1,3-propanediol", "OCCCO", ALCOHOLS),
    _e("1,4-butanediol", "OCCCCO", ALCOHOLS),
    _e("1-butanol", "CCCCO", ALCOHOLS, "n-BuOH"),
    _e("1-propanol", "CCCO", ALCOHOLS, "n-PrOH"),
    _e("2-butanol", "CCC(C)O", ALCOHOLS),
    _e("2-ethylhexanol", "CCCCC(CC)CO", ALCOHOLS, pubchem="2-ethyl-1-hexanol"),
    _e("isopropanol", "CC(C)O", ALCOHOLS, "IPA"),
    _e("2-methoxyethanol", "COCCO", ALCOHOLS),
    _e("2-pentanol", "CCCC(C)O", ALCOHOLS),
    _e("benzyl alcohol", "OCc1ccccc1", ALCOHOLS),
    _e("cyclohexanol", "OC1CCCCC1", ALCOHOLS),
    _e("ethanol", "CCO", ALCOHOLS, "EtOH"),
    _e("ethylene glycol", "OCCO", ALCOHOLS),
    _e("glycerol", "OCC(O)CO", ALCOHOLS),
    _e("isoamyl alcohol", "CC(C)CCO", ALCOHOLS),
    _e("isobutanol", "CC(C)CO", ALCOHOLS),
    _e("methanol", "CO", ALCOHOLS, "MeOH"),
    _e("tert-butanol", "CC(C)(C)O", ALCOHOLS, "t-BuOH"),
    _e("2,2,2-trifluoroethanol", "OCC(F)(F)F", ALCOHOLS, "TFE"),
    _e("water", "O", SMALL),
    _e("heavy water", "[2H]O[2H]", SMALL, "D2O", "deuterium oxide"),
    _e("carbon dioxide", "O=C=O", SMALL),
    _e("carbon disulfide", "S=C=S", SMALL),
    _e("benzene", "c1ccccc1", AROMATICS),
    _e("cumene", "CC(C)c1ccccc1", AROMATICS),
    _e("mesitylene", "Cc1cc(C)cc(C)c1", AROMATICS),
    _e("toluene", "Cc1ccccc1", AROMATICS),
    _e("o-xylene", "Cc1ccccc1C", AROMATICS),
    _e("m-xylene", "Cc1cccc(C)c1", AROMATICS),
    _e("p-xylene", "Cc1ccc(C)cc1", AROMATICS),
    _e("N,N-dimethylaniline", "CN(C)c1ccccc1", BASES),
    _e("pyridine", "c1ccncc1", BASES),
    _e("triethylamine", "CCN(CC)CC", BASES, "Et3N"),
    _e("N,N-diisopropylethylamine", "CCN(C(C)C)C(C)C", BASES, "DIPEA"),
    _e("1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1", BASES, "DBU"),
    _e("tris(hydroxymethyl)aminomethane", "NC(CO)(CO)CO", BASES, "TRIS",
       "tromethamine"),
    _e("acetonitrile", "CC#N", APROTIC, "MeCN"),
    _e("N,N-dimethylacetamide", "CC(=O)N(C)C", APROTIC, "DMA"),
    _e("N,N-dimethylformamide", "CN(C)C=O", APROTIC, "DMF"),
    _e("dimethyl sulfoxide", "CS(C)=O", APROTIC, "DMSO"),
    _e("N,N′-dimethylpropyleneurea", "CN1CCCN(C)C1=O", APROTIC, "DMPU",
       "1,3-dimethyl-3,4,5,6-tetrahydro-2(1H)-pyrimidinone"),
    _e("formamide", "NC=O", APROTIC),
    _e("N-methyl-2-pyrrolidone", "CN1CCCC1=O", APROTIC, "NMP"),
    _e("N-methylformamide", "CNC=O", APROTIC),
    _e("propionitrile", "CCC#N", APROTIC),
    _e("sulfolane", "O=S1(=O)CCCC1", APROTIC),
    _e("hexamethylphosphoramide", "CN(C)P(=O)(N(C)C)N(C)C", APROTIC, "HMPA"),
    _e("hexamethylphosphorous triamide", "CN(C)P(N(C)C)N(C)C", APROTIC, "HMPT"),
    _e("nitromethane", "C[N+](=O)[O-]", APROTIC),
    _e("amyl acetate", "CCCCCOC(C)=O", ESTERS, pubchem="pentyl acetate"),
    _e("butyl acetate", "CCCCOC(C)=O", ESTERS),
    _e("dimethyl carbonate", "COC(=O)OC", ESTERS),
    _e("ethyl acetate", "CCOC(C)=O", ESTERS, "EtOAc"),
    _e("ethyl formate", "CCOC=O", ESTERS),
    _e("ethyl lactate", "CCOC(=O)C(C)O", ESTERS),
    _e("ethyl propionate", "CCOC(=O)CC", ESTERS),
    _e("ethylene carbonate", "O=C1OCCO1", ESTERS),
    _e("isobutyl acetate", "CC(=O)OCC(C)C", ESTERS),
    _e("isopropyl acetate", "CC(=O)OC(C)C", ESTERS, "IPAc"),
    _e("methyl acetate", "COC(C)=O", ESTERS),
    _e("methyl formate", "COC=O", ESTERS),
    _e("methyl lactate", "COC(=O)C(C)O", ESTERS),
    _e("octyl acetate", "CCCCCCCCOC(C)=O", ESTERS),
    _e("propyl acetate", "CCCOC(C)=O", ESTERS),
    _e("propylene carbonate", "CC1COC(=O)O1", ESTERS),
    _e("tert-butyl acetate", "CC(=O)OC(C)(C)C", ESTERS),
    _e("1,2-dimethoxyethane", "COCCOC", ETHERS, "DME"),
    _e("1,4-dioxane", "C1COCCO1", ETHERS),
    _e("2-methyltetrahydrofuran", "CC1CCCO1", ETHERS, "2-MeTHF"),
    _e("anisole", "COc1ccccc1", ETHERS),
    _e("butyl carbitol", "CCCCOCCOCCO", ETHERS, pubchem="diethylene glycol monobutyl ether"),
    _e("cyclopentyl methyl ether", "COC1CCCC1", ETHERS, "CPME"),
    _e("dibutyl ether", "CCCCOCCCC", ETHERS),
    _e("diethyl ether", "CCOCC", ETHERS, "Et2O"),
    _e("diethylene glycol", "OCCOCCO", ETHERS),
    _e("diglyme", "COCCOCCOC", ETHERS),
    _e("diisopropyl ether", "CC(C)OC(C)C", ETHERS),
    _e("dimethyl ether", "COC", ETHERS),
    _e("diphenyl ether", "c1ccc(Oc2ccccc2)cc1", ETHERS),
    _e("ethoxybenzene", "CCOc1ccccc1", ETHERS),
    _e("methyl tert-butyl ether", "COC(C)(C)C", ETHERS, "MTBE"),
    _e("tert-amyl methyl ether", "CCC(C)(C)OC", ETHERS, "TAME"),
    _e("ethyl tert-butyl ether", "CCOC(C)(C)C", ETHERS, "ETBE"),
    _e("tetrahydrofuran", "C1CCOC1", ETHERS, "THF"),
    _e("triethylene glycol", "OCCOCCOCCO", ETHERS),
    _e("1,2,4-trichlorobenzene", "Clc1ccc(Cl)c(Cl)c1", HALO),
    _e("1,2-dichlorobenzene", "Clc1ccccc1Cl", HALO),
    _e("1,2-dichloroethane", "ClCCCl", HALO, "DCE"),
    _e("2,3,4-trifluorotoluene", "Cc1ccc(F)c(F)c1F", HALO),
    _e("carbon tetrachloride", "ClC(Cl)(Cl)Cl", HALO),
    _e("chlorobenzene", "Clc1ccccc1", HALO),
    _e("chloroform", "ClC(Cl)Cl", HALO),
    _e("dichloromethane", "ClCCl", HALO, "DCM"),
    _e("fluorobenzene", "Fc1ccccc1", HALO),
    _e("perfluoro(2-butyltetrahydrofuran)",
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)OC(F)(F)C(F)(F)C1(F)F", HALO,
       pubchem="335-36-4"),
    _e("perfluorocyclohexane", "FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F", HALO,
       pubchem="dodecafluorocyclohexane"),
    _e("perfluorohexane", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", HALO),
    _e("perfluorotoluene", "Fc1c(F)c(F)c(C(F)(F)F)c(F)c1F", HALO, pubchem="octafluorotoluene"),
    _e("trichloroacetonitrile", "N#CC(Cl)(Cl)Cl", HALO),
    _e("(trifluoromethyl)benzene", "FC(F)(F)c1ccccc1", HALO, pubchem="benzotrifluoride"),
    _e("2-methylpentane", "CCCC(C)C", ALKANES),
    _e("cis-decalin", "C1CC[C@H]2CCCC[C@H]2C1", ALKANES, pubchem="493-01-6"),
    _e("cyclohexane", "C1CCCCC1", ALKANES),
    _e("isooctane", "CC(C)CC(C)(C)C", ALKANES, pubchem="2,2,4-trimethylpentane"),
    _e("methylcyclohexane", "CC1CCCCC1", ALKANES),
    _e("methylcyclopentane", "CC1CCCC1", ALKANES),
    _e("heptane", "CCCCCCC", ALKANES),
    _e("hexane", "CCCCCC", ALKANES),
    _e("pentane", "CCCCC", ALKANES),
    _e("2-pentanone", "CCCC(C)=O", KETONES),
    _e("3-pentanone", "CCC(=O)CC", KETONES),
    _e("acetone", "CC(C)=O", KETONES),
    _e("cyclohexanone", "O=C1CCCCC1", KETONES),
    _e("cyclopentanone", "O=C1CCCC1", KETONES),
    _e("methyl ethyl ketone", "CCC(C)=O", KETONES, "MEK", "2-butanone"),
    _e("methyl isobutyl ketone", "CC(=O)CC(C)C", KETONES, "MIBK"),
    _e("TRIS hydrochloride", "[NH3+]C(CO)(CO)CO.[Cl-]", SALTS, "TRIS·HCl",
       "tromethamine hydrochloride"),
    _e("sodium chloride", "[Na+].[Cl-]", SALTS, "NaCl"),
    _e("sodium acetate", "CC(=O)[O-].[Na+]", SALTS, "NaOAc"),
    _e("potassium chloride", "[K+].[Cl-]", SALTS, "KCl"),
    _e("monopotassium phosphate", "OP(=O)(O)[O-].[K+]", SALTS, "KH2PO4",
       "potassium dihydrogen phosphate"),
    _e("lithium chloride", "[Li+].[Cl-]", SALTS, "LiCl"),
    _e("guanidine hydrochloride", "NC(N)=[NH2+].[Cl-]", SALTS, "Gu·HCl"),
    _e("disodium phosphate", "OP(=O)([O-])[O-].[Na+].[Na+]", SALTS, "Na2HPO4",
       "sodium phosphate dibasic"),
    _e("disodium EDTA", "OC(=O)CN(CCN(CC(=O)O)CC(=O)[O-])CC(=O)[O-].[Na+].[Na+]", SALTS,
       pubchem="139-33-3"),
    _e("ammonium sulfate", "[NH4+].[NH4+].[O-]S(=O)(=O)[O-]", SALTS, "(NH4)2SO4"),
)
# fmt: on

GROUPS: tuple[str, ...] = tuple(dict.fromkeys(e.group for e in TABLE))


def sketch_of(entry: Entry) -> Sketch:
    """The entry drawn with Clean's layout (``layout.py``; its ions side by side)."""
    from chemsuite.sketch import layout
    from chemsuite.sketch.edits import History

    sketch = convert.from_text(entry.smiles)
    assert sketch is not None, entry.smiles
    history = History(sketch)
    change = history.change("Clean")
    layout.clean(change)
    history.commit(change)
    return history.sketch


def place(change: Change, entry: Entry, near: Sketch | None = None) -> list[int]:
    """Add *entry* to the right of what is drawn (or at the origin). Returns its atoms."""
    piece = sketch_of(entry)
    sketch = change.sketch if near is None else near
    if sketch.is_empty():
        at = (0.0, 0.0)
    else:
        x0, y0, x1, y1 = sketch.bounds()
        px0, _, px1, _ = piece.bounds()
        at = (x1 + 1.5 + (px1 - px0) / 2, (y0 + y1) / 2)
    return insert(change, convert.placed(piece, at))
