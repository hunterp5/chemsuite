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
"""Every descriptor the Descriptors tool offers: key, column title, group, format and function.

Keys are stable identifiers (sent to worker processes and remembered between runs); titles are
the column names the results get. Functions take a sanitized RDKit molecule.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, inchi, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from . import medchem

TEXT = None  # ``decimals`` for descriptors whose value is text


@dataclass(frozen=True)
class Descriptor:
    key: str
    title: str
    group: str
    fn: Callable[[Chem.Mol], object]
    decimals: int | None = 2  # display decimals for numbers; TEXT for text values
    help: str = ""

    @property
    def is_text(self) -> bool:
        return self.decimals is TEXT


ESTIMATE = (
    " Estimated without a pKa model (basic pKa guessed from the nitrogen count); "
    "use for ranking, not as a measured value."
)


def _count(atomic_number: int) -> Callable[[Chem.Mol], int]:
    return lambda mol: sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_number)


def _murcko(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))


def _shape(fn: Callable[..., float]) -> Callable[[Chem.Mol], float]:
    """A 3D descriptor of the first (lowest-energy) conformer; NaN without 3D coordinates."""

    def compute(mol: Chem.Mol) -> float:
        if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
            raise ValueError("no 3D conformer")
        return fn(mol)

    compute.__name__ = fn.__name__
    return compute


SHAPE_GROUP = "3D Shape"
_SHAPE_HELP = " Needs 3D structures (Chemistry ▸ Generate Conformers); uses the first conformer."

_SHAPES = (
    ("PMI1", "PMI 1", rdMolDescriptors.CalcPMI1, 1, "Smallest principal moment of inertia."),
    ("PMI2", "PMI 2", rdMolDescriptors.CalcPMI2, 1, "Middle principal moment of inertia."),
    ("PMI3", "PMI 3", rdMolDescriptors.CalcPMI3, 1, "Largest principal moment of inertia."),
    (
        "NPR1",
        "NPR1",
        rdMolDescriptors.CalcNPR1,
        3,
        "Normalised PMI ratio I1/I3: 0 rod, 0.5 disc, 1 sphere (x of the PMI triangle).",
    ),
    (
        "NPR2",
        "NPR2",
        rdMolDescriptors.CalcNPR2,
        3,
        "Normalised PMI ratio I2/I3 (y of the PMI triangle).",
    ),
    ("PBF", "PBF", rdMolDescriptors.CalcPBF, 3, "Plane of best fit: mean distance from it (Å)."),
    (
        "RadiusOfGyration",
        "Radius of Gyration",
        rdMolDescriptors.CalcRadiusOfGyration,
        2,
        "Size of the 3D structure (Å).",
    ),
    ("Asphericity", "Asphericity", rdMolDescriptors.CalcAsphericity, 3, "0 for a sphere."),
    ("Eccentricity", "Eccentricity", rdMolDescriptors.CalcEccentricity, 3, "0 for a sphere."),
    (
        "SpherocityIndex",
        "Spherocity Index",
        rdMolDescriptors.CalcSpherocityIndex,
        3,
        "1 for a sphere, 0 for a flat molecule.",
    ),
    (
        "InertialShapeFactor",
        "Inertial Shape Factor",
        rdMolDescriptors.CalcInertialShapeFactor,
        5,
        "I2 / (I1 · I3).",
    ),
)


def _tpsa_sp(mol: Chem.Mol) -> float:
    return rdMolDescriptors.CalcTPSA(mol, includeSandP=True)


def _generic_murcko(mol: Chem.Mol) -> str:
    core = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(core)) if core.GetNumAtoms() else ""


_ITEMS: tuple[Descriptor, ...] = (
    # -- physicochemical --
    Descriptor("MolWt", "MW", "Physicochemical", Descriptors.MolWt, 2, "Average molecular weight."),
    Descriptor(
        "ExactMolWt",
        "Exact Mass",
        "Physicochemical",
        Descriptors.ExactMolWt,
        4,
        "Monoisotopic mass.",
    ),
    Descriptor("MolLogP", "cLogP", "Physicochemical", Crippen.MolLogP, 2, "Wildman-Crippen logP."),
    Descriptor(
        "CLOGD74",
        "cLogD 7.4 (est.)",
        "Physicochemical",
        medchem.clogd74,
        2,
        "cLogP corrected for ionisation at pH 7.4." + ESTIMATE,
    ),
    Descriptor(
        "ESOL",
        "LogS (ESOL)",
        "Physicochemical",
        medchem.esol_logs,
        2,
        "Delaney ESOL intrinsic aqueous solubility, log10(mol/L).",
    ),
    Descriptor(
        "TPSA",
        "TPSA",
        "Physicochemical",
        rdMolDescriptors.CalcTPSA,
        1,
        "Topological polar surface area (N, O).",
    ),
    Descriptor(
        "TPSA_SP",
        "TPSA (N,O,S,P)",
        "Physicochemical",
        _tpsa_sp,
        1,
        "Polar surface area counting polar S and P too, as SwissADME and the BOILED-Egg do.",
    ),
    Descriptor("MolMR", "Molar Refractivity", "Physicochemical", Crippen.MolMR, 2),
    Descriptor(
        "LabuteASA",
        "Labute ASA",
        "Physicochemical",
        rdMolDescriptors.CalcLabuteASA,
        1,
        "Labute approximate surface area.",
    ),
    Descriptor(
        "FractionCSP3",
        "Fsp3",
        "Physicochemical",
        rdMolDescriptors.CalcFractionCSP3,
        2,
        "Fraction of carbons that are sp3.",
    ),
    # -- drug-likeness --
    Descriptor(
        "QED",
        "QED",
        "Drug-likeness",
        QED.qed,
        3,
        "Quantitative estimate of drug-likeness (Bickerton 2012), 0-1.",
    ),
    Descriptor(
        "RO5_VIOLATIONS",
        "Ro5 Violations",
        "Drug-likeness",
        medchem.lipinski_violations,
        0,
        "Lipinski rule-of-5 violations (0-4).",
    ),
    Descriptor(
        "RO5_PASS",
        "Ro5 Pass",
        "Drug-likeness",
        medchem.ro5_pass,
        TEXT,
        "Yes when there are no rule-of-5 violations.",
    ),
    Descriptor(
        "CNS_MPO",
        "CNS MPO (est.)",
        "Drug-likeness",
        medchem.cns_mpo,
        2,
        "Wager 2010 CNS MPO score, 0-6 (≥ 4 favourable)." + ESTIMATE,
    ),
    Descriptor(
        "BBB_SCORE",
        "BBB Score (est.)",
        "Drug-likeness",
        medchem.bbb_score,
        2,
        "Gupta 2019 BBB score, 0-6 (≥ 4 CNS-like)." + ESTIMATE,
    ),
    Descriptor(
        "AB_MPS",
        "AB-MPS (est.)",
        "Drug-likeness",
        medchem.ab_mps,
        2,
        "AbbVie multiparameter score for beyond-rule-of-5 space (≤ 14 favourable)." + ESTIMATE,
    ),
    # -- counts --
    Descriptor("NumHDonors", "HBD", "Counts", Lipinski.NumHDonors, 0, "Lipinski H-bond donors."),
    Descriptor(
        "NumHAcceptors", "HBA", "Counts", Lipinski.NumHAcceptors, 0, "Lipinski H-bond acceptors."
    ),
    Descriptor("NHOHCount", "NH/OH Count", "Counts", Lipinski.NHOHCount, 0),
    Descriptor("NOCount", "N+O Count", "Counts", Lipinski.NOCount, 0),
    Descriptor("NumRotatableBonds", "Rotatable Bonds", "Counts", Lipinski.NumRotatableBonds, 0),
    Descriptor("HeavyAtomCount", "Heavy Atoms", "Counts", Lipinski.HeavyAtomCount, 0),
    Descriptor("NumHeteroatoms", "Heteroatoms", "Counts", Lipinski.NumHeteroatoms, 0),
    Descriptor(
        "NumValenceElectrons", "Valence Electrons", "Counts", Descriptors.NumValenceElectrons, 0
    ),
    Descriptor(
        "NET_CHARGE", "Net Charge", "Counts", medchem.net_formal_charge, 0, "Sum of formal charges."
    ),
    # -- rings --
    Descriptor("RingCount", "Rings", "Rings", rdMolDescriptors.CalcNumRings, 0),
    Descriptor(
        "NumAromaticRings", "Aromatic Rings", "Rings", rdMolDescriptors.CalcNumAromaticRings, 0
    ),
    Descriptor(
        "NumAliphaticRings", "Aliphatic Rings", "Rings", rdMolDescriptors.CalcNumAliphaticRings, 0
    ),
    Descriptor(
        "NumSaturatedRings", "Saturated Rings", "Rings", rdMolDescriptors.CalcNumSaturatedRings, 0
    ),
    Descriptor("NumSpiroAtoms", "Spiro Atoms", "Rings", rdMolDescriptors.CalcNumSpiroAtoms, 0),
    Descriptor(
        "NumBridgeheadAtoms",
        "Bridgehead Atoms",
        "Rings",
        rdMolDescriptors.CalcNumBridgeheadAtoms,
        0,
    ),
    # -- atoms --
    *(
        Descriptor(f"Count_{symbol}", f"{name} Atoms", "Atom Counts", _count(number), 0)
        for symbol, name, number in (
            ("C", "C", 6),
            ("N", "N", 7),
            ("O", "O", 8),
            ("S", "S", 16),
            ("P", "P", 15),
            ("F", "F", 9),
            ("Cl", "Cl", 17),
            ("Br", "Br", 35),
            ("I", "I", 53),
        )
    ),
    # -- complexity --
    Descriptor("BertzCT", "Bertz Complexity", "Complexity", Descriptors.BertzCT, 1),
    Descriptor("BalabanJ", "Balaban J", "Complexity", Descriptors.BalabanJ, 3),
    Descriptor(
        "HallKierAlpha", "Hall-Kier Alpha", "Complexity", rdMolDescriptors.CalcHallKierAlpha, 2
    ),
    Descriptor("Kappa1", "Kappa 1", "Complexity", rdMolDescriptors.CalcKappa1, 2),
    Descriptor("Kappa2", "Kappa 2", "Complexity", rdMolDescriptors.CalcKappa2, 2),
    Descriptor("Kappa3", "Kappa 3", "Complexity", rdMolDescriptors.CalcKappa3, 2),
    # -- charges --
    Descriptor(
        "MaxPartialCharge",
        "Max Partial Charge",
        "Partial Charges",
        Descriptors.MaxPartialCharge,
        3,
        "Gasteiger charges.",
    ),
    Descriptor(
        "MinPartialCharge",
        "Min Partial Charge",
        "Partial Charges",
        Descriptors.MinPartialCharge,
        3,
        "Gasteiger charges.",
    ),
    Descriptor(
        "MaxAbsPartialCharge",
        "Max Abs Partial Charge",
        "Partial Charges",
        Descriptors.MaxAbsPartialCharge,
        3,
        "Gasteiger charges.",
    ),
    Descriptor(
        "MinAbsPartialCharge",
        "Min Abs Partial Charge",
        "Partial Charges",
        Descriptors.MinAbsPartialCharge,
        3,
        "Gasteiger charges.",
    ),
    # -- identifiers --
    *(
        Descriptor(key, title, SHAPE_GROUP, _shape(fn), decimals, help + _SHAPE_HELP)
        for key, title, fn, decimals, help in _SHAPES
    ),
    Descriptor("SMILES", "Canonical SMILES", "Identifiers", Chem.MolToSmiles, TEXT),
    Descriptor("INCHIKEY", "InChIKey", "Identifiers", inchi.MolToInchiKey, TEXT),
    Descriptor("INCHI", "InChI", "Identifiers", inchi.MolToInchi, TEXT),
    Descriptor("FORMULA", "Formula", "Identifiers", rdMolDescriptors.CalcMolFormula, TEXT),
    Descriptor(
        "MURCKO",
        "Murcko Scaffold",
        "Identifiers",
        _murcko,
        TEXT,
        "Bemis-Murcko framework (rings and linkers) as SMILES.",
    ),
    Descriptor(
        "MURCKO_GENERIC",
        "Generic Scaffold",
        "Identifiers",
        _generic_murcko,
        TEXT,
        "Murcko framework with every atom carbon and every bond single.",
    ),
)

CATALOG: dict[str, Descriptor] = {d.key: d for d in _ITEMS}
GROUPS: tuple[str, ...] = tuple(dict.fromkeys(d.group for d in _ITEMS))
COMMON: tuple[str, ...] = (
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "HeavyAtomCount",
    "NumAromaticRings",
    "FractionCSP3",
    "QED",
)


def in_group(group: str) -> list[Descriptor]:
    return [d for d in _ITEMS if d.group == group]
