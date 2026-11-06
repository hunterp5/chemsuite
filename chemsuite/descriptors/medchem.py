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
"""Medicinal-chemistry scores used by the descriptor catalog (ported from MolManager).

* **CNS MPO**: Wager et al., ACS Chem. Neurosci. 2010 (doi:10.1021/cn100008c), Table 1.
* **BBB Score**: Gupta et al., J. Med. Chem. 2019 (doi:10.1021/acs.jmedchem.9b01220), Tables 1, 4.
* **AB-MPS**: DeGoey et al., J. Med. Chem. 2018 (doi:10.1021/acs.jmedchem.7b00717).
* **ESOL intrinsic log S**: Delaney, J. Chem. Inf. Comput. Sci. 2004 (doi:10.1021/ci034243x).

cLogD 7.4 and the scores that need a basic pKa (CNS MPO, BBB Score, AB-MPS) use MolGpKa's
predicted pKa values (``chemsuite.ionisation``). A molecule without a basic site is scored as
if its most basic pKa were 2, as MolManager did.
"""

from __future__ import annotations

import math

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from chemsuite.ionisation.profile import clogd, profile

_AROMATIC_ATOM = Chem.MolFromSmarts("a")


def _mono_dec(x: float, x_good: float, x_bad: float) -> float:
    """Decreasing desirability: 1 for x <= x_good, 0 for x >= x_bad, linear between."""
    if x <= x_good:
        return 1.0
    if x >= x_bad:
        return 0.0
    return 1.0 - (x - x_good) / (x_bad - x_good)


def _tpsa_hump(tpsa: float) -> float:
    """TPSA desirability from Wager 2010 Table 1 (piecewise linear)."""
    if tpsa <= 20.0:
        return 0.0
    if tpsa <= 40.0:
        return (tpsa - 20.0) / 20.0
    if tpsa <= 90.0:
        return 1.0
    if tpsa <= 120.0:
        return 1.0 - (tpsa - 90.0) / 30.0
    return 0.0


def lipinski_violations(mol: Chem.Mol) -> int:
    """Rule-of-5 violations (0-4) with RDKit's Lipinski H-bond definitions."""
    return (
        int(Descriptors.MolWt(mol) > 500.0)
        + int(Descriptors.MolLogP(mol) > 5.0)
        + int(Lipinski.NumHDonors(mol) > 5)
        + int(Lipinski.NumHAcceptors(mol) > 10)
    )


def ro5_pass(mol: Chem.Mol) -> str:
    return "Yes" if lipinski_violations(mol) == 0 else "No"


NO_BASE_PKA = 2.0  # the most basic pKa used for a molecule without a basic site


def basic_pka(mol: Chem.Mol) -> float:
    """MolGpKa's pKa of the most basic site (its conjugate acid); NO_BASE_PKA without one."""
    found = profile(mol)
    strongest = found.most_basic if found is not None else None
    return NO_BASE_PKA if strongest is None else strongest


def clogd74(mol: Chem.Mol) -> float:
    """cLogD at pH 7.4: Crippen cLogP corrected with MolGpKa's pKa values."""
    return clogd(mol, 7.4)


def esol_logs(mol: Chem.Mol) -> float:
    """ESOL intrinsic log10(S / mol L⁻¹) with Delaney's original coefficients."""
    heavy = mol.GetNumAtoms()
    aromatic = len(mol.GetSubstructMatches(_AROMATIC_ATOM)) / heavy if heavy else 0.0
    return (
        0.16
        - 0.63 * float(Crippen.MolLogP(mol))
        - 0.0062 * float(Descriptors.MolWt(mol))
        + 0.066 * float(Lipinski.NumRotatableBonds(mol))
        - 0.74 * aromatic
    )


def ab_mps(mol: Chem.Mol) -> float:
    """AbbVie MPS: |cLogD7.4 − 3| + aromatic rings + rotatable bonds (≤ 14 favourable)."""
    return (
        abs(clogd74(mol) - 3.0)
        + float(rdMolDescriptors.CalcNumAromaticRings(mol))
        + float(Lipinski.NumRotatableBonds(mol))
    )


def cns_mpo(mol: Chem.Mol) -> float:
    """CNS MPO (0-6) from Wager 2010 Table 1 desirability functions."""
    clogp = float(Crippen.MolLogP(mol))
    pka = basic_pka(mol)
    return (
        _mono_dec(clogp, 3.0, 5.0)
        + _mono_dec(clogd74(mol), 2.0, 4.0)
        + _mono_dec(float(Descriptors.MolWt(mol)), 360.0, 500.0)
        + _tpsa_hump(float(Descriptors.TPSA(mol)))
        + _mono_dec(float(Lipinski.NumHDonors(mol)), 0.5, 3.5)
        + _mono_dec(pka, 8.0, 10.0)
    )


# Gupta 2019 Table 1: the aromatic-ring term is stepwise; the rest are polynomials normalised by
# their local maximum on (low, high].
_BBB_AROMATIC = {0: 0.336376, 1: 0.816016, 2: 1.0, 3: 0.691115, 4: 0.199399}


def _bbb_piece(x: float, low: float, high: float, norm: float, coeffs: tuple[float, ...]) -> float:
    if x <= low or x > high:
        return 0.0
    acc = 0.0
    for coeff in coeffs:  # highest degree first
        acc = acc * x + coeff
    return min(1.0, max(0.0, acc / norm))


def bbb_score_from_descriptors(
    aromatic_rings: float,
    heavy_atoms: float,
    mol_wt: float,
    hbd: float,
    hba: float,
    tpsa: float,
    pka: float,
) -> float:
    """``P(AroR) + P(HA) + 1.5 P(MWHBN) + 2 P(TPSA) + 0.5 P(pKa)``; 4-6 is CNS-like."""
    mwhbn = (hbd + hba) / math.sqrt(mol_wt) if mol_wt > 0.0 else 0.0
    return (
        _BBB_AROMATIC.get(int(aromatic_rings), 0.0)
        + _bbb_piece(heavy_atoms, 5.0, 45.0, 0.624231, (0.0000443, -0.004556, 0.12775, -0.463))
        + 1.5 * _bbb_piece(mwhbn, 0.05, 0.45, 0.72258, (26.733, -31.495, 9.5202, -0.1358))
        + 2.0 * _bbb_piece(tpsa, 0.0, 120.0, 0.9598, (-0.0067, 0.9598))
        + 0.5
        * _bbb_piece(pka, 3.0, 11.0, 0.597488, (0.00045068, -0.016331, 0.18618, -0.71043, 0.8579))
    )


def bbb_score(mol: Chem.Mol) -> float:
    return bbb_score_from_descriptors(
        float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        float(mol.GetNumHeavyAtoms()),
        float(Descriptors.MolWt(mol)),
        float(Lipinski.NumHDonors(mol)),
        float(Lipinski.NumHAcceptors(mol)),
        float(Descriptors.TPSA(mol)),
        basic_pka(mol),
    )


def net_formal_charge(mol: Chem.Mol) -> int:
    return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
