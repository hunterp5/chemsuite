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
"""Colour schemes for protein atoms: RGBA in 0..1, one row per atom.

The default (:data:`ELEMENT`) is what medicinal chemists expect: heteroatoms by element (Jmol),
protein carbons in their chain's soft colour and ligand carbons in a strong green, so the
ligand stands out from the pocket around it.
"""

from __future__ import annotations

import numpy as np

from chemsuite.molview.elements import colours as element_colours
from chemsuite.molview.elements import rgb

from ..structure.model import (
    ADDITIVE,
    AMINO,
    GLYCAN,
    HELIX,
    ION,
    LIGAND,
    NUCLEIC,
    STRAND,
    WATER,
    Structure,
)
from ..structure.sequence import one_letter

ELEMENT = "element"
CHAIN = "chain"
SECONDARY = "secondary"
B_FACTOR = "b_factor"
HYDROPHOBICITY = "hydrophobicity"
KIND = "kind"
SCHEMES = {
    ELEMENT: "Element",
    CHAIN: "Chain",
    SECONDARY: "Secondary Structure",
    B_FACTOR: "B-factor",
    HYDROPHOBICITY: "Hydrophobicity",
    KIND: "Residue Kind",
}

#: soft chain colours for carbons and cartoons (Tableau-like, light enough for O and N to show)
CHAIN_PALETTE = np.array(
    [
        rgb(c)
        for c in (
            0xB0B7C3,
            0x8FB4D9,
            0xE3B58A,
            0xA9CF9A,
            0xD9A0A7,
            0xB9A6D6,
            0xC9B29B,
            0xE9B8D3,
            0xA7D3D0,
            0xD8D59A,
        )
    ],  # fmt: skip
    np.float32,
)
LIGAND_CARBON = np.array(rgb(0x2EC84A), np.float32)
_SS = {HELIX: rgb(0xE0557A), STRAND: rgb(0xF2C14E), 0: rgb(0xC8CCD2)}
_KIND = {
    AMINO: rgb(0xB0B7C3),
    NUCLEIC: rgb(0xE3B58A),
    LIGAND: rgb(0x2EC84A),
    WATER: rgb(0xFF4040),
    ION: rgb(0x9A6AD6),
    ADDITIVE: rgb(0xA0A0A0),
    GLYCAN: rgb(0x5FA8D3),
}
#: Kyte & Doolittle (1982) hydropathy
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip


def _rgba(rgb_rows: np.ndarray) -> np.ndarray:
    out = np.ones((len(rgb_rows), 4), np.float32)
    out[:, :3] = rgb_rows
    return out


def chain_colours(structure: Structure) -> np.ndarray:
    """RGB per chain."""
    return CHAIN_PALETTE[np.arange(structure.n_chains) % len(CHAIN_PALETTE)]


def colour(structure: Structure, scheme: str = ELEMENT) -> np.ndarray:
    """RGBA (n, 4) for every atom."""
    if scheme == CHAIN:
        return _rgba(chain_colours(structure)[structure.atom_chain])
    if scheme == SECONDARY:
        table = np.array([_SS[k] for k in (0, HELIX, STRAND)], np.float32)
        return _rgba(table[structure.ss[structure.residue_of]])
    if scheme == KIND:
        table = np.array([_KIND[k] for k in sorted(_KIND)], np.float32)
        return _rgba(table[structure.atom_kind])
    if scheme == B_FACTOR:
        return _rgba(_ramp(_scaled(structure.b)))
    if scheme == HYDROPHOBICITY:
        values = np.array([KYTE_DOOLITTLE.get(one_letter(str(n)), 0.0) for n in structure.res_name])
        return _rgba(_ramp((values[structure.residue_of] + 4.5) / 9.0, hydrophobic=True))
    out = _rgba(element_colours(structure.element))
    carbon = structure.element == 6
    polymer = np.isin(structure.atom_kind, (AMINO, NUCLEIC))
    out[carbon & polymer, :3] = chain_colours(structure)[structure.atom_chain[carbon & polymer]]
    out[carbon & (structure.atom_kind == LIGAND), :3] = LIGAND_CARBON
    return out


def _scaled(values: np.ndarray) -> np.ndarray:
    low, high = np.percentile(values, [5, 95]) if len(values) else (0.0, 1.0)
    return np.clip((values - low) / max(high - low, 1e-6), 0, 1)


def _ramp(t: np.ndarray, hydrophobic: bool = False) -> np.ndarray:
    """Blue - white - red (B-factor), or teal - white - orange (hydrophobicity)."""
    low, high = (rgb(0x2A9D8F), rgb(0xE07A2E)) if hydrophobic else (rgb(0x3B6FD4), rgb(0xD6453D))
    low, mid, high = np.array(low), np.ones(3), np.array(high)
    t = np.asarray(t, np.float32)[:, None]
    out = np.where(t < 0.5, low + (mid - low) * (t * 2), mid + (high - mid) * (t * 2 - 1))
    return out.astype(np.float32)
