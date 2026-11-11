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
"""What each residue is: polymer, ligand, water, ion, additive or glycan.

The file's entities decide first (polymer, water, branched sugars); a residue outside any
entity (a PDB file without SEQRES) is judged by its name. Among the rest, a lone metal or
halide atom is an ion, a buffer or crystallisation additive (gemmi's list plus the common ones
below) is an additive, and everything else is a ligand: the residues a medicinal chemist cares
about.
"""

from __future__ import annotations

from functools import cache

import gemmi

from .model import ADDITIVE, AMINO, GLYCAN, ION, LIGAND, NUCLEIC, WATER

#: additives gemmi does not list as buffers
ADDITIVES = frozenset(
    [
        "PO4",
        "BME",
        "DMS",
        "MPD",
        "PG4",
        "PGE",
        "1PE",
        "P6G",
        "PE4",
        "IMD",
        "TRS",
        "MES",
        "EPE",
        "FMT",
        "ACT",
        "CIT",
        "TLA",
        "PEG",
        "IPA",
        "EOH",
        "MRD",
        "BU3",
        "NO3",
        "SCN",
        "IOD",
        "AZI",
        "CO3",
        "MLI",
        "MLA",
        "SIN",
    ]
)
#: halides and metals that, alone in a residue, are ions
_ION_ELEMENTS = frozenset(
    [
        "F",
        "CL",
        "BR",
        "I",
        "LI",
        "NA",
        "K",
        "RB",
        "CS",
        "MG",
        "CA",
        "SR",
        "BA",
        "MN",
        "FE",
        "CO",
        "NI",
        "CU",
        "ZN",
        "CD",
        "HG",
        "PT",
        "AU",
        "AG",
        "PB",
        "AL",
        "GA",
        "TL",
    ]
)
_WATERS = frozenset({"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"})
_AMINO_KINDS = {gemmi.ResidueKind.AA, gemmi.ResidueKind.AAD, gemmi.ResidueKind.PAA,
                gemmi.ResidueKind.MAA}  # fmt: skip
_NUCLEIC_KINDS = {gemmi.ResidueKind.RNA, gemmi.ResidueKind.DNA}
_SUGAR_KINDS = {gemmi.ResidueKind.PYR, gemmi.ResidueKind.KET}


@cache
def _tabulated(name: str) -> gemmi.ResidueKind | None:
    info = gemmi.find_tabulated_residue(name)
    return info.kind if info is not None else None


def residue_kind(name: str, entity_type: gemmi.EntityType, heavy_elements: tuple[str, ...]) -> int:
    """The kind of one residue; ``heavy_elements`` are its non-hydrogen atoms' symbols
    (upper case)."""
    tabulated = _tabulated(name)
    if entity_type == gemmi.EntityType.Polymer:
        return NUCLEIC if tabulated in _NUCLEIC_KINDS else AMINO
    if entity_type == gemmi.EntityType.Water or name in _WATERS:
        return WATER
    if entity_type == gemmi.EntityType.Branched:
        return GLYCAN
    if entity_type == gemmi.EntityType.Unknown:
        if tabulated in _AMINO_KINDS:
            return AMINO
        if tabulated in _NUCLEIC_KINDS:
            return NUCLEIC
    return non_polymer_kind(name, heavy_elements)


def non_polymer_kind(name: str, heavy_elements: tuple[str, ...]) -> int:
    if len(heavy_elements) == 1 and heavy_elements[0] in _ION_ELEMENTS:
        return ION
    if name in ADDITIVES or _tabulated(name) == gemmi.ResidueKind.BUF:
        return ADDITIVE
    if _tabulated(name) in _SUGAR_KINDS:
        return GLYCAN
    return LIGAND
