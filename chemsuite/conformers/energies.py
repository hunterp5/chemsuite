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
"""Each conformer's energy, kept on the conformer itself.

A generated ensemble's conformers carry their force-field energy (kcal/mol), their energy above
the ensemble's lowest (ΔE) and their Boltzmann population at 298.15 K over the kept conformers.
ΔE and population are fixed when the ensemble is made, so a conformer kept alone later (a
superposed pose) still says where it stood in its ensemble. They are RDKit conformer properties
and travel with the conformer through copies, hydrogens and alignment; the table's blob keeps
them (:func:`chemsuite.core.mols.mol_to_blob`), a molfile does not.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

ENERGY = "Energy"
RELATIVE = "Relative Energy"
POPULATION = "Population"
FORCE_FIELD = "Force Field"
PROPS = (ENERGY, RELATIVE, POPULATION, FORCE_FIELD)

#: the temperature populations are worked out at (K)
TEMPERATURE = 298.15
#: the gas constant in kcal/(mol·K)
GAS_CONSTANT = 0.0019872


@dataclass(frozen=True)
class Row:
    """One conformer's energies; *number* counts from 1 in stored order."""

    number: int
    energy: float
    delta: float
    population: float  # a fraction: 0–1


def boltzmann(deltas: Sequence[float], temperature: float = TEMPERATURE) -> list[float]:
    """Boltzmann populations (fractions summing to 1) for energies in kcal/mol."""
    if not len(deltas):
        return []
    values = np.asarray(deltas, dtype=float)
    weights = np.exp(-(values - values.min()) / (GAS_CONSTANT * temperature))
    return (weights / weights.sum()).tolist()


def annotate(mol: Chem.Mol, energies: Sequence[float], force_field: str) -> None:
    """Store *energies* (kcal/mol, one per conformer in stored order) on *mol*'s conformers.
    Nothing is stored without a force field, or when any energy is missing."""
    conformers = list(mol.GetConformers())
    if not force_field or len(energies) != len(conformers) or not conformers:
        return
    if not all(math.isfinite(e) for e in energies):
        return
    lowest = min(energies)
    deltas = [e - lowest for e in energies]
    for conformer, energy, delta, population in zip(
        conformers, energies, deltas, boltzmann(deltas), strict=True
    ):
        conformer.SetDoubleProp(ENERGY, float(energy))
        conformer.SetDoubleProp(RELATIVE, float(delta))
        conformer.SetDoubleProp(POPULATION, float(population))
        conformer.SetProp(FORCE_FIELD, force_field)


def clear(mol: Chem.Mol) -> None:
    """Remove the stored energies from every conformer of *mol*."""
    for conformer in mol.GetConformers():
        for name in PROPS:
            if conformer.HasProp(name):
                conformer.ClearProp(name)


def has_energies(conformer: Chem.Conformer) -> bool:
    return conformer.HasProp(RELATIVE)


def row_of(conformer: Chem.Conformer, number: int) -> Row | None:
    """*conformer*'s energies (None when it has none)."""
    if not has_energies(conformer):
        return None
    return Row(
        number,
        conformer.GetDoubleProp(ENERGY),
        conformer.GetDoubleProp(RELATIVE),
        conformer.GetDoubleProp(POPULATION),
    )


def rows(mol: Chem.Mol | None) -> list[Row]:
    """Every conformer's energies in stored order; empty unless all of them carry energies."""
    if mol is None or not mol.GetNumConformers():
        return []
    out = [row_of(c, n) for n, c in enumerate(mol.GetConformers(), 1)]
    return [r for r in out if r is not None] if all(r is not None for r in out) else []


def force_field(mol: Chem.Mol | None) -> str:
    """The force field the energies come from ("" when there are none)."""
    if mol is None or not mol.GetNumConformers():
        return ""
    first = mol.GetConformers()[0]
    return first.GetProp(FORCE_FIELD) if first.HasProp(FORCE_FIELD) else ""


def delta_of(mol: Chem.Mol | None, position: int = 0) -> float:
    """ΔE (kcal/mol) of the conformer at *position*; NaN when it has none."""
    if mol is None or not 0 <= position < mol.GetNumConformers():
        return math.nan
    conformer = mol.GetConformers()[position]
    return conformer.GetDoubleProp(RELATIVE) if has_energies(conformer) else math.nan
