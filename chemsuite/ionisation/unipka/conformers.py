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

"""A microstate's atoms and conformers, made as Uni-pKa's training data were.

Upstream (``scripts/preprocess_pka.py``) reads the SMILES, adds hydrogens, and embeds the same
molecule once per seed 0, 1, … with RDKit's default embedding, then MMFF; a failed embedding
falls back to the 2D layout. A last conformer is the 2D layout itself. The model averages its
free energy over them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

N_3D = 10  # upstream's count (with the 2D layout, 11)


@dataclass(frozen=True)
class Atoms:
    symbols: tuple[str, ...]
    charges: tuple[int, ...]
    coords: np.ndarray  # (conformers, atoms, 3)


def _layout(smiles: str) -> np.ndarray:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.Compute2DCoords(mol)
    return mol.GetConformer().GetPositions()


def atoms(smiles: str, n_3d: int = N_3D) -> Atoms:
    """The atoms of the microstate *smiles*, hydrogens explicit, in *n_3d* 3D conformers and
    the 2D layout."""
    RDLogger.DisableLog("rdApp.*")
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        coords = []
        for seed in range(n_3d):
            if AllChem.EmbedMolecule(mol, randomSeed=seed) == 0:
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                    coords.append(mol.GetConformer().GetPositions())
                except (ValueError, RuntimeError):  # no MMFF parameters: upstream's 2D layout
                    coords.append(_layout(smiles))
            else:
                coords.append(_embedded_without_hydrogens(smiles, seed))
        coords.append(_layout(smiles))
    finally:
        RDLogger.EnableLog("rdApp.*")
    symbols = tuple(a.GetSymbol() for a in mol.GetAtoms())
    charges = tuple(a.GetFormalCharge() for a in mol.GetAtoms())
    return Atoms(symbols, charges, np.array(coords, dtype=np.float32))


def _embedded_without_hydrogens(smiles: str, seed: int) -> np.ndarray:
    """Upstream's second try: embed the heavy atoms, then place hydrogens; else the layout."""
    heavy = Chem.MolFromSmiles(smiles)
    if AllChem.EmbedMolecule(heavy, maxAttempts=5000, randomSeed=seed) != 0:
        return _layout(smiles)
    mol = Chem.AddHs(heavy, addCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except (ValueError, RuntimeError):
        return _layout(smiles)
    return mol.GetConformer().GetPositions()
