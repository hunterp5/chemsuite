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

"""chemsuite pKa: the free energy of every protonation microstate, so populations at any pH.

chemsuite's own weights (trained on pKaHub's data, docs/UNIPKA_PLAN.md) on the Uni-pKa
method (Luo et al., JACS Au 2024; Apache-2.0) and the Uni-Mol network (Zhou et al., ICLR 2023;
MIT). It is not the published Uni-pKa model, whose weights are non-commercial.
:mod:`.microstates` lists the microstates, :mod:`.conformers` builds each one's atoms,
:mod:`.model` gives its free energy, and :mod:`.thermo` turns free energies into macro pKa
values and populations. Numpy only.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from . import conformers, microstates, model, thermo
from .microstates import Microstate


@dataclass(frozen=True)
class Ensemble:
    """A molecule's microstates with their free energies (kT), and the model's mean pKa."""

    states: tuple[Microstate, ...]
    g: np.ndarray
    mean: float

    @property
    def protons(self) -> np.ndarray:
        return np.array([sum(m.hydrogens) for m in self.states])

    def populations(self, ph: float) -> np.ndarray:
        return thermo.populations(self.g, self.protons, ph, self.mean)

    def macro_pkas(self) -> list[tuple[int, float]]:
        """(charge of the protonated macrostate, pKa) for each adjacent pair, lowest first."""
        charges = np.array([m.charge for m in self.states])
        out = []
        for q in sorted(set(charges.tolist()))[1:]:
            if (charges == q - 1).any():
                out.append((q, thermo.macro_pka(self.g[charges == q], self.g[charges == q - 1],
                                                self.mean)))  # fmt: skip
        return sorted(out, key=lambda x: x[1])


class Predictor:
    """The folds of one weights archive, with free energies cached by microstate."""

    def __init__(self, path: Path, n_3d: int = conformers.N_3D) -> None:
        with np.load(path) as archive:
            self.mean = float(archive["mean"])
            folds = sorted({k.split("/", 1)[0] for k in archive.files if "/" in k})
            self.fold_means = (archive["fold_means"] if "fold_means" in archive.files
                               else np.full(len(folds), self.mean))  # fmt: skip
        self.nets = [model.Network.load(path, int(f[1:])) for f in folds]
        self.n_3d = n_3d
        self._cache: dict[str, float] = {}

    def free_energy(self, smiles: str) -> float:
        """The microstate's free energy (kT), averaged over conformers and folds."""
        if smiles not in self._cache:
            found = conformers.atoms(smiles, self.n_3d)
            g = [model.free_energy(net, list(found.symbols), list(found.charges), found.coords)
                 for net in self.nets]  # fmt: skip
            g = [float(np.mean(x)) for x in g]  # each fold over its conformers
            charge = int(sum(found.charges))
            self._cache[smiles] = thermo.common_mean(g, charge, self.fold_means, self.mean)
        return self._cache[smiles]

    def ensemble(
        self, mol: Chem.Mol, lower: int = microstates.LOWER, upper: int = microstates.UPPER
    ) -> Ensemble:
        """Every microstate of *mol* from charge *lower* to *upper*, with free energies."""
        states = tuple(m for pool in microstates.ensemble(mol, lower, upper).values() for m in pool)
        g = np.array([self.free_energy(m.smiles) for m in states])
        return Ensemble(states, g, self.mean)
