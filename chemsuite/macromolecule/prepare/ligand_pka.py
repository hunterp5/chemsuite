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

"""The microstate model behind a ligand's forms (:mod:`.ligand_sites`): chemsuite pKa.

:class:`LigandPka` is what :func:`.ligand_sites.find` needs: a molecule's protonation
microstates and each one's free energy in water (kT, as :mod:`chemsuite.ionisation.unipka`
gives it), with the model's mean pKa μ, so a form's free energy at a pH is
``G + ln10 · n · (pH − μ)`` kT. :class:`ModelStates` is chemsuite pKa; tests use a stub.
Without one, :mod:`.ligand_sites` uses MolGpKa's site pKa values instead.
"""

from __future__ import annotations

from typing import Protocol

from rdkit import Chem

from chemsuite.ionisation.unipka import Predictor, microstates
from chemsuite.ionisation.unipka.microstates import Microstate

MODEL_NAME = "chemsuite pKa"


class LigandPka(Protocol):
    name: str
    mean: float  # the model's mean pKa, μ

    def microstates(self, mol: Chem.Mol, lower: int, upper: int) -> dict[int, list[Microstate]]:
        """The microstates of *mol* by total charge (``lower`` and ``upper`` as
        :func:`chemsuite.ionisation.unipka.microstates.ensemble` takes them)."""
        ...

    def free_energy(self, mol: Chem.Mol) -> float:
        """The free energy (kT) of the microstate *mol* (heavy atoms, hydrogens counted)."""
        ...


class ModelStates:
    """chemsuite pKa through a :class:`~chemsuite.ionisation.unipka.Predictor`."""

    name = MODEL_NAME

    def __init__(self, predictor: Predictor) -> None:
        self.predictor = predictor
        self.mean = predictor.mean

    def microstates(self, mol: Chem.Mol, lower: int, upper: int) -> dict[int, list[Microstate]]:
        return microstates.ensemble(mol, lower, upper)

    def free_energy(self, mol: Chem.Mol) -> float:
        return self.predictor.free_energy(Chem.MolToSmiles(mol))
