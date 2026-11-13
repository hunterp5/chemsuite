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

"""From microstate free energies to macro pKa values and populations at a pH.

Uni-pKa trains on pKa = μ + [ln Σ_A e^(−G) − ln Σ_B e^(−G)] / ln 10, the sums over the
microstates of the protonated macrostate A and the deprotonated one B, G in kT and μ the
training set's mean pKa. That fixes every microstate's free energy at a pH up to one constant:

    E_i(pH) = G_i + ln 10 · n_i · (pH − μ)

with n_i its protons (any common offset cancels). Populations are the Boltzmann weights of E.
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

LN10 = float(np.log(10.0))


def energies(g: np.ndarray, protons: np.ndarray, ph: float, mean: float) -> np.ndarray:
    """Each microstate's free energy (kT) at *ph*."""
    return np.asarray(g, dtype=np.float64) + LN10 * np.asarray(protons) * (ph - mean)


def populations(g: np.ndarray, protons: np.ndarray, ph: float, mean: float) -> np.ndarray:
    """Each microstate's share at *ph* (they sum to 1)."""
    e = -energies(g, protons, ph, mean)
    return np.exp(e - logsumexp(e))


def macro_pka(g_acid: np.ndarray, g_base: np.ndarray, mean: float) -> float:
    """The pKa between the macrostate *g_acid* (one proton more) and *g_base*."""
    return mean + float(logsumexp(-np.asarray(g_acid)) - logsumexp(-np.asarray(g_base))) / LN10
