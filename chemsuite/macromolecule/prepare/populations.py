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

"""How much of each site's protonated form there is at each pH, and each site's pKa.

The energies are P2's tables (:mod:`.energy`, built with titration). Every pH of the grid is
solved at once, by one of:

- :func:`mean_field`: each site's state populations in the average field of the others,
  iterated to self-consistency. Exact for a lone site; fast enough to fit weights with.
- :func:`monte_carlo`: Metropolis sampling, one replica per pH (as in MCCE, Song et al.,
  *J. Comput. Chem.* 30, 2231–2247, 2009), with moves of single sites and of coupled pairs.

A site's pKa is the pH at which its protonated population crosses ½ (interpolated), held to
the ends of the grid when it does not cross.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from .energy import RT_LN10, Tables, Weights

KT = RT_LN10 / np.log(10)  # kcal/mol at 298 K
GRID = np.arange(0.0, 14.0 + 1e-9, 0.25)
MEAN_FIELD, MONTE_CARLO = "mean field", "monte carlo"
COUPLED = 2.0  # kcal/mol: pairs coupled this strongly also move together in Monte Carlo


@dataclass(frozen=True)
class Titration:
    grid: np.ndarray  # (P,)
    populations: list[np.ndarray]  # per site: (P, S)
    protonated: list[np.ndarray]  # per site: (P,), the share of its most protonated form
    pka: list[float]  # per site; NaN for a site that does not titrate

    def at(self, ph: float) -> list[float]:
        """Each site's protonated share at *ph*."""
        return [float(np.interp(ph, self.grid, curve)) for curve in self.protonated]


def titration(
    table: Tables,
    weights: Weights,
    grid: np.ndarray = GRID,
    method: str = MEAN_FIELD,
    cancel: CancelToken | None = None,
    seed: int = 0,
) -> Titration:
    single, pairs = table.energies(weights, 0.0)
    if method == MEAN_FIELD:
        found = mean_field(single, pairs, table.protons, grid, cancel)
    else:
        found = monte_carlo(single, pairs, table.protons, grid, cancel, seed=seed)
    protonated, pkas = [], []
    for pop, n in zip(found, table.protons, strict=True):
        if len(n) and n.max() > n.min():
            curve = pop[:, n == n.max()].sum(axis=1)
            protonated.append(curve)
            pkas.append(pka(grid, curve))
        else:
            protonated.append(np.zeros(len(grid)))
            pkas.append(float("nan"))
    return Titration(grid, found, protonated, pkas)


def pka(grid: np.ndarray, curve: np.ndarray) -> float:
    """Where *curve* (protonated share, falling with pH) first drops below ½."""
    below = np.flatnonzero(curve < 0.5)
    if not len(below):
        return float(grid[-1])
    k = int(below[0])
    if k == 0:
        return float(grid[0])
    a, b = curve[k - 1], curve[k]
    return float(grid[k - 1] + (a - 0.5) / (a - b) * (grid[k] - grid[k - 1]))


def _at_ph(single: list[np.ndarray], protons: list[np.ndarray], grid: np.ndarray) -> list:
    """Each site's own energies at every pH: (P, S)."""
    return [e[None, :] + RT_LN10 * grid[:, None] * n[None, :]
            for e, n in zip(single, protons, strict=True)]  # fmt: skip


def _links(n: int, pairs: dict) -> list[list[tuple[int, np.ndarray]]]:
    """Each site's neighbours with the pair energies seen from it: (S_i, S_j)."""
    links: list[list[tuple[int, np.ndarray]]] = [[] for _ in range(n)]
    for (i, j), m in pairs.items():
        links[i].append((j, m))
        links[j].append((i, m.T))
    return links


def _boltzmann(e: np.ndarray) -> np.ndarray:
    """Populations from energies (P, S), each row summing to 1."""
    w = np.exp(-(e - e.min(axis=1, keepdims=True)) / KT)
    return w / w.sum(axis=1, keepdims=True)


# -- mean field ------------------------------------------------------------------------------------
def mean_field(
    single: list[np.ndarray],
    pairs: dict,
    protons: list[np.ndarray],
    grid: np.ndarray = GRID,
    cancel: CancelToken | None = None,
    rounds: int = 200,
    tolerance: float = 1e-4,
) -> list[np.ndarray]:
    """Populations (P, S) per site in the mean field of the others, at every pH of *grid*.
    Each round updates only the sites that moved, or whose neighbours moved, in the last."""
    own = _at_ph(single, protons, grid)
    links = _links(len(single), pairs)
    pop = [_boltzmann(e) for e in own]
    active = {i for i, e in enumerate(own) if e.shape[1] > 1}
    for _ in range(rounds):
        check(cancel)
        moved = set()
        for i in sorted(active):
            field = own[i].copy()
            for j, m in links[i]:
                field += pop[j] @ m.T
            new = 0.5 * pop[i] + 0.5 * _boltzmann(field)
            if np.abs(new - pop[i]).max() >= tolerance:
                moved.add(i)
            pop[i] = new
        if not moved:
            break
        active = {j for i in moved for j, _m in links[i] if own[j].shape[1] > 1} | moved
    return pop


# -- Monte Carlo -----------------------------------------------------------------------------------
def monte_carlo(
    single: list[np.ndarray],
    pairs: dict,
    protons: list[np.ndarray],
    grid: np.ndarray = GRID,
    cancel: CancelToken | None = None,
    sweeps: int = 2000,
    burn: int = 400,
    seed: int = 0,
) -> list[np.ndarray]:
    """Populations (P, S) per site from Metropolis sampling, one replica per pH."""
    rng = np.random.default_rng(seed)
    own = _at_ph(single, protons, grid)
    links = _links(len(single), pairs)
    n_ph = len(grid)
    rows = np.arange(n_ph)
    state = np.stack([np.argmin(e, axis=1) for e in own], axis=1)  # (P, n)
    counts = [np.zeros_like(e) for e in own]
    moving = [i for i, e in enumerate(own) if e.shape[1] > 1]
    movable = set(moving)
    coupled = [(i, j, m) for (i, j), m in pairs.items()
               if i in movable and j in movable and np.ptp(m) > COUPLED]  # fmt: skip

    def local(i: int, s: np.ndarray) -> np.ndarray:
        e = own[i][rows, s]
        for j, m in links[i]:
            e = e + m[s, state[:, j]]
        return e

    for sweep in range(sweeps):
        if sweep % 50 == 0:
            check(cancel)
        for i in rng.permutation(moving).tolist():
            old = state[:, i]
            new = rng.integers(own[i].shape[1], size=n_ph)
            delta = local(i, new) - local(i, old)
            take = rng.random(n_ph) < np.exp(-np.clip(delta, -50, None) / KT)
            state[:, i] = np.where(take, new, old)
        for i, j, m in coupled:  # both at once, across a barrier single moves would not cross
            old_i, old_j = state[:, i].copy(), state[:, j].copy()
            before = local(i, old_i) + local(j, old_j) - m[old_i, old_j]
            state[:, i] = rng.integers(own[i].shape[1], size=n_ph)
            state[:, j] = rng.integers(own[j].shape[1], size=n_ph)
            after = local(i, state[:, i]) + local(j, state[:, j]) - m[state[:, i], state[:, j]]
            take = rng.random(n_ph) < np.exp(-np.clip(after - before, -50, None) / KT)
            state[:, i] = np.where(take, state[:, i], old_i)
            state[:, j] = np.where(take, state[:, j], old_j)
        if sweep >= burn:
            for i in moving:
                np.add.at(counts[i], (rows, state[:, i]), 1)
    out = []
    for i, e in enumerate(own):
        if i in movable:
            out.append(counts[i] / counts[i].sum(axis=1, keepdims=True))
        else:
            out.append(np.ones_like(e))
    return out
