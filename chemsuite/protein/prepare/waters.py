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

"""Pocket waters as sites: each water near a ligand takes hydrogens, in one of a few
orientations scored with everything around it.

A water whose O is within :data:`REACH` of a ligand's heavy atom becomes a site. Its
orientations (at most :data:`MAX_STATES`) are, in order:

- both hydrogens towards two polar neighbours within :data:`NEIGHBOURS` (each pair);
- one hydrogen towards one polar neighbour, the other turned about it in thirds;
- :data:`SPREAD` fixed orientations, so a water with few neighbours still has a choice.

Orientations closer than :data:`SAME` (degrees, both hydrogens) to one already listed are left
out. The geometry is TIP3P's (O–H 0.9572 Å, H–O–H 104.52°). Other waters keep no hydrogens and
are scored as both donor and acceptor (:mod:`.energy`'s ``water`` term).
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from ..structure.model import LIGAND, WATER
from . import sites
from .hydrogens import HYDROGEN, Protonated
from .sites import POLAR, Site, State

REACH = 5.0  # Å from a ligand heavy atom
NEIGHBOURS = 3.5  # Å: polar atoms a hydrogen may point at
MAX_STATES = 30
SPREAD = 12
SAME = 20.0  # degrees
OH, HOH = 0.9572, np.radians(104.52)


def find(p: Protonated) -> list[Site]:
    """A site for each water near a ligand that came without hydrogens."""
    s = p.structure
    active = s.active()
    x = s.coords.astype(np.float64)
    ligand = np.flatnonzero(active & (s.atom_kind == LIGAND) & (s.element != HYDROGEN))
    if not len(ligand):
        return []
    near_ligand = cKDTree(x[0, ligand])
    polar = np.flatnonzero(active & np.isin(s.element, list(POLAR)))
    polar_tree = cKDTree(x[0, polar])
    graph = sites.bond_graph(s)
    found = []
    for r in np.flatnonzero(s.res_kind == WATER).tolist():
        span = s.atoms_of(r)
        atoms = [i for i in range(span.start, span.stop) if active[i]]
        oxygen = [i for i in atoms if s.element[i] == 8]
        if len(oxygen) != 1 or len(atoms) != 1:
            continue  # it has hydrogens already, or is not a plain water
        o = oxygen[0]
        if not near_ligand.query_ball_point(x[0, o], REACH):
            continue
        partners = [int(polar[k]) for k in polar_tree.query_ball_point(x[0, o], NEIGHBOURS)
                    if polar[k] != o]  # fmt: skip
        directions = _orientations([_unit(x[0, j] - x[0, o]) for j in partners])
        states = tuple(_state(x[:, o], h1, h2, k) for k, (h1, h2) in enumerate(directions))
        site = sites.make_site(s, r, sites.WATER, np.array([o]), [o], list(states), graph)
        found.append(site)
    return found


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def _pair(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    """Two O–H directions at H–O–H, in the plane of *a* and *b*, about their bisector."""
    bisector = a + b
    if np.linalg.norm(bisector) < 1e-6:
        return None
    bisector = _unit(bisector)
    across = a - b - np.dot(a - b, bisector) * bisector
    if np.linalg.norm(across) < 1e-6:
        return None
    across = _unit(across)
    half = HOH / 2
    return (np.cos(half) * bisector + np.sin(half) * across,
            np.cos(half) * bisector - np.sin(half) * across)  # fmt: skip


def _around(first: np.ndarray, turn: float) -> np.ndarray:
    """A direction at H–O–H from *first*, turned by *turn* about it."""
    helper = np.array([1.0, 0.0, 0.0]) if abs(first[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = _unit(np.cross(first, helper))
    v = np.cross(first, u)
    side = np.cos(turn) * u + np.sin(turn) * v
    return np.cos(HOH) * first + np.sin(HOH) * side


def _spread(count: int) -> list[np.ndarray]:
    """*count* directions spread over the sphere (a Fibonacci lattice)."""
    k = np.arange(count) + 0.5
    phi = np.arccos(1 - 2 * k / count)
    theta = np.pi * (1 + 5**0.5) * k
    return [np.array([np.cos(t) * np.sin(f), np.sin(t) * np.sin(f), np.cos(f)])
            for f, t in zip(phi, theta, strict=True)]  # fmt: skip


def _orientations(towards: list[np.ndarray]) -> list[tuple[np.ndarray, np.ndarray]]:
    candidates = []
    for i in range(len(towards)):
        for j in range(i + 1, len(towards)):
            pair = _pair(towards[i], towards[j])
            if pair is not None:
                candidates.append(pair)
    thirds = (0.0, 2 * np.pi / 3, 4 * np.pi / 3)
    candidates += [(d, _around(d, turn)) for d in towards for turn in thirds]
    candidates += [(d, _around(d, 0.0)) for d in _spread(SPREAD)]
    out: list[tuple[np.ndarray, np.ndarray]] = []
    limit = np.cos(np.radians(SAME))
    for h1, h2 in candidates:
        if any(_alike(h1, h2, a, b, limit) for a, b in out):
            continue
        out.append((h1, h2))
        if len(out) == MAX_STATES:
            break
    return out


def _alike(h1: np.ndarray, h2: np.ndarray, a: np.ndarray, b: np.ndarray, limit: float) -> bool:
    """Whether two orientations put their hydrogens (in either order) within the limit."""
    same = np.dot(h1, a) > limit and np.dot(h2, b) > limit
    swapped = np.dot(h1, b) > limit and np.dot(h2, a) > limit
    return bool(same or swapped)


def _state(oxygen: np.ndarray, h1: np.ndarray, h2: np.ndarray, k: int) -> State:
    """Water with its hydrogens along *h1* and *h2* from the O at *oxygen* (k, 3)."""
    xyz = np.stack([oxygen, oxygen + OH * h1, oxygen + OH * h2], axis=1)  # (k, 3, 3)
    donor = np.full((3, 3), np.nan)
    donor[1:] = oxygen[0]
    return State(f"orientation {k + 1}", ("O", "H1", "H2"), ("", "O", "O"),
                 np.array([8, HYDROGEN, HYDROGEN]), xyz, donor, formal=(0, 0, 0),
                 partial=(0.0, 0.0, 0.0))  # fmt: skip
