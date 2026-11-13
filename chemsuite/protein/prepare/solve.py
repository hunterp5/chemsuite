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

"""The best state of every site at once, and the structure with those states.

Sites that touch form clusters (Reduce's cliques), each solved on its own: dead-end
elimination (Goldstein's singles criterion) first, then a depth-first branch-and-bound over the
states left, started from a greedy answer. A cluster too large to finish within
:data:`NODE_LIMIT` keeps the best answer found. Each site's *gap* is how much worse its next
best state is with every other site held, for marking the ambiguous ones later.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from ..structure import edit
from . import energy, sites
from .energy import Tables, Weights
from .hydrogens import HYDROGEN, Protonated

NODE_LIMIT = 200_000
FORBIDDEN = 1e6  # the energy of a state a site is held away from
#: forms recorded in :attr:`.hydrogens.Protonated.states` (the textbook ones are not)
RECORDED = frozenset(["HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYM", "TYM", "NH2"])


@dataclass(frozen=True)
class Choice:
    site: sites.Site
    state: int
    gap: float  # kcal/mol-like units, to the next best state of this site (inf when none)

    @property
    def label(self) -> str:
        return self.site.states[self.state].label


def optimise(
    p: Protonated,
    ph: float = 7.0,
    weights: Weights | None = None,
    cancel: CancelToken | None = None,
    titrate: str = sites.OFF,
) -> tuple[Protonated, list[Choice]]:
    """*p* with each flip, His state and rotor at its best; with *titrate* (``ALL``), each
    titratable group in its best protonation state at *ph* too."""
    found = sites.find_sites(p, titrate)
    check(cancel)
    table = energy.tables(p.structure, found, ph)
    check(cancel)
    picked = solve(table, weights or Weights())
    return apply(p, found, picked), choices(table, weights or Weights(), picked)


# -- solving ---------------------------------------------------------------------------------------
def solve(
    table: Tables,
    weights: Weights,
    ph: float | None = None,
    allowed: dict[int, np.ndarray] | None = None,
) -> list[int]:
    """The best state of every site; *allowed* holds some sites to some of their states (a
    mask per site)."""
    single, pairs = table.energies(weights, ph)
    for i, mask in (allowed or {}).items():
        single[i] = np.where(mask, single[i], FORBIDDEN)
    n = len(single)
    links: list[set[int]] = [set() for _ in range(n)]
    for i, j in pairs:
        links[i].add(j)
        links[j].add(i)
    picked = [int(np.argmin(e)) for e in single]
    seen = [False] * n
    for start in range(n):
        if seen[start]:
            continue
        cluster, stack = [], [start]
        seen[start] = True
        while stack:
            i = stack.pop()
            cluster.append(i)
            for j in links[i]:
                if not seen[j]:
                    seen[j] = True
                    stack.append(j)
        if len(cluster) > 1:
            for i, state in _Cluster(sorted(cluster), single, pairs).best().items():
                picked[i] = state
    return picked


def _pair(pairs: dict, i: int, j: int) -> np.ndarray | None:
    if (i, j) in pairs:
        return pairs[(i, j)]
    if (j, i) in pairs:
        return pairs[(j, i)].T
    return None


class _Cluster:
    def __init__(self, members: list[int], single: list[np.ndarray], pairs: dict) -> None:
        self.members = members
        self.single = {i: single[i] for i in members}
        self.links = {i: {j: m for j in members if (m := _pair(pairs, i, j)) is not None}
                      for i in members}  # fmt: skip
        self.alive = {i: np.arange(len(single[i])) for i in members}

    def energy(self, states: dict[int, int]) -> float:
        total = sum(float(self.single[i][s]) for i, s in states.items())
        for i, s in states.items():
            for j, m in self.links[i].items():
                if j > i and j in states:
                    total += float(m[s, states[j]])
        return total

    def eliminate(self) -> None:
        """Goldstein singles: drop state r of site i when some state t beats it whatever the
        neighbours do."""
        changed = True
        while changed:
            changed = False
            for i in self.members:
                alive = self.alive[i]
                if len(alive) < 2:
                    continue
                e = self.single[i][alive]
                gain = e[:, None] - e[None, :]  # (r, t)
                for j, m in self.links[i].items():
                    sub = m[np.ix_(alive, self.alive[j])]
                    gain = gain + (sub[:, None, :] - sub[None, :, :]).min(axis=2)
                np.fill_diagonal(gain, -np.inf)
                dead = (gain > 1e-9).any(axis=1)
                if dead.any() and not dead.all():
                    self.alive[i] = alive[~dead]
                    changed = True

    def greedy(self) -> dict[int, int]:
        states = {i: int(self.alive[i][np.argmin(self.single[i][self.alive[i]])])
                  for i in self.members}  # fmt: skip
        for _ in range(20):
            moved = False
            for i in self.members:
                e = self.single[i][self.alive[i]].copy()
                for j, m in self.links[i].items():
                    e += m[self.alive[i], states[j]]
                best = int(self.alive[i][np.argmin(e)])
                if best != states[i]:
                    states[i], moved = best, True
            if not moved:
                break
        return states

    def best(self) -> dict[int, int]:
        self.eliminate()
        best = self.greedy()
        bound = [self.energy(best)]
        found = [dict(best)]
        order = sorted(self.members, key=lambda i: -len(self.links[i]))
        position = {i: k for k, i in enumerate(order)}
        alive = {i: self.alive[i] for i in order}
        # each site's own energy, plus the least it can get from every site after it
        base = {i: self.single[i][alive[i]] + sum(
                    (m[np.ix_(alive[i], alive[j])].min(axis=1)
                     for j, m in self.links[i].items() if position[j] > position[i]),
                    np.zeros(len(alive[i])))
                for i in order}  # fmt: skip
        # what the sites already chosen give each state of the sites still open
        given = {i: np.zeros(len(alive[i])) for i in order}
        nodes = [0]

        def search(depth: int, states: dict[int, int], so_far: float) -> None:
            nodes[0] += 1
            if nodes[0] > NODE_LIMIT:
                return
            if depth == len(order):
                if so_far < bound[0] - 1e-9:
                    bound[0], found[0] = so_far, dict(states)
                return
            i = order[depth]
            rest = sum(float((base[k] + given[k]).min()) for k in order[depth + 1 :])
            options = np.argsort(self.single[i][alive[i]] + given[i])
            for k in options.tolist():
                r = int(alive[i][k])
                e = float(self.single[i][r] + given[i][k])
                if so_far + float(base[i][k] - self.single[i][r]) + e + rest >= bound[0] - 1e-9:
                    continue
                states[i] = r
                for j, m in self.links[i].items():
                    if position[j] > depth:
                        given[j] += m[r, alive[j]]
                search(depth + 1, states, so_far + e)
                for j, m in self.links[i].items():
                    if position[j] > depth:
                        given[j] -= m[r, alive[j]]
                del states[i]

        search(0, {}, 0.0)
        return found[0]


def choices(
    table: Tables, weights: Weights, picked: list[int], ph: float | None = None
) -> list[Choice]:
    single, pairs = table.energies(weights, ph)
    out = []
    for i, site in enumerate(table.sites):
        e = single[i].copy()
        for j in range(len(single)):
            m = _pair(pairs, i, j)
            if m is not None:
                e += m[:, picked[j]]
        others = np.delete(e, picked[i])
        gap = float(others.min() - e[picked[i]]) if len(others) else float("inf")
        out.append(Choice(site, picked[i], gap))
    return out


# -- writing the answer back -----------------------------------------------------------------------
def apply(p: Protonated, found: list[sites.Site], picked: list[int]) -> Protonated:
    """*p* with the state *picked* for each site: atoms moved, hydrogens added or removed and
    formal charges set (other atoms keep their uids)."""
    s = p.structure
    coords = s.coords.copy()
    charge = s.charge.copy()
    states = dict(p.states)
    gone: list[int] = []
    new: list[tuple[int, str, int]] = []  # residue, name, parent uid
    new_xyz: list[np.ndarray] = []
    for site, k in zip(found, picked, strict=True):
        state = site.states[k]
        by_name = {str(s.name[i]): int(i) for i in site.moving}
        span = s.atoms_of(site.residue)
        residue_atoms = {str(s.name[i]): i for i in range(span.start, span.stop)}
        for m, name in enumerate(state.names):
            if name in by_name:
                coords[:, by_name[name]] = state.xyz[:, m]
            else:
                parent = residue_atoms.get(state.parents[m], -1) if state.parents[m] else -1
                if parent >= 0:
                    new.append((site.residue, name, int(s.uid[parent])))
                    new_xyz.append(state.xyz[:, m])
        gone += [i for name, i in by_name.items()
                 if name not in state.names and s.element[i] == HYDROGEN]  # fmt: skip
        for name, q in zip(state.names, state.formal, strict=False):
            if name in residue_atoms:
                charge[residue_atoms[name]] = q
        if state.form in RECORDED:
            states[site.residue] = state.form
    out = replace(s, coords=coords, charge=charge)
    if gone:
        mask = np.zeros(out.n_atoms, dtype=bool)
        mask[gone] = True
        out = edit.delete(out, mask)
    if new:
        parents = out.indices(np.array([u for _r, _n, u in new]))
        out = edit.add_atoms(
            out,
            residues=np.array([r for r, _n, _u in new]),
            names=np.array([n for _r, n, _u in new]),
            elements=np.full(len(new), HYDROGEN),
            xyz=np.stack(new_xyz, axis=1),
            bonds=np.stack([parents, out.n_atoms + np.arange(len(new))], axis=1),
        )
    return replace(p, structure=out, states=states)
