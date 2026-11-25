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

"""Seeded 2D layouts of compound graphs: a force layout, its multilevel form, and packing.

**Force layout.** Fruchterman, T. M. J. and Reingold, E. M. (1991). "Graph drawing by
force-directed placement." *Softw. Pract. Exp.* 21, 1129–1164. Edges pull their ends together
by d²/k, nodes push each other apart by k²/d, and each step moves a node at most the current
temperature, which cools. Repulsion acts only between nodes closer than a cutoff (found with a
k-d tree), so a step costs about n log n rather than n².

**Multilevel.** Walshaw, C. (2003). "A multilevel algorithm for force-directed graph-drawing."
*J. Graph Algorithms Appl.* 7, 253–285. The graph is coarsened by merging matched pairs of
neighbours (and free leaves into their neighbour) until it is small, laid out there, and each
level is then placed from the one above and refined with local forces. A node's weight (how
many nodes it stands for) scales the push it gives. This lays out trees of hundreds of
thousands of nodes, as TMAP's layouts (OGDF's FM³: Hachul, S. and Jünger, M. (2004),
"Drawing large graphs with a potential-field-based multilevel algorithm," *Graph Drawing*,
LNCS 3383, 285–295) do, without a compiled dependency.

How far the push reaches decides the picture. A short reach leaves branches free to pile into
clumps and fold over each other; an unlimited one drives everything to a hollow ring with
stretched edges. Sixteen edge lengths on the coarse levels, where the branches are arranged
and nodes are few, and eight on the finest spread a tree over the page as FM³ does, at a cost
that grows with the square of the reach. Each finer level starts every node on its group's
place moved toward the groups its edges lead to, so the nodes of a group start on their
neighbors' sides rather than in a heap.

**Packing.** Components are laid out alone and placed in rows, largest first.

Every layout is seeded, so the same graph always gets the same picture.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

Progress = Callable[[int, int], None]

SEED = 42
COARSEST = 64  # nodes at which coarsening stops
MIN_SHRINK = 0.9  # stop coarsening once a level keeps more than this share of its nodes
CUTOFF = 8.0  # repulsion reaches this many edge lengths on the finest level
COARSE_CUTOFF = 16.0  # and this many on the coarser ones
FINE_STEPS = 30
COARSE_STEPS = 60
FIRST_STEPS = 300
COOLING = 0.93
TOWARD = 0.5  # how far a placed node moves toward its neighbors outside its group
NOISE = 0.3  # and the jitter that parts the rest, in edge lengths
GAP = 2.0  # edge lengths between packed components


# -- one level ----------------------------------------------------------------------------------
def _pairs(pos: np.ndarray, reach: float) -> tuple[np.ndarray, np.ndarray]:
    if math.isinf(reach):  # small graphs: every pair
        return np.triu_indices(len(pos), 1)
    pairs = cKDTree(pos).query_pairs(reach, output_type="ndarray")
    if not len(pairs):
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    return pairs[:, 0], pairs[:, 1]


def _forces(
    pos: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    weight: np.ndarray,
    k: float,
    reach: float,
) -> np.ndarray:
    n = len(pos)
    disp = np.zeros((n, 2))
    i, j = _pairs(pos, reach)
    if len(i):
        delta = pos[i] - pos[j]
        d2 = np.maximum((delta**2).sum(axis=1), 1e-12)
        push = (k * k / d2)[:, None] * delta  # k²/d along the unit vector
        for axis in (0, 1):
            disp[:, axis] += np.bincount(i, push[:, axis] * weight[j], minlength=n)
            disp[:, axis] -= np.bincount(j, push[:, axis] * weight[i], minlength=n)
    if len(a):
        delta = pos[a] - pos[b]
        d = np.sqrt(np.maximum((delta**2).sum(axis=1), 1e-12))
        pull = (d / k)[:, None] * delta  # d²/k along the unit vector
        for axis in (0, 1):
            disp[:, axis] -= np.bincount(a, pull[:, axis], minlength=n)
            disp[:, axis] += np.bincount(b, pull[:, axis], minlength=n)
    return disp


def force(
    pos: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    *,
    weight: np.ndarray | None = None,
    k: float = 1.0,
    cutoff: float = CUTOFF,
    steps: int = FIRST_STEPS,
    temperature: float | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """*pos* (``n × 2``) after *steps* force steps; edges ``(a, b)`` want length *k*.
    Repulsion reaches ``cutoff · k`` (``inf``: every pair)."""
    pos = np.array(pos, dtype=np.float64)
    if len(pos) < 2:
        return pos
    weight = np.ones(len(pos)) if weight is None else np.asarray(weight, dtype=np.float64)
    t = temperature if temperature is not None else k * math.sqrt(len(pos)) / 2
    reach = cutoff * k
    for _ in range(steps):
        check(cancel)
        disp = _forces(pos, a, b, weight, k, reach)
        length = np.sqrt((disp**2).sum(axis=1))
        scale = np.minimum(length, t) / np.maximum(length, 1e-12)
        pos += disp * scale[:, None]
        t *= COOLING
    return pos


# -- coarsening ---------------------------------------------------------------------------------
def _coarsen(
    n: int, a: np.ndarray, b: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, int]:
    """``(group per node, groups)``: matched neighbours share a group, and every free leaf
    joins its neighbour's."""
    partner = np.full(n, -1, dtype=np.int64)
    for _ in range(3):  # each free node bids for its free edge of lowest random key
        free = (partner[a] < 0) & (partner[b] < 0)
        ea, eb = a[free], b[free]
        if not len(ea):
            break
        key = rng.random(len(ea))
        best = np.full(n, np.inf)
        np.minimum.at(best, ea, key)
        np.minimum.at(best, eb, key)
        won = (best[ea] == key) & (best[eb] == key)
        partner[ea[won]], partner[eb[won]] = eb[won], ea[won]
    nodes = np.arange(n)
    lead = np.where(partner >= 0, np.minimum(nodes, partner), nodes)
    degree = np.bincount(a, minlength=n) + np.bincount(b, minlength=n)
    leaf_a = (degree[a] == 1) & (partner[a] < 0)  # a is a free leaf: it joins b
    leaf_b = (degree[b] == 1) & (partner[b] < 0) & ~leaf_a
    lead[a[leaf_a]] = lead[b[leaf_a]]
    lead[b[leaf_b]] = lead[a[leaf_b]]
    _, group = np.unique(lead, return_inverse=True)
    return group, int(group.max()) + 1


def multilevel(
    n: int,
    a: np.ndarray,
    b: np.ndarray,
    *,
    seed: int = SEED,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """``n × 2`` positions of a connected graph whose edges want length 1."""
    rng = np.random.default_rng(seed)
    if n <= 1:
        return np.zeros((n, 2))
    levels: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []  # (group, a, b) per level
    weights = [np.ones(n)]
    na, ea, eb = n, a, b
    while na > COARSEST:
        check(cancel)
        group, ng = _coarsen(na, ea, eb, rng)
        if ng > MIN_SHRINK * na:
            break
        levels.append((group, ea, eb))
        ca, cb = group[ea], group[eb]
        keep = ca != cb
        ca, cb = np.minimum(ca[keep], cb[keep]), np.maximum(ca[keep], cb[keep])
        _, first = np.unique(ca * ng + cb, return_index=True)
        ea, eb, na = ca[first], cb[first], ng
        weights.append(np.bincount(group, weights=weights[-1], minlength=ng))
    steps = len(levels) + 1
    k = math.sqrt(n / na)  # a coarse node stands for about n / na nodes
    angle = rng.random(na) * 2 * math.pi
    spread = k * math.sqrt(na) * np.sqrt(rng.random(na))
    pos = np.column_stack([spread * np.cos(angle), spread * np.sin(angle)])
    small = na <= COARSEST
    pos = force(
        pos,
        ea,
        eb,
        weight=weights[-1],
        k=k,
        cutoff=math.inf if small else COARSE_CUTOFF,
        steps=FIRST_STEPS if small else COARSE_STEPS,
        cancel=cancel,
    )
    if progress is not None:
        progress(1, steps)
    for done, (group, fa, fb) in enumerate(reversed(levels)):
        nf = len(group)
        k = math.sqrt(n / nf)
        pos = _place(pos, group, fa, fb, k, rng)
        pos = force(
            pos,
            fa,
            fb,
            weight=weights[len(levels) - done - 1],
            k=k,
            cutoff=CUTOFF if nf == n else COARSE_CUTOFF,
            steps=FINE_STEPS if nf == n else COARSE_STEPS,
            temperature=2 * k,
            cancel=cancel,
        )
        if progress is not None:
            progress(done + 2, steps)
    return pos


def _place(
    coarse: np.ndarray,
    group: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    k: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """A finer level's positions from the coarser one: each node starts on its group's
    position, moved part of the way toward the groups its edges lead out to, so the nodes of
    a group start on the sides their neighbors are on rather than in a heap."""
    n = len(group)
    start = coarse[group]
    out = group[a] != group[b]
    oa, ob = a[out], b[out]
    count = np.bincount(oa, minlength=n) + np.bincount(ob, minlength=n)
    toward = np.zeros((n, 2))
    for axis in (0, 1):
        toward[:, axis] = np.bincount(oa, start[ob, axis], minlength=n) + np.bincount(
            ob, start[oa, axis], minlength=n
        )
    has = count > 0
    shift = np.zeros((n, 2))
    shift[has] = toward[has] / count[has, None] - start[has]
    return start + TOWARD * shift + rng.normal(scale=NOISE * k, size=(n, 2))


# -- many components ----------------------------------------------------------------------------
def layout(
    n: int,
    a: np.ndarray,
    b: np.ndarray,
    labels: np.ndarray,
    *,
    seed: int = SEED,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """``n × 2`` positions: each component (*labels*, 0 the largest) laid out alone with
    :func:`multilevel`, then :func:`pack`\\ ed."""
    pos = np.zeros((n, 2))
    if n == 0:
        return pos
    count = int(labels.max()) + 1
    order = np.argsort(labels, kind="stable")
    bounds = np.searchsorted(labels[order], np.arange(count + 1))
    edge_label = labels[a]
    edge_order = np.argsort(edge_label, kind="stable")
    edge_bounds = np.searchsorted(edge_label[edge_order], np.arange(count + 1))
    sizes = np.diff(bounds)
    local = np.empty(n, dtype=np.int64)
    for c in np.flatnonzero(sizes > COARSEST):  # large components one by one
        check(cancel)
        nodes = order[bounds[c] : bounds[c + 1]]
        edges = edge_order[edge_bounds[c] : edge_bounds[c + 1]]
        local[nodes] = np.arange(len(nodes))
        pos[nodes] = multilevel(
            len(nodes), local[a[edges]], local[b[edges]], seed=seed + int(c), cancel=cancel
        )
        if progress is not None:
            progress(int(bounds[c + 1]), n)
    small = (sizes[labels] > 1) & (sizes[labels] <= COARSEST)
    if small.any():
        pos[small] = _small_components(a, b, labels, small, seed, cancel)
    if progress is not None:
        progress(n, n)
    return pack(pos, labels)


def _small_components(
    a: np.ndarray,
    b: np.ndarray,
    labels: np.ndarray,
    small: np.ndarray,
    seed: int,
    cancel: CancelToken | None,
) -> np.ndarray:
    """Positions of the nodes in *small* components, laid out together in one force run: each
    component starts in a cell of its own, farther from the others than repulsion reaches."""
    rows = np.flatnonzero(small)
    local = np.full(len(labels), -1, dtype=np.int64)
    local[rows] = np.arange(len(rows))
    keep = small[a]  # an edge lies within one component
    ea, eb = local[a[keep]], local[b[keep]]
    _, cell = np.unique(labels[rows], return_inverse=True)
    side = math.ceil(math.sqrt(int(cell.max()) + 1))
    spacing = 4 * math.sqrt(COARSEST) + CUTOFF * 4
    centre = np.column_stack([cell % side, cell // side]) * spacing
    rng = np.random.default_rng(seed)
    pos = centre + rng.normal(scale=1.0, size=(len(rows), 2))
    return force(pos, ea, eb, steps=FIRST_STEPS, temperature=math.sqrt(COARSEST) / 2, cancel=cancel)


def pack(pos: np.ndarray, labels: np.ndarray, gap: float = GAP) -> np.ndarray:
    """Move each component (*labels*, 0 the largest) so they sit in rows, largest first, the
    rows about as wide as the whole is tall."""
    if not len(pos):
        return pos
    count = int(labels.max()) + 1
    lo = np.full((count, 2), np.inf)
    hi = np.full((count, 2), -np.inf)
    for axis in (0, 1):
        np.minimum.at(lo[:, axis], labels, pos[:, axis])
        np.maximum.at(hi[:, axis], labels, pos[:, axis])
    size = hi - lo + gap
    width = max(math.sqrt(float((size[:, 0] * size[:, 1]).sum())) * 1.2, float(size[:, 0].max()))
    corner = np.zeros((count, 2))
    x = y = row = 0.0
    for c in range(count):  # the labels run largest first already
        if x > 0 and x + size[c, 0] > width:
            x, y, row = 0.0, y - row, 0.0
        corner[c] = (x, y - size[c, 1])
        x += size[c, 0]
        row = max(row, size[c, 1])
    return pos - lo[labels] + corner[labels]
