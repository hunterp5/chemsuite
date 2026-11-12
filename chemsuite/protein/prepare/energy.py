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

"""The score of each site state: contacts with the fixed atoms and with the other sites.

The energy is a weighted sum of terms (:data:`TERMS`), kept apart so the weights can be fitted
(``scripts/review/prepare_review.py fit``):

- ``hbond``: minus the strength of each hydrogen bond, a polar H to an acceptor (O, or an N
  with no H): 1 from 1.5 to 2.2 Å, falling to 0 at 2.7 Å, times an angle factor that is 1 for
  a straight D–H···A (≥150°) and 0 below 90°;
- ``clash``: overlap² of an H with an H (closer than 2.0 Å) or with a heavy atom (2.4 Å; 1.6 Å
  for a polar H on an acceptor);
- ``lone``: two acceptors without hydrogens facing each other (closer than 3.1 Å, lone pair
  on lone pair), as an amide O on a carbonyl O after a wrong flip;
- ``water``: minus each contact of an acceptor without hydrogens with a water O (waters carry
  no hydrogens here, so they may donate);
- ``coulomb``: charge-charge energy, 332·qᵢqⱼ/(εr) with ε = 4r, between a HIP ring (+½ on
  each N) and the charged groups within :data:`COULOMB_REACH` (a carboxylate's charge shared by
  its two O, Arg's by its three N, metal ions +2), and between two HIP;
- ``metal``: an acceptor without hydrogens within :data:`METAL_REACH` of a metal ion counts
  −1, a carbon there +1 (a His ring or an amide must point an N or O at the metal);
- ``burial``: for HIP, the heavy atoms within :data:`BURIAL_REACH` of the ring per hundred:
  a buried charge loses its water;
- priors: ``flip`` (1 for a flipped state, so a flip must win by a margin, as in Reduce),
  ``hip`` and ``hid`` (1 for those His states), and ``ph`` (1.364·(pH − 6.5) for HIP, the free
  energy of protonating a His of solution pKa 6.5).

Atoms within :data:`~.sites.EXCLUDED_BONDS` bonds of a site's anchors are not scored
against it.
"""

from __future__ import annotations

from dataclasses import astuple, dataclass, fields

import numpy as np
from scipy.spatial import cKDTree

from ..structure.model import ION, WATER, Structure
from .hydrogens import HYDROGEN
from .sites import POLAR, Site, State

TERMS = ("hbond", "clash", "lone", "water", "flip", "hip", "hid", "ph", "coulomb", "burial",
         "metal")  # fmt: skip
HBOND, CLASH, LONE, WATER_TERM, FLIP, HIP, HID, PH, COULOMB, BURIAL, METAL = range(len(TERMS))
METAL_REACH = 2.8
REACH = 4.0  # Å: atoms further than this from every state do not touch the site
COULOMB_REACH = 12.0
BURIAL_REACH = 10.0
HIS_PKA = 6.5
#: groups whose charge is shared by their atoms, by residue and atom names
SHARED = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "ARG": ("NE", "NH1", "NH2"),
          "HIS": ("ND1", "NE2")}  # fmt: skip
ION_CHARGE = {3: 1, 11: 1, 19: 1, 37: 1, 55: 1}  # other metal ions count as +2
HALIDES = (9, 17, 35, 53)


@dataclass(frozen=True)
class Weights:
    """Fitted on half the neutron set (2026-11-12, docs/prepare/p1-review.md)."""

    hbond: float = 1.0
    clash: float = 160.0
    lone: float = 16.0
    water: float = 1.0
    flip: float = 0.0
    hip: float = 3.0
    hid: float = 0.8
    ph: float = 1.0
    coulomb: float = 1.0
    burial: float = 1.0
    metal: float = 2.0

    def vector(self) -> np.ndarray:
        assert tuple(f.name for f in fields(self)) == TERMS
        return np.array(astuple(self), dtype=np.float64)


@dataclass(frozen=True)
class Batch:
    """Atoms of several states (or of the environment, as one state), padded to one size."""

    xyz: np.ndarray  # (S, M, 3)
    hydrogen: np.ndarray  # (S, M) bool
    polar_h: np.ndarray
    acceptor: np.ndarray
    lone: np.ndarray
    water: np.ndarray
    donor: np.ndarray  # (S, M, 3)
    valid: np.ndarray  # (S, M)
    charge: np.ndarray  # (S, M)
    metal: np.ndarray  # (S, M) bool: a metal ion

    def take(self, atoms: np.ndarray) -> Batch:
        """The environment's atoms *atoms* (a batch of one state)."""
        arrays = (self.xyz, self.hydrogen, self.polar_h, self.acceptor, self.lone, self.water,
                  self.donor, self.valid, self.charge, self.metal)  # fmt: skip
        return Batch(*(v[:, atoms] for v in arrays))


def of_states(states: tuple[State, ...]) -> Batch:
    size = max(len(st.names) for st in states)
    n = len(states)
    out = {k: np.zeros((n, size), dtype=bool) for k in ("hydrogen", "polar", "acc", "valid")}
    xyz = np.zeros((n, size, 3))
    donor = np.full((n, size, 3), np.nan)
    charge = np.zeros((n, size))
    for k, st in enumerate(states):
        m = len(st.names)
        xyz[k, :m] = st.xyz[0]
        donor[k, :m] = st.donor
        out["hydrogen"][k, :m] = st.elements == HYDROGEN
        out["polar"][k, :m] = st.polar_h
        out["acc"][k, :m] = st.acceptor
        out["valid"][k, :m] = True
        if st.tautomer == "HIP":
            for ring_n in ("ND1", "NE2"):
                charge[k, st.names.index(ring_n)] = 0.5
    water = np.zeros((n, size), dtype=bool)
    return Batch(xyz, out["hydrogen"], out["polar"], out["acc"], out["acc"].copy(), water, donor,
                 out["valid"], charge, np.zeros((n, size), dtype=bool))  # fmt: skip


@dataclass(frozen=True)
class Environment:
    atoms: np.ndarray  # structure indices of the fixed atoms
    batch: Batch
    tree: cKDTree


def environment(s: Structure, sites: list[Site]) -> Environment:
    """Every atom shown that no site moves."""
    moving = np.zeros(s.n_atoms, dtype=bool)
    for site in sites:
        moving[site.moving] = True
    atoms = np.flatnonzero(s.active() & ~moving)
    hydrogen = s.element == HYDROGEN
    h_on = np.zeros(s.n_atoms, dtype=np.int64)
    heavy_degree = np.zeros(s.n_atoms, dtype=np.int64)
    donor = np.full((s.n_atoms, 3), np.nan)
    polar_h = np.zeros(s.n_atoms, dtype=bool)
    x = s.xyz().astype(np.float64)
    for a, b in s.bonds.tolist():
        for h, p in ((a, b), (b, a)):
            if hydrogen[h] and not hydrogen[p]:
                h_on[p] += 1
                if int(s.element[p]) in POLAR:
                    polar_h[h] = True
                    donor[h] = x[p]
            elif not hydrogen[h] and not hydrogen[p]:
                heavy_degree[h] += 1
    element = s.element
    acceptor = (element == 8) | ((element == 7) & (h_on == 0) & (heavy_degree <= 2))
    water = s.atom_kind == WATER
    lone = acceptor & (h_on == 0) & ~water
    metal = (s.atom_kind == ION) & (element > 2) & ~np.isin(element, HALIDES)
    batch = Batch(x[atoms][None], hydrogen[atoms][None], polar_h[atoms][None],
                  acceptor[atoms][None], lone[atoms][None], water[atoms][None],
                  donor[atoms][None], np.ones((1, len(atoms)), dtype=bool),
                  group_charges(s)[atoms][None], metal[atoms][None])  # fmt: skip
    return Environment(atoms, batch, cKDTree(x[atoms]))


def group_charges(s: Structure) -> np.ndarray:
    """Each atom's formal charge, a group's shared among its atoms, and metal ions charged."""
    q = s.charge.astype(np.float64)
    for name, atoms in SHARED.items():
        for r in np.flatnonzero(s.res_name == name).tolist():
            span = s.atoms_of(r)
            group = [i for i in range(span.start, span.stop) if s.name[i] in atoms]
            if group:
                q[group] = q[group].sum() / len(group)
    ions = np.flatnonzero((s.atom_kind == ION) & (s.element > 2) & (q == 0)
                          & ~np.isin(s.element, HALIDES))  # fmt: skip
    q[ions] = [ION_CHARGE.get(int(e), 2) for e in s.element[ions]]
    return q


def coulomb(a: Batch, b: Batch) -> np.ndarray:
    """Charge-charge energy of every state of *a* with every state of *b*, ε = 4r: (Sa, Sb)."""
    pa = (slice(None), None, slice(None), None)
    pb = (None, slice(None), None, slice(None))
    d = np.maximum(np.linalg.norm(b.xyz[pb] - a.xyz[pa], axis=-1), 2.5)
    return _sum(a.charge[pa] * b.charge[pb] * 83.0 / d**2)


# -- terms ---------------------------------------------------------------------------------
def contact(a: Batch, b: Batch) -> np.ndarray:
    """The contact terms of every state of *a* with every state of *b*: (Sa, Sb, len(TERMS))."""
    out = np.zeros((a.xyz.shape[0], b.xyz.shape[0], len(TERMS)))
    if not a.xyz.shape[1] or not b.xyz.shape[1]:
        return out
    pa = (slice(None), None, slice(None), None)  # (Sa, 1, Ma, 1)
    pb = (None, slice(None), None, slice(None))  # (1, Sb, 1, Mb)
    d = np.linalg.norm(b.xyz[pb] - a.xyz[pa], axis=-1)  # (Sa, Sb, Ma, Mb)
    both = a.valid[pa] & b.valid[pb]
    out[..., HBOND] = -(_hbond(a, b) + _hbond(b, a).T)
    ha, hb = a.hydrogen[pa], b.hydrogen[pb]
    polar_pair = (a.polar_h[pa] & b.acceptor[pb]) | (b.polar_h[pb] & a.acceptor[pa])
    limit = np.where(ha & hb, 2.0, np.where(polar_pair, 1.6, 2.4))
    touching = both & (ha | hb)
    out[..., CLASH] = _sum(np.where(touching, np.clip(limit - d, 0, None) ** 2, 0))
    lone = both & a.lone[pa] & b.lone[pb] & ~a.water[pa] & ~b.water[pb]
    out[..., LONE] = _sum(np.where(lone, np.clip(3.1 - d, 0, None), 0))
    wet = both & ((a.lone[pa] & ~a.water[pa] & b.water[pb]) | (b.lone[pb] & ~b.water[pb]
                                                                & a.water[pa]))  # fmt: skip
    near = np.clip((3.6 - d) / 0.4, 0, 1) * (d > 2.4)
    out[..., WATER_TERM] = -_sum(np.where(wet, near, 0))
    bound = both & (d < METAL_REACH)
    sign_a = np.where(a.lone, -1.0, np.where(a.hydrogen, 0.0, 1.0))
    sign_b = np.where(b.lone, -1.0, np.where(b.hydrogen, 0.0, 1.0))
    out[..., METAL] = _sum(
        np.where(bound & b.metal[pb], sign_a[pa], 0) + np.where(bound & a.metal[pa], sign_b[pb], 0)
    )
    return out


def _hbond(a: Batch, b: Batch) -> np.ndarray:
    """Hydrogen bonds from the polar H of *a* to the acceptors of *b*: (Sa, Sb)."""
    pa = (slice(None), None, slice(None), None)
    pb = (None, slice(None), None, slice(None))
    pairs = a.valid[pa] & b.valid[pb] & a.polar_h[pa] & b.acceptor[pb]
    if not pairs.any():
        return np.zeros((a.xyz.shape[0], b.xyz.shape[0]))
    d = np.linalg.norm(b.xyz[pb] - a.xyz[pa], axis=-1)
    reach = np.clip((2.7 - d) / 0.5, 0, 1) * (d >= 1.5)
    h = a.xyz[pa]
    to_donor = a.donor[pa] - h
    to_acceptor = b.xyz[pb] - h
    with np.errstate(invalid="ignore", divide="ignore"):
        lengths = np.linalg.norm(to_donor, axis=-1) * np.linalg.norm(to_acceptor, axis=-1)
        cos = (to_donor * to_acceptor).sum(-1) / lengths
    theta = np.degrees(np.arccos(np.clip(np.nan_to_num(cos, nan=1.0), -1, 1)))
    straight = np.clip((theta - 90.0) / 60.0, 0, 1) ** 2
    return _sum(np.where(pairs, reach * straight, 0))


def _sum(v: np.ndarray) -> np.ndarray:
    return v.sum(axis=(2, 3))


def priors(site: Site, ph: float) -> np.ndarray:
    out = np.zeros((len(site.states), len(TERMS)))
    for k, st in enumerate(site.states):
        out[k, FLIP] = st.flipped
        out[k, HIP] = st.tautomer == "HIP"
        out[k, HID] = st.tautomer == "HID"
        out[k, PH] = 1.364 * (ph - HIS_PKA) if st.tautomer == "HIP" else 0.0
    return out


# -- every table a structure needs ---------------------------------------------------------
@dataclass(frozen=True)
class Tables:
    sites: list[Site]
    single: list[np.ndarray]  # per site: (S, T)
    pairs: dict[tuple[int, int], np.ndarray]  # (i, j), i < j: (Si, Sj, T)

    def energies(self, weights: Weights) -> tuple[list[np.ndarray], dict]:
        w = weights.vector()
        return [t @ w for t in self.single], {k: t @ w for k, t in self.pairs.items()}


def _charged(site: Site, batch: Batch, env: Environment) -> np.ndarray:
    """The Coulomb and burial terms of a site whose states carry charges."""
    out = np.zeros((len(site.states), len(TERMS)))
    if not batch.charge.any():
        return out
    centre = batch.xyz[batch.valid].mean(axis=0)
    near = [j for j in env.tree.query_ball_point(centre, COULOMB_REACH)
            if int(env.atoms[j]) not in site.excluded]  # fmt: skip
    charged = np.array([j for j in near if env.batch.charge[0, j] != 0], dtype=np.int64)
    if len(charged):
        out[:, COULOMB] = coulomb(batch, env.batch.take(charged))[:, 0]
    heavy = sum(1 for j in env.tree.query_ball_point(centre, BURIAL_REACH)
                if not env.batch.hydrogen[0, j])  # fmt: skip
    out[:, BURIAL] = (batch.charge.sum(axis=1) > 0) * heavy / 100.0
    return out


def tables(s: Structure, sites: list[Site], ph: float = 7.0) -> Tables:
    env = environment(s, sites)
    batches = [of_states(site.states) for site in sites]
    single = []
    for site, batch in zip(sites, batches, strict=True):
        points = batch.xyz[batch.valid]
        near = sorted({j for found in env.tree.query_ball_point(points, REACH) for j in found})
        keep = np.array([j for j in near if int(env.atoms[j]) not in site.excluded],
                        dtype=np.int64)  # fmt: skip
        terms = contact(batch, env.batch.take(keep))[:, 0] if len(keep) else 0.0
        single.append(priors(site, ph) + terms + _charged(site, batch, env))
    pairs = {}
    if sites:
        centres = np.array([b.xyz[b.valid].mean(axis=0) for b in batches])
        radius = np.array([np.linalg.norm(b.xyz[b.valid] - c, axis=-1).max()
                           for b, c in zip(batches, centres, strict=True)])  # fmt: skip
        tree = cKDTree(centres)
        charged = [b.charge.any() for b in batches]
        for i, j in sorted(tree.query_pairs(2 * radius.max() + COULOMB_REACH)):
            if sites[i].excluded.intersection(sites[j].moving.tolist()):
                continue  # too close in the bond graph
            apart = np.linalg.norm(centres[i] - centres[j])
            touching = apart <= radius[i] + radius[j] + REACH
            both = charged[i] and charged[j] and apart <= COULOMB_REACH
            if not touching and not both:
                continue
            table = (
                contact(batches[i], batches[j])
                if touching
                else np.zeros((len(sites[i].states), len(sites[j].states), len(TERMS)))
            )
            if both:
                table[..., COULOMB] = coulomb(batches[i], batches[j])
            if np.abs(table).max() > 1e-9:
                pairs[(i, j)] = table
    return Tables(sites, single, pairs)
