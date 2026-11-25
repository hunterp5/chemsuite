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

"""Missing loops: the residues of a gap in the deposited sequence built onto both its ends
(docs/LOOPS_PLAN.md).

A gap of *m* residues between residue *a* (whose C is open) and residue *b* (whose N is
open) is closed by cyclic coordinate descent (Canutescu & Dunbrack, Protein Sci. 2003, 12,
963) on Ramachandran-sampled backbones:

1. The free torsions are ψ of *a*, φ and ψ of each new residue, and φ of *b*; every ω is
   trans but a share of those before a Pro (:data:`CIS_SHARE`), which are cis, and bonds and
   angles are Engh & Huber's (Acta Cryst. A 1991, 47, 392), *b*'s copy aside, which takes
   *b*'s own. *a*'s O is rebuilt in its new peptide plane.
2. :data:`CANDIDATES` backbones are built forward from *a*, each new residue's φ/ψ drawn from
   its class's grid (general, Gly, Pro, pre-Pro, after Lovell et al., Proteins 2003, 50, 437;
   ``data/ramachandran.json``, from ``scripts/build_ramachandran.py``).
3. Each is closed by CCD onto *b*'s N, CA and C, every torsion turned by the closed-form best
   angle or a share of it (:data:`STEPS`), the largest that leaves its residue's φ/ψ allowed
   (or no worse). Closures under :data:`CLOSED` Å RMSD over the three atoms are kept.
4. They are scored: heavy-atom overlaps (backbone and a CB) with the rest of the structure
   and within the loop, the Ramachandran log-probability measured on the closed loop, and a
   penalty for each outlier.
5. For :data:`ROUNDS` rounds, the best few are nudged at two torsions and closed again, and
   the new closures join the pool. While no loop in it is clear of overlaps (under
   :data:`GOOD`) with no outlier, another pool is drawn, up to :data:`POOLS`.
6. The best few that differ from each other get their side chains
   (:func:`.repair.build_missing`, turned clear) and are scored again with them; the best is
   kept.

The random draws are seeded from the gap's residue labels, so preparing again gives the same
loop. A gap whose ends are too far apart for its residues, or that no candidate closes, is left
to the caps. Each model of a multi-model file is closed on its own anchors.
"""

from __future__ import annotations

import json
import zlib
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

from ..structure import edit
from ..structure.model import AMINO, Structure
from .hydrogens import HYDROGEN
from .inspect import ends
from .repair import _free_icode, build_missing
from .templates import Library, place

#: the longest gap built (residues): longer ones are capped
LONGEST = 12
CANDIDATES = 300
ROUNDS = 5  # of refinement: the best closed loops nudged and closed again
PARENTS, CHILDREN = 10, 30  # loops nudged per round, and nudged copies of each
NUDGE = 20.0  # degrees, the spread of a nudge to one torsion
POOLS = 3  # draws and refinements, the next only while none holds a GOOD loop
GOOD = 1.5  # Å² of backbone overlap under which a loop with no outlier is good enough
SWEEPS = 150  # CCD passes over every torsion, at most
STALL = 10  # passes with no new closure (after the first) after which CCD stops
NUDGE_SWEEPS = 40  # for a nudged loop, which starts closed
STEPS = (1.0, 0.5, 0.25, 0.1)  # the shares of a CCD turn tried, largest first
#: Å, RMSD over b's N, CA and C at which a candidate counts as closed. Ideal geometry with
#: trans peptides cannot always do better on a short, tight gap (a two-residue helical turn
#: rebuilt from the deposit's own torsions stops at 0.13 Å); Minimize takes up the rest.
CLOSED = 0.25
TIGHT = 0.1  # Å: CCD goes on turning a candidate until it is this close (or stalls)
CLOSURE_WEIGHT = 100.0  # per Å² of closure RMSD, so a tighter closure scores better
CIS_SHARE = 0.15  # of the candidates drawn with a cis peptide before each Pro
CIS_COST = 3.0  # per cis peptide: about -ln(0.05), cis-Pro's share of Pro peptides
RESCORED = 6  # the best closed backbones whose side chains are built and scored
DISTINCT = 0.5  # Å CA RMSD by which each of those differs from the others
# Engh & Huber 1991: Å and degrees
N_CA, CA_C, C_N, C_O = 1.458, 1.525, 1.329, 1.231
N_CA_C, CA_C_N, C_N_CA, CA_C_O = 111.2, 116.2, 121.7, 120.1
CB_BOND, N_CA_CB, C_N_CA_CB = 1.53, 110.5, -123.0  # the pseudo-CB of a non-Gly residue
REACH = 3.8  # Å, the most one residue adds to the CA–CA span
CLASH = 3.0  # Å: heavy atoms closer than this overlap
OVERLAP_WEIGHT = 10.0  # per Å² of overlap, against the summed log-probability
OUTLIER = 100.0  # per Ramachandran outlier, measured on the closed loop
DATA = Path(__file__).parent / "data" / "ramachandran.json"
GENERAL, GLY, PRO, PREPRO = "general", "gly", "pro", "prepro"


@dataclass(frozen=True)
class Gap:
    before: int  # residue a (its C open)
    after: int  # residue b (its N open)
    names: tuple[str, ...]  # the missing residues, from the deposited sequence
    first: int  # the first one's place in the entity's sequence


@dataclass(frozen=True)
class Grid:
    logp: np.ndarray  # (bins, bins) natural-log probability over φ (rows) and ψ, from -180°
    step: float  # degrees
    allowed: float  # the level under which a φ/ψ is an outlier

    def at(self, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
        """Log-probability at φ/ψ in radians (any shape)."""
        bins = self.logp.shape[0]
        i = (np.floor((np.degrees(phi) + 180.0) / self.step).astype(np.int64)) % bins
        j = (np.floor((np.degrees(psi) + 180.0) / self.step).astype(np.int64)) % bins
        return self.logp[i, j]


@cache
def grids() -> dict[str, Grid]:
    classes = json.loads(DATA.read_text(encoding="utf-8"))["classes"]
    return {name: Grid(np.array(g["logp"]), float(g["step"]), float(g["allowed"]))
            for name, g in classes.items()}  # fmt: skip


def rama_class(name: str, following: str | None) -> str:
    if name == "GLY":
        return GLY
    if name == "PRO":
        return PRO
    return PREPRO if following == "PRO" else GENERAL


def gaps(s: Structure, longest: int = LONGEST) -> list[Gap]:
    """The gaps of 2 to *longest* residues inside a chain, between two amino acids with the
    backbone to build on (a without an OXT)."""
    found = ends(s)
    out = []
    for a in sorted(found.open_c):
        b = a + 1
        if b not in found.open_n or s.chain_of[a] != s.chain_of[b]:
            continue
        e = int(s.entity[a])
        if s.seqres[a] < 0 or s.entity[b] != e or not 0 <= e < len(s.entity_id):
            continue
        m = int(s.seqres[b] - s.seqres[a] - 1)
        if not 2 <= m <= longest:
            continue
        if any(n not in _heavy(s, a) for n in ("N", "CA", "C")) or "OXT" in _heavy(s, a):
            continue
        if any(n not in _heavy(s, b) for n in ("N", "CA", "C")):
            continue
        sequence = s.sequence_of(e)
        names = tuple(str(sequence[k]) for k in range(int(s.seqres[a]) + 1, int(s.seqres[b])))
        out.append(Gap(a, b, names, int(s.seqres[a]) + 1))
    return out


def _heavy(s: Structure, r: int) -> dict[str, int]:
    span = s.atoms_of(r)
    return {str(s.name[i]): i for i in range(span.start, span.stop) if s.element[i] != HYDROGEN}


def fill_gaps(
    s: Structure,
    library: Library,
    longest: int = LONGEST,
    cancel: CancelToken | None = None,
) -> tuple[Structure, list[str]]:
    """*s* with each gap of 2 to *longest* residues built (side chains too), where one closes.
    Returns the labels of the residues built."""
    labels = [(s.residue_label(g.before), s.residue_label(g.after)) for g in gaps(s, longest)]
    built: list[str] = []
    for before, after in labels:
        check(cancel)
        where = {s.residue_label(r): r for r in range(s.n_residues)}
        gap = next(g for g in gaps(s, longest) if g.before == where[before])
        done = _fill(s, gap, library, f"{before}|{after}", cancel)
        if done is not None:
            s, new = done
            built += new
    return s, built


def _fill(
    s: Structure, gap: Gap, library: Library, seed: str, cancel: CancelToken | None
) -> tuple[Structure, list[str]] | None:
    per_model = []
    for k in range(s.n_models):
        check(cancel)
        rng = np.random.default_rng([zlib.crc32(seed.encode()), k])
        options = _close(s, gap, k, rng)
        if not options:
            return None
        per_model.append(options)
    # side chains of model 0's best few decide; the other models take their best backbone
    best = None
    for backbone, score in per_model[0]:
        xyz = np.stack([backbone] + [options[0][0] for options in per_model[1:]])
        trial, new = _insert(s, gap, xyz, library, cancel)
        total = score + OVERLAP_WEIGHT * _side_chain_overlap(trial, new)
        if best is None or total < best[0]:
            best = (total, trial, new)
    if best is None:
        return None
    _total, out, new = best
    return out, [out.residue_label(r) for r in new]


# -- closing one gap ---------------------------------------------------------------------------
def _close(s: Structure, gap: Gap, k: int, rng: np.random.Generator) -> list[tuple]:
    """The :data:`RESCORED` best distinct closed backbones for *gap* in model *k*: each a
    ((m + 1, 4, 3) array of N, CA, C, O — a's rebuilt O first, then the new residues; score)
    pair. Empty when no candidate closes."""
    x = s.coords[k].astype(np.float64)
    a, b = _heavy(s, gap.before), _heavy(s, gap.after)
    m = len(gap.names)
    fixed_start = np.stack([x[a["N"]], x[a["CA"]], x[a["C"]]])
    target = np.stack([x[b["N"]], x[b["CA"]], x[b["C"]]])
    if np.linalg.norm(x[a["CA"]] - x[b["CA"]]) > REACH * (m + 1):
        return []
    names = [str(s.res_name[gap.before]), *gap.names, str(s.res_name[gap.after])]
    after_b = _next_name(s, gap.after)
    classes = [rama_class(n, f) for n, f in zip(names, [*names[1:], after_b], strict=True)]
    g = grids()
    # the fixed halves of the anchors' φ/ψ (None where a neighbour is missing)
    phi_a = _phi(s, x, gap.before)
    psi_b = _psi(s, x, gap.after)

    # the peptides before the new residues and b: a Pro's may be cis
    pro = np.array([n == "PRO" for n in names[1:]])

    def attempt(angles: np.ndarray, omega: np.ndarray, sweeps: int) -> tuple[np.ndarray, ...]:
        """The torsions, ω, backbones, scores and goodness of those of *angles* that close."""
        chain = _build(fixed_start, angles, target, omega)  # (K, 3(m + 2), 3)
        chain, angles, rmsd = _ccd(chain, angles, target, classes, phi_a, psi_b, g, sweeps)
        closed = rmsd < CLOSED
        chain, angles, omega, rmsd = chain[closed], angles[closed], omega[closed], rmsd[closed]
        chain[:, -3:] = target  # b as it is: the angles at the junction are measured to it
        rama, outliers = _rama(_torsions(chain), classes, phi_a, psi_b, g)
        backbone = _with_oxygens(chain, x[b["N"]], m)  # (K, m + 1, 4, 3)
        overlap = _overlap(s, gap, backbone, names, x)
        score = OVERLAP_WEIGHT * overlap - rama + OUTLIER * outliers + CLOSURE_WEIGHT * rmsd**2
        score += CIS_COST * np.sum(np.cos(omega) > 0, axis=1)
        return angles, omega, backbone, score, (overlap < GOOD) & (outliers == 0)

    def pool() -> tuple[np.ndarray, ...]:
        """Closed loops from :data:`CANDIDATES` draws, refined for :data:`ROUNDS` rounds."""
        torsions = _draw(rng, [g[c] for c in classes[1:-1]], CANDIDATES)  # (K, 2m) φ1 ψ1 …
        psi_a0 = _psi(s, x, gap.before)
        psi_a = np.full((CANDIDATES, 1), np.pi if psi_a0 is None else psi_a0)
        psi_a[1:] += rng.normal(0.0, np.radians(30.0), (CANDIDATES - 1, 1))
        phi_b = rng.uniform(-np.pi, np.pi, (CANDIDATES, 1))
        cis = pro & (rng.random((CANDIDATES, len(pro))) < CIS_SHARE)
        omega = np.where(cis, 0.0, np.pi)
        # ψa, φ1, ψ1, …, φm, ψm, φb
        found = attempt(np.concatenate([psi_a, torsions, phi_b], axis=1), omega, SWEEPS)
        for _round in range(ROUNDS):
            angles, omega, score = found[0], found[1], found[3]
            if not len(score):
                break
            # the best so far, each nudged at two torsions and closed again
            parents = np.argsort(score, kind="stable")[:PARENTS]
            children = np.repeat(angles[parents], CHILDREN, axis=0)
            rows = np.arange(len(children))
            for _turn in range(2):
                columns = rng.integers(children.shape[1], size=len(children))
                children[rows, columns] += rng.normal(0.0, np.radians(NUDGE), len(children))
            more = attempt(children, np.repeat(omega[parents], CHILDREN, axis=0), NUDGE_SWEEPS)
            found = tuple(np.concatenate([old, new]) for old, new in zip(found, more, strict=True))
        return found

    # pools until one holds a loop that is clear and allowed, or POOLS of them
    found = pool()
    for _pool in range(POOLS - 1):
        if found[4].any():
            break
        found = tuple(np.concatenate([old, new]) for old, new in zip(found, pool(), strict=True))
    _angles, _omega, backbone, score, _good = found
    # the best, then each next best whose CAs are DISTINCT from every one kept
    kept: list[int] = []
    for i in np.argsort(score, kind="stable"):
        ca = backbone[i, 1:, 1]
        if all(_rmsd(ca[None], backbone[j, 1:, 1])[0] > DISTINCT for j in kept):
            kept.append(int(i))
            if len(kept) == RESCORED:
                break
    return [(backbone[i], float(score[i])) for i in kept]


def _next_name(s: Structure, r: int) -> str | None:
    nxt = r + 1
    if nxt < s.n_residues and s.chain_of[nxt] == s.chain_of[r] and s.res_kind[nxt] == AMINO:
        return str(s.res_name[nxt])
    return None


def _phi(s: Structure, x: np.ndarray, r: int) -> float | None:
    prev = r - 1
    if prev < 0 or s.chain_of[prev] != s.chain_of[r] or "C" not in _heavy(s, prev):
        return None
    here = _heavy(s, r)
    c_prev = x[_heavy(s, prev)["C"]]
    if np.linalg.norm(c_prev - x[here["N"]]) > 2.0:
        return None
    return float(_dihedral(c_prev, x[here["N"]], x[here["CA"]], x[here["C"]]))


def _psi(s: Structure, x: np.ndarray, r: int) -> float | None:
    """ψ from the next residue's N, or, before a gap, from this residue's own O."""
    here = _heavy(s, r)
    nxt = r + 1
    if nxt < s.n_residues and s.chain_of[nxt] == s.chain_of[r] and "N" in _heavy(s, nxt):
        n_next = x[_heavy(s, nxt)["N"]]
        if np.linalg.norm(n_next - x[here["C"]]) < 2.0:
            return float(_dihedral(x[here["N"]], x[here["CA"]], x[here["C"]], n_next))
    if "O" in here:
        o = _dihedral(x[here["N"]], x[here["CA"]], x[here["C"]], x[here["O"]])
        return float(np.angle(np.exp(1j * (o + np.pi))))
    return None


def _dihedral(p0, p1, p2, p3) -> np.ndarray:
    """Dihedral angles (radians) of points with any leading shape."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    return np.arctan2(np.sum(_cross(b1, v) * w, axis=-1), np.sum(v * w, axis=-1))


def _draw(rng: np.random.Generator, classes: list[Grid], count: int) -> np.ndarray:
    """(count, 2 × len(classes)) φ/ψ pairs in radians, each drawn from its class's grid."""
    out = np.empty((count, 2 * len(classes)))
    for i, grid in enumerate(classes):
        p = np.exp(grid.logp).ravel()
        cell = rng.choice(p.size, size=count, p=p / p.sum())
        bins = grid.logp.shape[0]
        phi = -180.0 + grid.step * (cell // bins + rng.random(count))
        psi = -180.0 + grid.step * (cell % bins + rng.random(count))
        out[:, 2 * i], out[:, 2 * i + 1] = np.radians(phi), np.radians(psi)
    return out


def _build(
    start: np.ndarray,
    angles: np.ndarray,
    end: np.ndarray | None = None,
    omega: np.ndarray | None = None,
) -> np.ndarray:
    """N, CA, C of a and every residue after it, b's copy last (K, 3(m + 2), 3), from a's
    N, CA, C (*start*) and the torsions (K, 2m + 2: ψa, φ1, ψ1, …, φm, ψm, φb). With *end*
    (b's N, CA, C), b's copy takes b's own N–CA and CA–C lengths and N–CA–C angle, so that
    it can land on b exactly. *omega* (K, m + 1) are the peptides' ω, trans without it."""
    k, t = angles.shape
    residues = t // 2 + 1  # a, the m new, b
    chain = np.empty((k, 3 * residues, 3))
    chain[:, :3] = start
    for r in range(1, residues):
        n_prev, ca_prev, c_prev = chain[:, 3 * r - 3], chain[:, 3 * r - 2], chain[:, 3 * r - 1]
        psi_prev, phi = angles[:, 2 * r - 2], angles[:, 2 * r - 1]
        n_ca, ca_c, n_ca_c = N_CA, CA_C, np.radians(N_CA_C)
        if end is not None and r == residues - 1:
            n_ca, ca_c = np.linalg.norm(end[1] - end[0]), np.linalg.norm(end[2] - end[1])
            n_ca_c = _angle(end[0], end[1], end[2])
        n = place(n_prev, ca_prev, c_prev, C_N, np.radians(CA_C_N), psi_prev)
        w = np.pi if omega is None else omega[:, r - 1]
        ca = place(ca_prev, c_prev, n, n_ca, np.radians(C_N_CA), w)
        c = place(c_prev, n, ca, ca_c, n_ca_c, phi)
        chain[:, 3 * r], chain[:, 3 * r + 1], chain[:, 3 * r + 2] = n, ca, c
    return chain


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """The angle a–b–c in radians."""
    u, w = a - b, c - b
    return float(np.arccos(np.clip(np.dot(u, w) / np.linalg.norm(u) / np.linalg.norm(w), -1, 1)))


def _ccd(chain, angles, target, classes, phi_a, psi_b, g, sweeps: int = SWEEPS):
    """CCD on every candidate at once (those still open), each torsion's turn cut back where
    the whole would take its residue's φ/ψ out of the allowed region (:func:`_allowed_step`).
    Stops after *sweeps* passes, or :data:`STALL` passes with no new closure."""
    t = angles.shape[1]
    chain, angles = chain.copy(), angles.copy()
    rmsd = _rmsd(chain[:, -3:], target)
    closed, since = 0, 0
    for _sweep in range(sweeps):
        open_ = np.flatnonzero(rmsd >= TIGHT)
        now = len(rmsd) - len(open_)
        since = 0 if now > closed or not now else since + 1
        closed = now
        if not len(open_) or since >= STALL:
            break
        sub, turns = chain[open_], angles[open_]
        for j in range(t):
            residue = (j + 1) // 2
            # ψ of residue r turns about CA→C and moves from the next N on; φ about N→CA
            # and moves from its C on
            if j % 2 == 0:
                axis_from, axis_to, moving = 3 * residue + 1, 3 * residue + 2, 3 * residue + 3
            else:
                axis_from, axis_to, moving = 3 * residue, 3 * residue + 1, 3 * residue + 2
            origin = sub[:, axis_to]
            axis = origin - sub[:, axis_from]
            axis /= np.linalg.norm(axis, axis=-1, keepdims=True)
            alpha = _best_turn(sub[:, -3:], target, origin, axis)
            alpha = _allowed_step(turns, j, alpha, classes, phi_a, psi_b, g)
            sub[:, moving:] = _rotate(sub[:, moving:], origin, axis, alpha)
            turns[:, j] = np.angle(np.exp(1j * (turns[:, j] + alpha)))
        chain[open_], angles[open_] = sub, turns
        rmsd[open_] = _rmsd(sub[:, -3:], target)
    return chain, angles, rmsd


def _cross(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """u × v over the last axis, broadcast (np.cross without its overhead on small arrays)."""
    u, v = np.broadcast_arrays(u, v)
    return np.stack([u[..., 1] * v[..., 2] - u[..., 2] * v[..., 1],
                     u[..., 2] * v[..., 0] - u[..., 0] * v[..., 2],
                     u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]], axis=-1)  # fmt: skip


def _best_turn(moving, fixed, origin, axis) -> np.ndarray:
    """The turn about *axis* through *origin* (per candidate) that brings the *moving* points
    (K, 3, 3) closest to *fixed* (3, 3) (Canutescu & Dunbrack's closed form)."""
    u = axis[:, None, :]
    r = moving - origin[:, None, :]
    r = r - np.sum(r * u, axis=-1, keepdims=True) * u
    size = np.linalg.norm(r, axis=-1, keepdims=True)
    r_hat = r / np.maximum(size, 1e-9)
    s_hat = _cross(u, r_hat)
    f = fixed[None] - origin[:, None, :]
    a = np.sum(size[..., 0] * np.sum(f * r_hat, axis=-1), axis=1)
    b = np.sum(size[..., 0] * np.sum(f * s_hat, axis=-1), axis=1)
    return np.arctan2(b, a)


def _allowed_step(angles, j, alpha, classes, phi_a, psi_b, g) -> np.ndarray:
    """*alpha*, or the largest of its :data:`STEPS` that keeps torsion *j*'s residue in its
    class's allowed region, or no worse than it is when it is outside (0 when none does). A
    residue whose other angle is unknown is not held."""
    residue = (j + 1) // 2
    grid = g[classes[residue]]
    t = angles.shape[1]
    if j % 2 == 0:  # ψ: its φ is the torsion before, or a's own
        other = angles[:, j - 1] if j > 0 else phi_a
        pair = (lambda v, o=other: (o, v)) if other is not None else None
    else:  # φ: its ψ is the torsion after, or b's own
        other = angles[:, j + 1] if j + 1 < t else psi_b
        pair = (lambda v, o=other: (v, o)) if other is not None else None
    if pair is None:
        return alpha
    chosen = np.zeros_like(alpha)
    left = np.ones(len(alpha), dtype=bool)
    # an outlier may move as long as it gets no worse
    floor = np.minimum(grid.at(*pair(angles[:, j])), grid.allowed)
    for scale in STEPS:
        trial = angles[:, j] + scale * alpha
        ok = left & (grid.at(*pair(trial)) >= floor)
        chosen[ok] = scale * alpha[ok]
        left &= ~ok
    return chosen


def _rotate(p: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: np.ndarray) -> np.ndarray:
    """Points *p* (K, n, 3) turned by *angle* (K,) about the line through *origin* (K, 3)
    along the unit *axis* (K, 3), per candidate (Rodrigues)."""
    o, u = origin[:, None, :], axis[:, None, :]
    cos, sin = np.cos(angle)[:, None, None], np.sin(angle)[:, None, None]
    v = p - o
    dot = np.sum(v * u, axis=-1, keepdims=True)
    return o + v * cos + _cross(u, v) * sin + u * dot * (1 - cos)


def _rmsd(moving: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean(np.sum((moving - fixed[None]) ** 2, axis=-1), axis=-1))


def _torsions(chain: np.ndarray) -> np.ndarray:
    """The torsions of *chain* (K, 3(m + 2), 3) as :func:`_build` takes them: ψa, φ1, ψ1, …,
    φm, ψm, φb."""
    n, ca, c = chain[:, 0::3], chain[:, 1::3], chain[:, 2::3]
    psi = _dihedral(n[:, :-1], ca[:, :-1], c[:, :-1], n[:, 1:])  # residues 0 … m
    phi = _dihedral(c[:, :-1], n[:, 1:], ca[:, 1:], c[:, 1:])  # residues 1 … m + 1
    out = np.empty((chain.shape[0], 2 * psi.shape[1]))
    out[:, 0::2], out[:, 1::2] = psi, phi
    return out


def _rama(angles, classes, phi_a, psi_b, g) -> tuple[np.ndarray, np.ndarray]:
    """Per candidate, the summed log-probability of every residue's φ/ψ whose two angles are
    known, and how many of those are outliers."""
    t = angles.shape[1]
    total = np.zeros(angles.shape[0])
    outliers = np.zeros(angles.shape[0], dtype=np.int64)
    for r, name in enumerate(classes):
        phi = angles[:, 2 * r - 1] if r > 0 else phi_a
        psi = angles[:, 2 * r] if 2 * r < t else psi_b
        if phi is None or psi is None:
            continue
        p = g[name].at(np.broadcast_to(phi, total.shape), np.broadcast_to(psi, total.shape))
        total += p
        outliers += p < g[name].allowed
    return total, outliers


def _with_oxygens(chain: np.ndarray, n_b: np.ndarray, m: int) -> np.ndarray:
    """(K, m + 1, 4, 3): a's rebuilt O alone in row 0 (N, CA, C as they are), then each new
    residue's N, CA, C, O, every O in the plane of its peptide (b's real N after the last)."""
    k = chain.shape[0]
    out = np.empty((k, m + 1, 4, 3))
    for r in range(m + 1):
        n, ca, c = chain[:, 3 * r], chain[:, 3 * r + 1], chain[:, 3 * r + 2]
        n_next = chain[:, 3 * r + 3] if r < m else np.broadcast_to(n_b, n.shape)
        psi = _dihedral(n, ca, c, n_next)
        o = place(n, ca, c, C_O, np.radians(CA_C_O), psi + np.pi)
        out[:, r] = np.stack([n, ca, c, o], axis=1)
    return out


def _overlap(s: Structure, gap: Gap, backbone, names, x) -> np.ndarray:
    """Per candidate, the summed squared overlap (Å²) of the new atoms (backbone, a pseudo-CB,
    a's rebuilt O) with every other heavy atom (waters, and the anchors' side chains beyond
    CB, too), and with each other two or more residues apart (the anchors' backbone and CB
    included)."""
    count, m1 = backbone.shape[0], backbone.shape[1]
    points, residue = [], []
    for r in range(m1):
        n, ca, c, o = (backbone[:, r, i] for i in range(4))
        if r == 0:
            points.append(o[:, None])
            residue.append([0])
            continue
        atoms = [n, ca, c, o]
        if names[r] != "GLY":
            atoms.append(place(c, n, ca, CB_BOND, np.radians(N_CA_CB), np.radians(C_N_CA_CB)))
        points.append(np.stack(atoms, axis=1))
        residue.append([r] * len(atoms))
    new = np.concatenate(points, axis=1)  # (K, P, 3)
    owner = np.concatenate(residue)
    # the rest of the structure, the anchors' side chains beyond CB included
    near = ("N", "CA", "C", "O", "CB")
    keep = (s.element != HYDROGEN) & ~(np.isin(s.residue_of, [gap.before, gap.after])
                                       & np.isin(s.name, near))  # fmt: skip
    total = np.zeros(count)
    if keep.any() and count:
        tree = cKDTree(x[keep])
        d, _ = tree.query(new.reshape(-1, 3), k=8, distance_upper_bound=CLASH)
        over = np.clip(CLASH - d, 0.0, None)
        over[~np.isfinite(d)] = 0.0
        total += np.sum(over**2, axis=1).reshape(count, -1).sum(axis=1)
    # within the loop, with the anchors' backbone and CB (a's O aside) as residues 0 and m + 1
    a, b = _heavy(s, gap.before), _heavy(s, gap.after)
    fixed = [a[n] for n in ("N", "CA", "C", "CB") if n in a]
    fixed_owner = [0] * len(fixed) + [m1] * len([n for n in near if n in b])
    fixed += [b[n] for n in near if n in b]
    inner = np.concatenate([new, np.broadcast_to(x[fixed], (count, len(fixed), 3))], axis=1)
    owner = np.concatenate([owner, fixed_owner])
    apart = np.triu(np.abs(owner[:, None] - owner[None, :]) >= 2, 1)
    d = np.linalg.norm(inner[:, :, None] - inner[:, None, :], axis=-1)
    over = np.clip(CLASH - d, 0.0, None) * apart
    return total + np.sum(over**2, axis=(1, 2))


# -- into the structure ------------------------------------------------------------------------
def _insert(
    s: Structure, gap: Gap, xyz: np.ndarray, library: Library, cancel: CancelToken | None
) -> tuple[Structure, list[int]]:
    """*s* with a's O moved and the gap's residues put in (backbone from *xyz*, (models,
    m + 1, 4, 3)), then their side chains built. Returns it and the new residues' indices."""
    m = len(gap.names)
    o = _heavy(s, gap.before).get("O")
    coords = s.coords.copy()
    if o is not None:
        coords[:, o] = xyz[:, 0, 3]
    s = s.with_coords(coords)
    chain = int(s.chain_of[gap.before])
    first_seq = int(s.res_seq[gap.before]) + 1
    fits = int(s.res_seq[gap.after]) - int(s.res_seq[gap.before]) - 1 >= m
    new = []
    for i in range(m):
        at = gap.before + 1 + i  # b's place, as each residue goes in before it
        n = s.n_atoms
        c_prev = _heavy(s, at - 1)["C"]
        bonds = [[c_prev, n], [n, n + 1], [n + 1, n + 2], [n + 2, n + 3]]
        order = [1, 1, 1, 2]
        if i == m - 1:
            bonds.append([n + 2, _heavy(s, at)["N"]])
            order.append(1)
        seq = first_seq + i if fits else int(s.res_seq[gap.before])
        s = edit.insert_residue(s, at, chain, gap.names[i], seq, np.array(["N", "CA", "C", "O"]),
                                np.array([7, 6, 6, 8]), xyz[:, i + 1], np.array(bonds),
                                np.array(order), icode="" if fits else _free_icode(s, chain, seq),
                                seqres=gap.first + i, entity=int(s.entity[gap.before]))  # fmt: skip
        new.append(at)
    s = build_missing(s, library, cancel, only=frozenset(new))
    return s, new


def _side_chain_overlap(s: Structure, new: list[int]) -> float:
    """The summed squared overlap (Å², first model) of the side-chain atoms of the residues
    *new* with heavy atoms of residues that are not theirs or next to them in the chain."""
    x = s.coords[0].astype(np.float64)
    heavy = s.element != HYDROGEN
    mine = np.isin(s.residue_of, new) & heavy
    mine &= ~np.isin(s.name, ["N", "CA", "C", "O"])
    if not mine.any():
        return 0.0
    tree = cKDTree(x[heavy])
    owner = s.residue_of[heavy]
    total = 0.0
    for i in np.flatnonzero(mine):
        r = int(s.residue_of[i])
        for j in tree.query_ball_point(x[i], CLASH):
            if abs(int(owner[j]) - r) >= 2 or s.chain_of[owner[j]] != s.chain_of[r]:
                d = float(np.linalg.norm(tree.data[j] - x[i]))
                total += (CLASH - d) ** 2
    return total
