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
"""Secondary structure when the file does not give it (a predicted or prepared model).

A reduced DSSP (Kabsch & Sander, 1983) in numpy: backbone hydrogen bonds from the
electrostatic energy of C=O against N-H (the H placed opposite the previous residue's C=O),
helices from two consecutive i -> i+4 turns, strands from bridges (parallel or antiparallel)
with a bridged neighbour. 3-10 and pi helices, bends and turns are left as coil: the cartoon
only draws helices and strands.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .model import AMINO, COIL, HELIX, STRAND, Structure

#: DSSP's energy constant (q1 * q2 * f = 0.42 * 0.20 * 332, kcal/mol * Å) and its
#: hydrogen-bond cut-off (kcal/mol)
_Q = 27.888
_CUTOFF = -0.5
#: CA–CA distance beyond which two residues cannot share a backbone hydrogen bond
_REACH = 9.0
#: longest C(i-1)–N(i) for residues to count as neighbours in the chain
_PEPTIDE = 2.5


def assign(structure: Structure, model: int = 0) -> np.ndarray:
    """Secondary structure (COIL, HELIX, STRAND) of every residue."""
    ss = np.full(structure.n_residues, COIL, dtype=np.uint8)
    backbone = _backbone(structure, model)
    if backbone is None:
        return ss
    residues, n, ca, c, o = backbone
    count = len(residues)
    # residue k follows k - 1 in the chain when the peptide bond is there
    follows = np.zeros(count, dtype=bool)
    follows[1:] = np.linalg.norm(n[1:] - c[:-1], axis=1) < _PEPTIDE
    h = n.copy()
    step = c[:-1] - o[:-1]
    step /= np.maximum(np.linalg.norm(step, axis=1, keepdims=True), 1e-6)
    h[1:][follows[1:]] = (n[1:] + step)[follows[1:]]
    hbond = _hbonds(n, h, ca, c, o, follows)

    def bonded(i: np.ndarray, j: np.ndarray) -> np.ndarray:
        ok = (i >= 0) & (j >= 0) & (i < count) & (j < count)
        out = np.zeros(len(i), dtype=bool)
        out[ok] = [(int(a), int(b)) in hbond for a, b in zip(i[ok], j[ok], strict=True)]
        return out

    index = np.arange(count)
    connected = _runs(follows)
    # alpha helix: 4-turns at i - 1 and i -> residues i..i+3
    turn = bonded(index, index + 4) & _same_run(connected, index, index + 4)
    start = np.zeros(count, dtype=bool)
    start[1:] = turn[1:] & turn[:-1]
    helix = np.zeros(count, dtype=bool)
    for offset in range(4):
        helix[np.minimum(index[start] + offset, count - 1)] = True
    # bridges between residues i and j at least three apart
    bridged = np.zeros(count, dtype=bool)
    if hbond:
        pairs = cKDTree(ca).query_pairs(_REACH, output_type="ndarray")
        pairs = pairs[np.abs(pairs[:, 0] - pairs[:, 1]) > 2]
        i, j = pairs[:, 0], pairs[:, 1]
        parallel = (bonded(i - 1, j) & bonded(j, i + 1)) | (bonded(j - 1, i) & bonded(i, j + 1))
        anti = (bonded(i, j) & bonded(j, i)) | (bonded(i - 1, j + 1) & bonded(j - 1, i + 1))
        hit = parallel | anti
        bridged[i[hit]] = True
        bridged[j[hit]] = True
    # a strand is a bridged residue with a bridged neighbour (a ladder, not a lone bridge)
    neighbour = np.zeros(count, dtype=bool)
    neighbour[1:] |= bridged[:-1] & follows[1:]
    neighbour[:-1] |= bridged[1:] & follows[1:]
    strand = bridged & neighbour & ~helix
    ss[residues[helix]] = HELIX
    ss[residues[strand]] = STRAND
    return ss


def _backbone(
    structure: Structure, model: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Amino-acid residues with all of N, CA, C, O, and those atoms' positions."""
    active = structure.active()
    xyz = structure.xyz(model)
    amino = structure.atom_kind == AMINO
    found = {}
    for name in ("N", "CA", "C", "O"):
        atoms = np.flatnonzero(active & amino & (structure.name == name))
        residues, first = np.unique(structure.residue_of[atoms], return_index=True)
        found[name] = dict(zip(residues.tolist(), atoms[first].tolist(), strict=True))
    residues = sorted(set.intersection(*(set(v) for v in found.values())))
    if len(residues) < 5:
        return None
    pick = {name: xyz[[atoms[r] for r in residues]].astype(np.float64)
            for name, atoms in found.items()}  # fmt: skip
    return np.array(residues), pick["N"], pick["CA"], pick["C"], pick["O"]


def _hbonds(
    n: np.ndarray, h: np.ndarray, ca: np.ndarray, c: np.ndarray, o: np.ndarray, has_h: np.ndarray
) -> set[tuple[int, int]]:
    """``(i, j)`` when the C=O of residue i accepts a hydrogen bond from the N-H of j."""
    pairs = cKDTree(ca).query_pairs(_REACH, output_type="ndarray")
    if not len(pairs):
        return set()
    both = np.concatenate([pairs, pairs[:, ::-1]])
    i, j = both[:, 0], both[:, 1]
    keep = (np.abs(i - j) > 1) & has_h[j]
    i, j = i[keep], j[keep]
    r_on = np.linalg.norm(o[i] - n[j], axis=1)
    r_ch = np.linalg.norm(c[i] - h[j], axis=1)
    r_oh = np.linalg.norm(o[i] - h[j], axis=1)
    r_cn = np.linalg.norm(c[i] - n[j], axis=1)
    with np.errstate(divide="ignore"):
        energy = _Q * (1 / r_on + 1 / r_ch - 1 / r_oh - 1 / r_cn)
    hit = energy < _CUTOFF
    return set(zip(i[hit].tolist(), j[hit].tolist(), strict=True))


def _runs(follows: np.ndarray) -> np.ndarray:
    """A run number per residue: residues in one unbroken stretch of chain share it."""
    return np.cumsum(~follows)


def _same_run(runs: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    ok = j < len(runs)
    out = np.zeros(len(i), dtype=bool)
    out[ok] = runs[i[ok]] == runs[j[ok]]
    return out
