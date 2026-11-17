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
"""Field superposition: electrostatic potential overlap plus shape.

Each molecule's potential is a sum over atoms of ``q / r``, with ``1/r`` approximated by three
Gaussians (as Good, Hodgkin and Richards, J. Chem. Inf. Comput. Sci. 1992, 32, 188, did; the
coefficients below are our own least-squares fit over 1–8 Å, within 3.3% there). The overlap
of two such potentials is then an analytic sum of Gaussian integrals, and the Carbó index
``F_AB / sqrt(F_AA F_BB)`` (−1…1) is the field similarity. Shape is the first-order Gaussian
overlap volume (Grant and Pickup, J. Phys. Chem. 1995, 99, 3503) as a Tanimoto.

The rigid pose maximizes ``(1 − w) · shape + w · (field + 1) / 2``. It starts from the best
shape poses of the best-fitting conformers and from those poses turned 180° about each of the
probe's principal axes (where shape alone cannot tell, e.g. a pyridine ring flipped), then
refines with L-BFGS. Hydrogens' charges are added to the atoms they are bonded to.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdPartialCharges
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties
from scipy.optimize import minimize

from . import shape
from .model import FIELD_SIMILARITY, SHAPE, SHAPE_TANIMOTO, Fit, Settings, apply_transform

#: 1/r ≈ Σ a_k exp(−b_k r²), r in Å (fitted over 1–8 Å, relative error under 3.3%)
POTENTIAL_A = np.array([0.84916, 0.36174, 0.23143])
POTENTIAL_B = np.array([0.70954, 0.12055, 0.0098609])
#: Gaussian height for the shape volume (as ROCS and PubChem)
SHAPE_P = 2.7
#: conformers (best by shape first) that are refined by the field
REFINED_CONFORMERS = 3

# ∫ a_k e^{−b_k |r−A|²} a_l e^{−b_l |r−B|²} dr = c_kl e^{−e_kl d²}
_SUM = POTENTIAL_B[:, None] + POTENTIAL_B[None, :]
_FIELD_C = (POTENTIAL_A[:, None] * POTENTIAL_A[None, :] * (np.pi / _SUM) ** 1.5).ravel()
_FIELD_E = (POTENTIAL_B[:, None] * POTENTIAL_B[None, :] / _SUM).ravel()


class Model:
    """A molecule as heavy atoms: positions, charges (hydrogens' included) and shape widths."""

    def __init__(self, mol: Chem.Mol, conf_id: int, charges: str) -> None:
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        per_atom = atom_charges(mol, charges)
        merged = per_atom.copy()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                for neighbor in atom.GetNeighbors():
                    merged[neighbor.GetIdx()] += per_atom[atom.GetIdx()]
                    break
        table = Chem.GetPeriodicTable()
        radii = np.array([table.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in heavy])
        self.atoms = heavy
        self.xyz = mol.GetConformer(conf_id).GetPositions()[heavy]
        self.charges = merged[heavy]
        self.alpha = np.pi * (3.0 * SHAPE_P / (4.0 * np.pi)) ** (2.0 / 3.0) / radii**2


def atom_charges(mol: Chem.Mol, kind: str) -> np.ndarray:
    """Partial charges per atom: MMFF94 (where it has parameters) or Gasteiger."""
    with rdBase.BlockLogs():
        if kind == "mmff":
            props = MMFFGetMoleculeProperties(mol)
            if props is not None:
                return np.array([props.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms())])
        work = Chem.Mol(mol)
        rdPartialCharges.ComputeGasteigerCharges(work)
        values = [a.GetDoubleProp("_GasteigerCharge") for a in work.GetAtoms()]
    return np.nan_to_num(np.array(values, dtype=np.float64))


def _distances2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)


def shape_overlap(a: Model, a_xyz: np.ndarray, b: Model, b_xyz: np.ndarray) -> float:
    sum_ = a.alpha[:, None] + b.alpha[None, :]
    k = (SHAPE_P * SHAPE_P) * (np.pi / sum_) ** 1.5
    e = a.alpha[:, None] * b.alpha[None, :] / sum_
    return float((k * np.exp(-e * _distances2(a_xyz, b_xyz))).sum())


def field_overlap(a: Model, a_xyz: np.ndarray, b: Model, b_xyz: np.ndarray) -> float:
    d2 = _distances2(a_xyz, b_xyz)[..., None]
    kernel = (_FIELD_C * np.exp(-_FIELD_E * d2)).sum(axis=2)
    return float(a.charges @ kernel @ b.charges)


def _rotation(vector: np.ndarray) -> np.ndarray:
    """The rotation matrix of a rotation vector (Rodrigues)."""
    angle = float(np.linalg.norm(vector))
    if angle < 1e-12:
        return np.eye(3)
    k = vector / angle
    cross = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * cross @ cross


class Scorer:
    """The objective for one probe conformer against the reference."""

    def __init__(self, probe: Model, reference: Model, weight: float) -> None:
        self.probe = probe
        self.reference = reference
        self.weight = weight
        self.ref_shape = shape_overlap(reference, reference.xyz, reference, reference.xyz)
        self.ref_field = field_overlap(reference, reference.xyz, reference, reference.xyz)
        self.probe_shape = shape_overlap(probe, probe.xyz, probe, probe.xyz)
        self.probe_field = field_overlap(probe, probe.xyz, probe, probe.xyz)

    def scores(self, xyz: np.ndarray) -> tuple[float, float]:
        """(shape Tanimoto, field Carbó index) of the probe at *xyz*."""
        overlap = shape_overlap(self.probe, xyz, self.reference, self.reference.xyz)
        tanimoto = overlap / (self.ref_shape + self.probe_shape - overlap)
        norm = np.sqrt(self.ref_field * self.probe_field)
        if norm <= 1e-12:  # no charges: the field says nothing
            return tanimoto, 0.0
        carbo = field_overlap(self.probe, xyz, self.reference, self.reference.xyz) / norm
        return tanimoto, carbo

    def objective(self, tanimoto: float, carbo: float) -> float:
        return (1.0 - self.weight) * tanimoto + self.weight * (carbo + 1.0) / 2.0


def _transform(rotation: np.ndarray, center: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """4×4 of ``x → R (x − c) + c + t``."""
    out = np.eye(4)
    out[:3, :3] = rotation
    out[:3, 3] = center + shift - rotation @ center
    return out


def _moved(xyz: np.ndarray, transform: np.ndarray) -> np.ndarray:
    return xyz @ transform[:3, :3].T + transform[:3, 3]


def _starts(xyz: np.ndarray, start: np.ndarray) -> list[np.ndarray]:
    """*start*, then *start* followed by a half turn about each principal axis of the posed
    probe, through its centroid."""
    posed = _moved(xyz, start)
    if len(posed) < 3:
        return [start]
    center = posed.mean(axis=0)
    _values, axes = np.linalg.eigh(np.cov((posed - center).T))
    out = [start]
    for axis in axes.T:
        flip = _transform(_rotation(np.pi * axis), center, np.zeros(3))
        out.append(flip @ start)
    return out


def refine(scorer: Scorer, start: np.ndarray) -> tuple[np.ndarray, float, float, float]:
    """The best pose near *start*: (transform, objective, shape Tanimoto, field Carbó)."""
    xyz = scorer.probe.xyz
    posed = _moved(xyz, start)
    center = posed.mean(axis=0)

    def pose(params: np.ndarray) -> np.ndarray:
        return _transform(_rotation(params[:3]), center, params[3:]) @ start

    def loss(params: np.ndarray) -> float:
        return -scorer.objective(*scorer.scores(_moved(xyz, pose(params))))

    found = minimize(loss, np.zeros(6), method="L-BFGS-B", options={"maxiter": 200})
    transform = pose(found.x)
    tanimoto, carbo = scorer.scores(_moved(xyz, transform))
    return transform, scorer.objective(tanimoto, carbo), tanimoto, carbo


def align(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> list[Fit]:
    """The best field fit of each of the conformers that fit the shape best."""
    by_shape = shape.align(probe, reference, Settings(method=SHAPE))
    by_shape.sort(key=lambda f: -f.scores[SHAPE_TANIMOTO])
    reference_model = Model(reference, 0, settings.charges)
    fits = []
    for start in by_shape[:REFINED_CONFORMERS]:
        model = Model(probe, start.conf_id, settings.charges)
        scorer = Scorer(model, reference_model, settings.field_weight)
        tried = [refine(scorer, pose) for pose in _starts(model.xyz, start.transform)]
        transform, objective, _tanimoto, carbo = max(tried, key=lambda t: t[1])
        work = Chem.Mol(probe)
        apply_transform(work.GetConformer(start.conf_id), transform)
        tanimoto = shape.score(work, reference, start.conf_id)
        scores = {FIELD_SIMILARITY: carbo, SHAPE_TANIMOTO: tanimoto}
        fits.append(Fit(start.conf_id, transform, scores, rank=objective))
    return fits
