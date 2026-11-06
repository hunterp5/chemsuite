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

"""MolGpKa: per-site pKa from a graph-convolutional network, run with numpy.

A port of the inference path of MolGpKa (Pan et al., J. Chem. Inf. Model. 2021,
doi:10.1021/acs.jcim.1c00075; https://github.com/Xundrug/MolGpKa, MIT). The weights and the
site templates ship in ``data/`` (see ``data/LICENSE-MolGpKa.txt``); the PyTorch weights were
converted with ``scripts/convert_molgpka_weights.py``, so neither torch nor PyTorch Geometric
is needed.

The port reproduces the upstream behaviour exactly, quirks included, because the weights were
trained on it:

* Molecules are uncharged and re-read from SMILES, then hydrogens are made explicit.
* Each site is predicted on its own: the graph is the whole molecule, and two features mark
  the site (its own flag, and every atom's path length to it).
* An acid is predicted at its acidic hydrogen and reported on the heavy atom carrying it.
* Unknown elements count as iodine and unknown hybridisations (including hydrogen's) as
  sp3d2, as upstream's one-hot encoding did.
* The donor and acceptor features are always 0: upstream compared atom indices with a list
  of match tuples, which never matched, so the trained weights expect 0.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from importlib import resources

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

ACID = "acid"
BASE = "base"

_SYMBOLS = ("C", "H", "O", "N", "S", "Cl", "F", "Br", "P", "I")
_HYBRIDISATIONS = (
    Chem.HybridizationType.SP,
    Chem.HybridizationType.SP2,
    Chem.HybridizationType.SP3,
    Chem.HybridizationType.SP3D,
    Chem.HybridizationType.SP3D2,
)
_RING_SIZES = (3, 4, 5, 6, 7, 8)
_N_FEATURES = 29
_BN_EPS = 1e-5  # torch.nn.BatchNorm1d's default


@dataclass(frozen=True)
class Site:
    """One ionisable site of a prepared molecule (see :func:`prepare`)."""

    atom: int  # heavy-atom index in the prepared molecule without hydrogens
    kind: str  # ACID or BASE
    pka: float


# -- templates -------------------------------------------------------------------------------
@cache
def _templates() -> tuple[tuple[Chem.Mol, tuple[int, ...], str], ...]:
    """(pattern, indices of the site atoms in the pattern, ACID/BASE) for each template."""
    text = resources.files(__package__).joinpath("data/molgpka_sites.tsv").read_text("utf-8")
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    out = []
    for row in rows[1:]:
        _number, smarts, index, kind = (cell.strip() for cell in row)
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            raise ValueError(f"bad MolGpKa template: {smarts}")
        out.append(
            (pattern, tuple(int(i) for i in index.split(",")), ACID if kind == "A" else BASE)
        )
    return tuple(out)


def site_atoms(mol_h: Chem.Mol, kind: str) -> list[int]:
    """Atoms (in the molecule with explicit hydrogens) the *kind* model is asked about.

    For acids these are hydrogens; for bases, heavy atoms. Order follows upstream.
    """
    pairs: list[list[int]] = []
    singles: set[int] = set()
    for pattern, index, template_kind in _templates():
        if template_kind != kind:
            continue
        for match in mol_h.GetSubstructMatches(pattern):
            if kind == ACID and len(index) == 2:
                pairs.append([match[index[0]], match[index[1]]])
            else:
                singles.update(match[i] for i in index)
    atoms = [a for pair in pairs for a in pair] + sorted(singles)
    return list(dict.fromkeys(atoms))


# -- the network -----------------------------------------------------------------------------
@cache
def _weights(kind: str) -> dict[str, np.ndarray]:
    path = resources.files(__package__).joinpath(f"data/molgpka_{kind}.npz")
    with path.open("rb") as file, np.load(file) as archive:
        weights = {name: archive[name].astype(np.float32) for name in archive.files}
    for layer in range(1, 6):  # fold eval-mode batch norm into a scale and a shift
        bn = f"bn{layer}."
        scale = weights[bn + "weight"] / np.sqrt(weights[bn + "running_var"] + _BN_EPS)
        weights[bn + "scale"] = scale
        weights[bn + "shift"] = weights[bn + "bias"] - weights[bn + "running_mean"] * scale
    return weights


def _atom_features(mol_h: Chem.Mol) -> np.ndarray:
    """The 27 site-independent features of every atom (upstream column order)."""
    ring = mol_h.GetRingInfo()
    rows = np.zeros((mol_h.GetNumAtoms(), _N_FEATURES), dtype=np.float32)
    for atom in mol_h.GetAtoms():
        i = atom.GetIdx()
        row = rows[i]
        symbol = atom.GetSymbol()
        row[_SYMBOLS.index(symbol) if symbol in _SYMBOLS else len(_SYMBOLS) - 1] = 1
        row[10] = atom.GetDegree()
        hybrid = atom.GetHybridization()
        row[11 + (_HYBRIDISATIONS.index(hybrid) if hybrid in _HYBRIDISATIONS else 4)] = 1
        row[16] = atom.GetValence(Chem.ValenceType.IMPLICIT)
        row[17] = atom.GetIsAromatic()
        for k, size in enumerate(_RING_SIZES):
            row[18 + k] = ring.IsAtomInRingOfSize(i, size)
        # 24, 25: donor / acceptor, always 0 (see the module docstring)
        row[26] = atom.GetFormalCharge()
    return rows


def _propagation(mol_h: Chem.Mol) -> np.ndarray:
    """GCN's symmetric-normalised adjacency with self loops, D^-1/2 (A + I) D^-1/2."""
    adjacency = Chem.GetAdjacencyMatrix(mol_h).astype(np.float32)
    adjacency[np.diag_indices_from(adjacency)] = 1.0
    inv_sqrt = 1.0 / np.sqrt(adjacency.sum(axis=1))
    return adjacency * inv_sqrt[:, None] * inv_sqrt[None, :]


def _forward(kind: str, features: np.ndarray, propagation: np.ndarray) -> np.ndarray:
    """pKa for each site; *features* is (sites, atoms, 29)."""
    w = _weights(kind)
    x = features
    for layer in range(1, 6):
        x = propagation @ (x @ w[f"conv{layer}.weight"]) + w[f"conv{layer}.bias"]
        np.maximum(x, 0.0, out=x)
        x = x * w[f"bn{layer}.scale"] + w[f"bn{layer}.shift"]
    gate = x @ w["att.gate_nn.weight"].T + w["att.gate_nn.bias"]  # (sites, atoms, 1)
    gate = np.exp(gate - gate.max(axis=1, keepdims=True))
    gate /= gate.sum(axis=1, keepdims=True)
    x = (gate * x).sum(axis=1)  # attention pooling over atoms
    x = np.maximum(x @ w["fc2.weight"].T + w["fc2.bias"], 0.0)
    x = np.maximum(x @ w["fc3.weight"].T + w["fc3.bias"], 0.0)
    return (x @ w["fc4.weight"].T + w["fc4.bias"])[:, 0]


def _predict(kind: str, mol_h: Chem.Mol, atoms: list[int], base: np.ndarray) -> np.ndarray:
    if not atoms:
        return np.zeros(0, dtype=np.float32)
    distances = Chem.GetDistanceMatrix(mol_h)
    features = np.repeat(base[None, :, :], len(atoms), axis=0)
    for s, site in enumerate(atoms):
        path = distances[site] + 1.0  # upstream counts the atoms on the path, ends included
        path[distances[site] > 1e6] = 0.0  # another fragment: no path
        path[site] = 0.0
        features[s, :, 27] = path
        features[s, site, 28] = 1.0
    return _forward(kind, features, _propagation(mol_h))


# -- molecules -------------------------------------------------------------------------------
def prepare(mol: Chem.Mol) -> Chem.Mol | None:
    """The molecule MolGpKa sees: largest fragment, uncharged, re-read from SMILES.

    Returns None when RDKit cannot process it.
    """
    try:
        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        return Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    except (RuntimeError, ValueError):
        return None


def predict_sites(prepared: Chem.Mol) -> list[Site]:
    """Every ionisable site of a prepared molecule (see :func:`prepare`) with its pKa.

    Sites are ordered by atom, bases before acids on the same atom. An atom can be both: an
    aniline nitrogen is a weak base and, through its hydrogens, a very weak acid. (Upstream's
    protonation tool keeps only the acid there; the pKa values themselves are the same.) Of
    several equivalent acidic hydrogens on one atom, the last one counts, as upstream.
    """
    mol_h = Chem.AddHs(prepared)
    base_features = _atom_features(mol_h)
    found: dict[tuple[int, str], Site] = {}
    bases = site_atoms(mol_h, BASE)
    for atom, pka in zip(bases, _predict(BASE, mol_h, bases, base_features), strict=True):
        found[atom, BASE] = Site(atom, BASE, float(pka))
    acids = site_atoms(mol_h, ACID)
    for hydrogen, pka in zip(acids, _predict(ACID, mol_h, acids, base_features), strict=True):
        heavy = mol_h.GetAtomWithIdx(hydrogen).GetNeighbors()[0].GetIdx()
        found[heavy, ACID] = Site(heavy, ACID, float(pka))
    # Hydrogens come after the heavy atoms in AddHs, so heavy indices carry over unchanged.
    order = {BASE: 0, ACID: 1}
    return sorted(found.values(), key=lambda site: (site.atom, order[site.kind]))
