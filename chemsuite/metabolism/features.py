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

"""The 43 columns the sites-of-metabolism model reads: 34 per atom, then the isoform.

Ported from ``som-benchmark``'s ``features.py`` (Will de Snoo, MIT; ``basic_atom_block``
and ``charge_block``, read 2026-11-18), which is the benchmark that chose this feature
set over a graph network and over Hückel descriptors — see docs/SOM_PLAN.md. The licence
is in ``data/LICENSE-SOMBENCH.txt``.

Three things here are load-bearing, and all three are "whatever training did":

* **The column order is the model's input layout.** :data:`FEATURE_NAMES` is the
  contract between this module and ``data/som_gbm.npz``; a column inserted anywhere but
  the end silently feeds every tree the wrong number.
* **Heavy atoms only, in canonical SMILES order.** The training labels are positions in
  ``Chem.MolToSmiles`` output order over a molecule with implicit hydrogens, so a
  molecule is featurized exactly as ``Chem.MolFromSmiles(smiles)`` gives it. Hydrogens a
  caller has made explicit are not features of their own;
  :func:`chemsuite.metabolism.sites.predict_sites` strips them and maps back.
* **Non-finite values become 0.0.** Gasteiger charges diverge on a few valences and come
  back NaN or infinite; upstream's ``_safe`` guard is kept, because a NaN here would
  reach the trees as a missing value and take the other branch.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.EState import EStateIndices

from . import ISOFORMS

#: Elements with a column of their own; anything else lands in ``el_other``.
ELEMENTS = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I")
#: Hybridizations with a column of their own; S, SP3D and unspecified are ``hyb_other``.
HYBRIDS = (
    Chem.HybridizationType.SP,
    Chem.HybridizationType.SP2,
    Chem.HybridizationType.SP3,
)

#: The per-atom columns, before the isoform block.
ATOM_FEATURE_NAMES = (
    *(f"el_{e}" for e in ELEMENTS),
    "el_other",
    "hyb_sp",
    "hyb_sp2",
    "hyb_sp3",
    "hyb_other",
    "n_hydrogens",
    "heavy_degree",
    "formal_charge",
    "is_aromatic",
    "in_ring",
    "min_ring_size",
    "n_rings_at_atom",
    "aromatic_neighbours",
    "hetero_neighbours",
    "has_double_bond",
    "has_triple_bond",
    "is_benzylic_or_allylic_CH",
    "is_terminal",
    "gasteiger_charge",
    "estate_index",
    "crippen_logp",
    "crippen_mr",
    "neighbour_charge_max",
    "neighbour_charge_min",
    "neighbour_charge_mean",
)

#: Every column the model reads, in order: the atom block, then one isoform one-hot.
FEATURE_NAMES = (*ATOM_FEATURE_NAMES, *(f"isoform_{i}" for i in ISOFORMS))

#: Columns in the atom block (34) and in the model's input (43).
N_ATOM_FEATURES = len(ATOM_FEATURE_NAMES)
N_FEATURES = len(FEATURE_NAMES)


def _safe(value: float) -> float:
    """*value* as a float, or 0.0 when it is NaN or infinite (upstream's guard)."""
    number = float(value)
    return number if np.isfinite(number) else 0.0


def _one_hot(value: object, options: tuple) -> list[float]:
    """A one-hot over *options* with a trailing "other" column."""
    row = [0.0] * (len(options) + 1)
    row[options.index(value) if value in options else -1] = 1.0
    return row


def _gasteiger(mol: Chem.Mol) -> list[float]:
    """Gasteiger partial charges, computed on a copy so *mol* keeps no properties."""
    work = Chem.Mol(mol)
    AllChem.ComputeGasteigerCharges(work)
    return [_safe(atom.GetDoubleProp("_GasteigerCharge")) for atom in work.GetAtoms()]


def _ring_context(mol: Chem.Mol) -> tuple[list[int], list[int]]:
    """Per atom: its smallest ring's size (0 if in none) and how many rings it is in."""
    rings = mol.GetRingInfo().AtomRings()
    smallest, count = [], []
    for atom in mol.GetAtoms():
        sizes = [len(ring) for ring in rings if atom.GetIdx() in ring]
        smallest.append(min(sizes) if sizes else 0)
        count.append(len(sizes))
    return smallest, count


def atom_block(mol: Chem.Mol) -> np.ndarray:
    """The 34 per-atom columns of *mol*, one row per heavy atom in its own order.

    *mol* must be sanitized and have its hydrogens implicit, as
    ``Chem.MolFromSmiles`` returns it.
    """
    smallest_ring, ring_count = _ring_context(mol)
    charges = _gasteiger(mol)
    estate = [_safe(value) for value in EStateIndices(mol)]
    # The per-atom Crippen contributions have no public name in RDKit; the private one
    # is what every user of them calls, upstream included.
    crippen = rdMolDescriptors._CalcCrippenContribs(mol)

    rows = []
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        neighbors = list(atom.GetNeighbors())
        aromatic_neighbors = sum(1 for n in neighbors if n.GetIsAromatic())
        hetero_neighbors = sum(1 for n in neighbors if n.GetSymbol() not in ("C", "H"))
        bonds = list(atom.GetBonds())
        has_double = any(b.GetBondType() == Chem.BondType.DOUBLE for b in bonds)
        has_triple = any(b.GetBondType() == Chem.BondType.TRIPLE for b in bonds)
        # A C-H next to an aromatic ring: the easy hydrogen abstraction, and the one
        # feature in the block that names a mechanism rather than describing the atom.
        benzylic = (
            atom.GetSymbol() == "C"
            and not atom.GetIsAromatic()
            and atom.GetTotalNumHs() > 0
            and aromatic_neighbors > 0
        )
        neighbor_charges = [charges[n.GetIdx()] for n in neighbors] or [0.0]
        logp, molar_refractivity = crippen[index]
        rows.append(
            _one_hot(atom.GetSymbol(), ELEMENTS)
            + _one_hot(atom.GetHybridization(), HYBRIDS)
            + [
                float(atom.GetTotalNumHs()),
                float(atom.GetDegree()),
                float(atom.GetFormalCharge()),
                float(atom.GetIsAromatic()),
                float(atom.IsInRing()),
                float(smallest_ring[index]),
                float(ring_count[index]),
                float(aromatic_neighbors),
                float(hetero_neighbors),
                float(has_double),
                float(has_triple),
                float(benzylic),
                float(atom.GetDegree() == 1),
                charges[index],
                estate[index],
                _safe(logp),
                _safe(molar_refractivity),
                max(neighbor_charges),
                min(neighbor_charges),
                float(np.mean(neighbor_charges)),
            ]
        )
    return np.asarray(rows, dtype=np.float64).reshape(mol.GetNumAtoms(), N_ATOM_FEATURES)


#: Radii of the atom environments the applicability domain is measured over.
DOMAIN_RADII = (0, 1, 2)


def atom_environments(mol: Chem.Mol) -> list[tuple[tuple[int, int], ...]]:
    """Per atom, its Morgan environment identifiers as ``(radius, identifier)`` pairs.

    The applicability domain is built and queried from these: an atom whose environments
    all appear in the training set is one the model has seen the like of. RDKit reports
    an environment only where it adds atoms to the smaller one, so an atom near the end
    of a chain has fewer than :data:`DOMAIN_RADII` entries — hence pairs rather than one
    identifier per radius.
    """
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=max(DOMAIN_RADII))
    output = rdFingerprintGenerator.AdditionalOutput()
    output.AllocateBitInfoMap()
    generator.GetSparseCountFingerprint(mol, additionalOutput=output)

    found: list[list[tuple[int, int]]] = [[] for _ in range(mol.GetNumAtoms())]
    for identifier, entries in output.GetBitInfoMap().items():
        for atom, radius in entries:
            found[atom].append((int(radius), int(identifier)))
    return [tuple(sorted(pairs)) for pairs in found]


def isoform_one_hot(isoform: str) -> np.ndarray:
    """The nine isoform columns for *isoform*, which must be one of :data:`ISOFORMS`."""
    if isoform not in ISOFORMS:
        raise ValueError(f"unknown isoform: {isoform!r}")
        # (A silent all-zero row would predict from a state the model never saw.)
    return np.asarray([float(isoform == i) for i in ISOFORMS])


def design_matrix(block: np.ndarray, isoforms: tuple[str, ...]) -> np.ndarray:
    """*block* repeated once per isoform, each copy with its one-hot appended.

    Rows run isoform-major: every atom of ``isoforms[0]``, then of ``isoforms[1]``, …
    which is the order :func:`chemsuite.metabolism.sites.predict_sites` reads back.
    """
    atoms = block.shape[0]
    stacked = np.empty((atoms * len(isoforms), N_FEATURES), dtype=np.float64)
    for position, isoform in enumerate(isoforms):
        rows = slice(position * atoms, (position + 1) * atoms)
        stacked[rows, :N_ATOM_FEATURES] = block
        stacked[rows, N_ATOM_FEATURES:] = isoform_one_hot(isoform)
    return stacked
