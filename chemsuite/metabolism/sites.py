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

"""Which atoms of a molecule a CYP isoform attacks, ranked.

The model is scored for one molecule and one or more isoforms at a time: the atom block
is computed once and repeated per isoform, so nine isoforms cost barely more than one.

Three things a caller should know:

* **Scores are per atom and per isoform**, each the model's probability that this atom is
  a *primary* site of metabolism for that isoform. They are not calibrated against each
  other across molecules, and they do not sum to one — ranking within a molecule is what
  was measured and what the tool shows.
* **Symmetry-equivalent atoms score identically**, so ties are broken by the lower atom
  index. Without that the displayed top site would wobble between runs of the same
  molecule.
* **Hydrogens are not sites.** The model is heavy-atom only: explicit hydrogens are
  removed before scoring, and the atom indices come back in the caller's own numbering.

Accuracy, measured molecule-grouped on ``PRIMARY_SOM`` labels, is in docs/SOM_PLAN.md:
for a molecule the model has not seen, the true site is the top-ranked atom about 59% of
the time and in the top three about 81%.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from importlib import resources

import numpy as np
from rdkit import Chem

from . import ISOFORMS
from .features import atom_block, atom_environments, design_matrix
from .trees import REQUIRED, Ensemble, probabilities

#: The vendored model, about 250 KB of tree arrays with the domain index beside them.
ARCHIVE = "som_gbm.npz"
#: Property holding an atom's index in the caller's molecule while hydrogens are removed.
_INDEX = "_somAtomIndex"


@dataclass(frozen=True)
class MetaboliteSite:
    """One atom's score for one isoform."""

    atom: int  # index in the molecule as given
    isoform: str
    score: float  # the model's probability, 0 to 1
    domain: float  # share of this atom's environments seen in training, 0 to 1

    @property
    def inside_domain(self) -> bool:
        """Every one of the atom's environments was seen while training."""
        return self.domain >= 1.0


@dataclass(frozen=True)
class Model:
    """The ensemble and the applicability-domain index it was trained with."""

    ensemble: Ensemble
    domain: dict[int, np.ndarray]


@cache
def model() -> Model:
    """The vendored model, loaded once per process."""
    path = resources.files(__package__) / "data" / ARCHIVE
    with resources.as_file(path) as file, np.load(file) as archive:
        return Model(
            ensemble=Ensemble.from_arrays({name: archive[name] for name in REQUIRED}),
            domain={
                int(name.removeprefix("domain_r")): archive[name]
                for name in archive.files
                if name.startswith("domain_r")
            },
        )


def _heavy(mol: Chem.Mol) -> tuple[Chem.Mol, list[int]]:
    """*mol* with every explicit hydrogen removed, and each kept atom's index in *mol*.

    ``RemoveAllHs``, not ``RemoveHs``: the latter keeps a hydrogen that carries
    information of its own (an isotope, a charge, a stereo marker), and a hydrogen left
    in the graph would change its neighbor's degree and so its features. Training
    molecules came from SMILES with implicit hydrogens throughout.
    """
    work = Chem.Mol(mol)
    for atom in work.GetAtoms():
        atom.SetIntProp(_INDEX, atom.GetIdx())
    work = Chem.RemoveAllHs(work)
    return work, [atom.GetIntProp(_INDEX) for atom in work.GetAtoms()]


def _domain(mol: Chem.Mol, known: dict[int, np.ndarray]) -> np.ndarray:
    """Per atom, the share of its Morgan environments the training set contains."""
    share = np.ones(mol.GetNumAtoms())
    for atom, pairs in enumerate(atom_environments(mol)):
        if not pairs:
            continue
        seen = [
            radius in known and bool(np.isin(identifier, known[radius]))
            for radius, identifier in pairs
        ]
        share[atom] = sum(seen) / len(seen)
    return share


def predict_sites(
    mol: Chem.Mol, isoforms: tuple[str, ...] = ISOFORMS
) -> dict[str, list[MetaboliteSite]]:
    """Every heavy atom of *mol* scored for each of *isoforms*, best site first.

    The lists are ranked by score, ties broken by the lower atom index.
    """
    if not isoforms:
        raise ValueError("no isoform to predict for")
    unknown = [name for name in isoforms if name not in ISOFORMS]
    if unknown:
        raise ValueError(f"unknown isoform(s): {', '.join(unknown)}")

    heavy, indices = _heavy(mol)
    if heavy.GetNumAtoms() == 0:
        return {isoform: [] for isoform in isoforms}

    loaded = model()
    block = atom_block(heavy)
    scores = probabilities(design_matrix(block, isoforms), loaded.ensemble)
    domain = _domain(heavy, loaded.domain)

    atoms = heavy.GetNumAtoms()
    ranked: dict[str, list[MetaboliteSite]] = {}
    for position, isoform in enumerate(isoforms):
        column = scores[position * atoms : (position + 1) * atoms]
        order = np.lexsort((np.arange(atoms), -column))
        ranked[isoform] = [
            MetaboliteSite(
                atom=indices[int(i)],
                isoform=isoform,
                score=float(column[i]),
                domain=float(domain[i]),
            )
            for i in order
        ]
    return ranked


def top_site(mol: Chem.Mol, isoform: str) -> MetaboliteSite | None:
    """The highest-scoring atom of *mol* for *isoform*, or None if it has no atoms."""
    ranked = predict_sites(mol, (isoform,))[isoform]
    return ranked[0] if ranked else None
