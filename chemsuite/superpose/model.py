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
"""What a superposition takes and gives: the reference, the settings and the scores."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase

#: the methods, in the order the dialog lists them
MCS = "mcs"
CORE = "core"
O3A = "o3a"
SHAPE = "shape"
COLOR = "shape_color"
FIELD = "field"
METHODS = (MCS, CORE, O3A, SHAPE, COLOR, FIELD)

#: the score columns of each method; the first is the one the best conformer is chosen by
CORE_RMSD = "Core RMSD (Å)"
CORE_ATOMS = "Core Atoms"
O3A_SCORE = "O3A Score"
RMSD = "RMSD (Å)"
SHAPE_TANIMOTO = "Shape Tanimoto"
COLOR_TANIMOTO = "Color Tanimoto"
COMBO = "Combo"
FIELD_SIMILARITY = "Field Similarity"
SCORES: dict[str, tuple[str, ...]] = {
    MCS: (CORE_RMSD, CORE_ATOMS),
    CORE: (CORE_RMSD, CORE_ATOMS),
    O3A: (O3A_SCORE, RMSD),
    SHAPE: (SHAPE_TANIMOTO,),
    COLOR: (COMBO, SHAPE_TANIMOTO, COLOR_TANIMOTO),
    FIELD: (FIELD_SIMILARITY, SHAPE_TANIMOTO),
}
#: scores where lower is better (the rest: higher is better)
LOWER_IS_BETTER = frozenset({CORE_RMSD, RMSD})
#: decimals shown for each score
DECIMALS = {CORE_ATOMS: 0}

#: "each row's conformers onto its first": the scores it writes
MEAN_RMSD = "Mean Conformer RMSD (Å)"
MAX_RMSD = "Max Conformer RMSD (Å)"

#: the private molecule property naming the reference asset a superposed molecule was fitted to
REFERENCE_PROP = "_superposeReference"
#: the kind of that asset, whose payload is the reference's MOL block
REFERENCE_KIND = "superpose-reference"
#: blobs that keep that property (chemsuite's usual blobs keep none)
_PICKLE = Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps


@dataclass(frozen=True)
class Settings:
    """How to superpose.

    *core* is a SMARTS or SMILES (the ``core`` method); the ``mcs_*`` fields tune the
    maximum common substructure. *o3a_typing* is ``"mmff"`` or ``"crippen"``. *color_weight*
    is the share of the pharmacophore features in the shape + features optimization (0–1).
    *charges* is ``"gasteiger"`` or ``"mmff"``; *field_weight* the share of the field in the
    field method (0–1). *conformers* is how many to generate for a row that has only 2D.
    *heavy_only* fits a molecule's own conformers on heavy atoms only.
    """

    method: str = COLOR
    core: str = ""
    mcs_ring_matches_ring: bool = True
    mcs_complete_rings: bool = True
    mcs_compare_elements: bool = True
    mcs_timeout: int = 10
    o3a_typing: str = "mmff"
    color_weight: float = 0.5
    charges: str = "gasteiger"
    field_weight: float = 0.5
    conformers: int = 30
    heavy_only: bool = True

    def __post_init__(self) -> None:
        if self.method not in METHODS:
            raise ValueError(f"unknown method: {self.method!r}")
        if self.o3a_typing not in ("mmff", "crippen"):
            raise ValueError(f"unknown O3A typing: {self.o3a_typing!r}")
        if self.charges not in ("gasteiger", "mmff"):
            raise ValueError(f"unknown charges: {self.charges!r}")
        if not 0.0 <= self.color_weight <= 1.0 or not 0.0 <= self.field_weight <= 1.0:
            raise ValueError("weights must be between 0 and 1")
        if self.conformers < 1:
            raise ValueError("conformers must be at least 1")

    @property
    def scores(self) -> tuple[str, ...]:
        return SCORES[self.method]


@dataclass(frozen=True)
class Fit:
    """One conformer's best pose: the 4×4 transform that moves it onto the reference, and its
    scores (named as in :data:`SCORES`). *rank*, when given, orders fits instead of the
    method's first score (the field method ranks by what it optimized)."""

    conf_id: int
    transform: np.ndarray
    scores: dict[str, float]
    rank: float | None = None

    def key(self, method: str) -> float:
        """Larger is better."""
        if self.rank is not None:
            return self.rank
        name = SCORES[method][0]
        value = self.scores[name]
        return -value if name in LOWER_IS_BETTER else value


def best(fits: list[Fit], method: str) -> Fit | None:
    return max(fits, key=lambda f: f.key(method)) if fits else None


def has_3d(mol: Chem.Mol | None) -> bool:
    """Whether *mol* has a conformer with real 3D coordinates (a 2D molfile has a flat one)."""
    return mol is not None and mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()


def first_conformer(mol: Chem.Mol) -> Chem.Mol:
    """A copy of *mol* with only its first conformer, numbered 0."""
    out = Chem.Mol(mol)
    conf = Chem.Conformer(mol.GetConformer(mol.GetConformers()[0].GetId()))
    out.RemoveAllConformers()
    conf.SetId(0)
    out.AddConformer(conf, assignId=False)
    return out


def blob_with_reference(mol: Chem.Mol, reference_id: str) -> bytes:
    """*mol*'s blob, naming the reference asset it was fitted to (see :func:`reference_id_of`)."""
    out = Chem.Mol(mol)
    out.ClearComputedProps()
    out.SetProp(REFERENCE_PROP, reference_id)
    return out.ToBinary(_PICKLE)


def reference_id_of(mol: Chem.Mol | None) -> str | None:
    """The id of the reference asset *mol* was superposed onto, if it was."""
    if mol is None or not mol.HasProp(REFERENCE_PROP):
        return None
    return mol.GetProp(REFERENCE_PROP) or None


def reference_block(mol: Chem.Mol) -> bytes:
    """The reference asset's payload: *mol*'s first conformer as a MOL block."""
    return Chem.MolToMolBlock(first_conformer(mol)).encode()


def reference_from_block(payload: bytes) -> Chem.Mol | None:
    """The reference molecule read back from :func:`reference_block` (None if damaged)."""
    with rdBase.BlockLogs():
        return Chem.MolFromMolBlock(payload.decode(errors="replace"), removeHs=False)


def apply_transform(conf: Chem.Conformer, transform: np.ndarray) -> None:
    """Move *conf* by a 4×4 rigid transform."""
    xyz = conf.GetPositions()
    moved = xyz @ transform[:3, :3].T + transform[:3, 3]
    for i, (x, y, z) in enumerate(moved):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
