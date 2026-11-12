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
"""Chain sequences against the full (SEQRES) sequence, so unmodelled loops show as gaps.

Each polymer chain gives a :class:`ChainSequence`: one letter per position of its entity's
full sequence and, for each position, the residue modelled there (or -1). A chain whose file
had no SEQRES gets the residues it has, in order.

Everything outside the polymers comes as :class:`KindRow` s, one per kind (ligands, glycans,
additives, ions, waters), each residue in its own position, so a single water can be found,
selected and deleted like a residue of a chain.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

import gemmi
import numpy as np

from .model import ADDITIVE, AMINO, GLYCAN, ION, LIGAND, NUCLEIC, WATER, Structure


@cache
def one_letter(name: str) -> str:
    """``HIS`` -> ``H``; a modified residue gives its parent's letter (``SEP`` -> ``S``),
    anything else ``X``."""
    info = gemmi.find_tabulated_residue(name)
    if info is None:
        return "X"
    code = info.one_letter_code
    return code.upper() if code.strip() else "X"


@dataclass(frozen=True)
class ChainSequence:
    chain: int  # chain index in the structure
    names: np.ndarray  # <U5 residue name at each position
    residue_at: np.ndarray  # int32: residue index modelled at each position, -1 if missing

    @property
    def letters(self) -> str:
        return "".join(one_letter(str(n)) for n in self.names)

    @property
    def modelled(self) -> np.ndarray:
        return self.residue_at >= 0


def chain_sequences(structure: Structure) -> list[ChainSequence]:
    out = []
    polymer = np.isin(structure.res_kind, (AMINO, NUCLEIC))
    for c in range(structure.n_chains):
        span = structure.residues_of(c)
        residues = np.arange(span.start, span.stop)[polymer[span]]
        if not len(residues):
            continue
        out.append(_chain(structure, c, residues))
    return out


def _chain(structure: Structure, chain: int, residues: np.ndarray) -> ChainSequence:
    entities = structure.entity[residues]
    entity = int(np.bincount(entities[entities >= 0]).argmax()) if (entities >= 0).any() else -1
    full = structure.sequence_of(entity) if entity >= 0 else np.zeros(0, "<U5")
    seqres = structure.seqres[residues]
    placed = (seqres >= 0) & (seqres < len(full)) & (entities == entity)
    if len(full) and placed.any():
        residue_at = np.full(len(full), -1, dtype=np.int32)
        # the first residue placed at a position wins (microheterogeneity keeps one)
        positions, first = np.unique(seqres[placed], return_index=True)
        residue_at[positions] = residues[placed][first]
        return ChainSequence(chain, full.copy(), residue_at)
    return ChainSequence(chain, structure.res_name[residues].copy(), residues.astype(np.int32))


#: the kind rows, in the order they are shown, with their plural titles
KIND_ROWS = {LIGAND: "Ligands", GLYCAN: "Glycans", ADDITIVE: "Additives", ION: "Ions",
             WATER: "Waters"}  # fmt: skip


@dataclass(frozen=True)
class KindRow:
    """The residues of one kind outside the polymers, by chain and then number."""

    kind: int
    names: np.ndarray  # <U5 residue name at each position
    residue_at: np.ndarray  # int32: the residue at each position (every one is modelled)

    @property
    def title(self) -> str:
        return KIND_ROWS[self.kind]

    @property
    def modelled(self) -> np.ndarray:
        return np.ones(len(self.residue_at), bool)


def kind_rows(structure: Structure) -> list[KindRow]:
    """A row for each kind of non-polymer residue the structure has, in :data:`KIND_ROWS`
    order."""
    out = []
    for kind in KIND_ROWS:
        residues = np.flatnonzero(structure.res_kind == kind)
        if not len(residues):
            continue
        order = np.lexsort((structure.icode[residues], structure.res_seq[residues],
                            structure.chain_of[residues]))  # fmt: skip
        residues = residues[order].astype(np.int32)
        out.append(KindRow(kind, structure.res_name[residues].copy(), residues))
    return out
