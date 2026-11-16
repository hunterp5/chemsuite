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
"""A protein structure as immutable numpy arrays (struct of arrays).

Atoms, residues and chains are each a set of parallel arrays. Residues are runs of atoms and
chains runs of residues, given by ``res_start`` and ``chain_start`` (one longer than the
number of residues or chains), so the atoms of residue ``i`` are
``res_start[i]:res_start[i + 1]``.

Every atom has a ``uid`` that survives edits: selections, measurements and restraints refer to
atoms by uid, never by position. All alternate locations are kept; :meth:`Structure.active`
says which atoms are shown. NMR models (and later trajectory frames) share one set of atom
arrays and differ only in ``coords[k]``.

A structure is never changed in place: an edit builds a new one that shares the arrays it did
not touch (:func:`dataclasses.replace`, :meth:`Structure.subset`).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, fields, replace
from functools import cached_property

import numpy as np

# Residue kinds (``res_kind``), decided by :mod:`.classify`.
AMINO = 0  # amino acid in a polymer, standard or modified
NUCLEIC = 1  # nucleotide in a polymer
LIGAND = 2
WATER = 3
ION = 4  # a single metal or halide atom
ADDITIVE = 5  # buffer, cryoprotectant, crystallisation additive
GLYCAN = 6  # sugar attached to the protein
KIND_NAMES = {
    AMINO: "Protein",
    NUCLEIC: "Nucleic Acid",
    LIGAND: "Ligand",
    WATER: "Water",
    ION: "Ion",
    ADDITIVE: "Additive",
    GLYCAN: "Glycan",
}

# Secondary structure (``ss``).
COIL = 0
HELIX = 1
STRAND = 2

# Bond orders: 1, 2, 3, and this for aromatic.
AROMATIC = 4

#: atom names of the protein backbone (with its hydrogens and the C-terminal oxygen)
BACKBONE = ("N", "CA", "C", "O", "OXT", "H", "H2", "H3", "HA", "HA2", "HA3")
#: and of the nucleic-acid backbone
NUCLEIC_BACKBONE = ("P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'",
                    "O2'", "C1'")  # fmt: skip


def _frozen(array: np.ndarray) -> np.ndarray:
    array.flags.writeable = False
    return array


@dataclass(frozen=True, eq=False)
class Structure:
    # atoms (n)
    uid: np.ndarray  # int64
    coords: np.ndarray  # float32 (k, n, 3): k models
    element: np.ndarray  # uint8 atomic number (0 unknown)
    name: np.ndarray  # <U4
    altloc: np.ndarray  # <U1, "" when the atom has one location
    occupancy: np.ndarray  # float32
    b: np.ndarray  # float32
    charge: np.ndarray  # int8
    # residues (r)
    res_start: np.ndarray  # int32 (r + 1)
    res_name: np.ndarray  # <U5
    res_seq: np.ndarray  # int32 author numbering
    icode: np.ndarray  # <U1 insertion code
    res_kind: np.ndarray  # uint8, the kinds above
    ss: np.ndarray  # uint8 secondary structure
    seqres: np.ndarray  # int32: position in the entity's full sequence, -1 if not in one
    label_asym: np.ndarray  # <U4 mmCIF label_asym_id (subchain)
    entity: np.ndarray  # int16: index into the entities, -1 if none
    # chains (c), by author chain name
    chain_start: np.ndarray  # int32 (c + 1)
    chain_name: np.ndarray  # <U4
    # entities (e): the full (SEQRES) sequence of each, as residue names
    entity_id: np.ndarray  # <U8
    entity_seq: np.ndarray  # <U5, all sequences one after the other
    entity_seq_start: np.ndarray  # int32 (e + 1)
    # bonds (m)
    bonds: np.ndarray  # int32 (m, 2), atom indices, i < j
    order: np.ndarray  # int8: 1, 2, 3 or AROMATIC
    next_uid: int = 0
    meta: Mapping[str, str] = field(default_factory=dict)  # id, title, method, resolution

    def __post_init__(self) -> None:
        for value in self._arrays():
            _frozen(value)

    def _arrays(self) -> list[np.ndarray]:
        values = (self.__dict__[f.name] for f in fields(self))
        return [v for v in values if isinstance(v, np.ndarray)]

    # sizes
    @property
    def n_atoms(self) -> int:
        return len(self.uid)

    @property
    def n_residues(self) -> int:
        return len(self.res_start) - 1

    @property
    def n_chains(self) -> int:
        return len(self.chain_start) - 1

    @property
    def n_models(self) -> int:
        return len(self.coords)

    @property
    def nbytes(self) -> int:
        return sum(v.nbytes for v in self._arrays())

    # derived per-atom and per-residue indices (computed once)
    @cached_property
    def residue_of(self) -> np.ndarray:
        """The residue index of every atom."""
        sizes = np.diff(self.res_start)
        return _frozen(np.repeat(np.arange(self.n_residues, dtype=np.int32), sizes))

    @cached_property
    def chain_of(self) -> np.ndarray:
        """The chain index of every residue."""
        sizes = np.diff(self.chain_start)
        return _frozen(np.repeat(np.arange(self.n_chains, dtype=np.int32), sizes))

    @cached_property
    def atom_chain(self) -> np.ndarray:
        return _frozen(self.chain_of[self.residue_of])

    @cached_property
    def atom_kind(self) -> np.ndarray:
        return _frozen(self.res_kind[self.residue_of])

    @cached_property
    def index_of_uid(self) -> dict[int, int]:
        return {int(u): i for i, u in enumerate(self.uid)}

    def atoms_of(self, residue: int) -> slice:
        return slice(int(self.res_start[residue]), int(self.res_start[residue + 1]))

    def residues_of(self, chain: int) -> slice:
        return slice(int(self.chain_start[chain]), int(self.chain_start[chain + 1]))

    def indices(self, uids: np.ndarray) -> np.ndarray:
        """Atom indices of these uids (-1 for a uid no longer in the structure)."""
        lookup = self.index_of_uid
        return np.array([lookup.get(int(u), -1) for u in uids], dtype=np.int64)

    def sequence_of(self, entity: int) -> np.ndarray:
        start, end = self.entity_seq_start[entity], self.entity_seq_start[entity + 1]
        return self.entity_seq[start:end]

    def residue_label(self, residue: int) -> str:
        """``A/HIS 57`` (with the insertion code, if any)."""
        chain = self.chain_name[self.chain_of[residue]]
        return f"{chain}/{self.res_name[residue]} {self.res_seq[residue]}{self.icode[residue]}"

    def atom_label(self, atom: int) -> str:
        alt = f" ({self.altloc[atom]})" if self.altloc[atom] else ""
        return f"{self.residue_label(int(self.residue_of[atom]))}/{self.name[atom]}{alt}"

    # alternate locations
    @cached_property
    def default_altlocs(self) -> dict[int, str]:
        """For each residue with alternate locations, the one shown unless the user picks
        another: the highest mean occupancy, then the first letter."""
        has = np.flatnonzero(self.altloc != "")
        residues = self.residue_of[has]
        choice: dict[int, str] = {}
        for residue in np.unique(residues):
            atoms = has[residues == residue]
            letters = self.altloc[atoms]
            choice[int(residue)] = max(
                np.unique(letters),
                key=lambda x: (float(self.occupancy[atoms][letters == x].mean()), -ord(x)),
            )
        return choice

    def active(self, altlocs: Mapping[int, str] | None = None) -> np.ndarray:
        """Atoms shown: those with one location, and the chosen location of the others.

        ``altlocs`` maps a residue index to the letter the user picked for it."""
        mask = self.altloc == ""
        if mask.all():
            return mask
        chosen = dict(self.default_altlocs)
        if altlocs:
            chosen.update(altlocs)
        has = np.flatnonzero(~mask)
        letters = np.array([chosen.get(int(r), "") for r in self.residue_of[has]], dtype="<U1")
        mask[has] = self.altloc[has] == letters
        return mask

    def xyz(self, model: int = 0) -> np.ndarray:
        return self.coords[model]

    # building new structures
    def with_coords(self, coords: np.ndarray) -> Structure:
        coords = np.asarray(coords, dtype=np.float32)
        if coords.ndim == 2:
            coords = coords[None]
        return replace(self, coords=coords)

    def subset(self, keep: np.ndarray) -> Structure:
        """A structure with only the atoms where ``keep`` is True (uids kept); residues and
        chains left without atoms go too."""
        keep = np.asarray(keep, dtype=bool)
        if keep.all():
            return self
        atoms = np.flatnonzero(keep)
        new_index = np.full(self.n_atoms, -1, dtype=np.int64)
        new_index[atoms] = np.arange(len(atoms))
        # residues that keep at least one atom
        res_count = np.bincount(self.residue_of[atoms], minlength=self.n_residues)
        residues = np.flatnonzero(res_count)
        res_start = np.concatenate([[0], np.cumsum(res_count[residues])]).astype(np.int32)
        chain_count = np.bincount(self.chain_of[residues], minlength=self.n_chains)
        chains = np.flatnonzero(chain_count)
        chain_start = np.concatenate([[0], np.cumsum(chain_count[chains])]).astype(np.int32)
        both = keep[self.bonds].all(axis=1) if len(self.bonds) else np.zeros(0, dtype=bool)
        bonds = new_index[self.bonds[both]].astype(np.int32).reshape(-1, 2)
        return replace(
            self,
            uid=self.uid[atoms],
            coords=self.coords[:, atoms],
            element=self.element[atoms],
            name=self.name[atoms],
            altloc=self.altloc[atoms],
            occupancy=self.occupancy[atoms],
            b=self.b[atoms],
            charge=self.charge[atoms],
            res_start=res_start,
            res_name=self.res_name[residues],
            res_seq=self.res_seq[residues],
            icode=self.icode[residues],
            res_kind=self.res_kind[residues],
            ss=self.ss[residues],
            seqres=self.seqres[residues],
            label_asym=self.label_asym[residues],
            entity=self.entity[residues],
            chain_start=chain_start,
            chain_name=self.chain_name[chains],
            bonds=bonds,
            order=self.order[both],
        )

    def to_bytes(self) -> bytes:
        from .store import pack  # the store reads this module; imported when first needed

        return pack(self)


def empty() -> Structure:
    """A structure without atoms."""
    return Structure(
        uid=np.zeros(0, np.int64),
        coords=np.zeros((1, 0, 3), np.float32),
        element=np.zeros(0, np.uint8),
        name=np.zeros(0, "<U4"),
        altloc=np.zeros(0, "<U1"),
        occupancy=np.zeros(0, np.float32),
        b=np.zeros(0, np.float32),
        charge=np.zeros(0, np.int8),
        res_start=np.zeros(1, np.int32),
        res_name=np.zeros(0, "<U5"),
        res_seq=np.zeros(0, np.int32),
        icode=np.zeros(0, "<U1"),
        res_kind=np.zeros(0, np.uint8),
        ss=np.zeros(0, np.uint8),
        seqres=np.zeros(0, np.int32),
        label_asym=np.zeros(0, "<U4"),
        entity=np.zeros(0, np.int16),
        chain_start=np.zeros(1, np.int32),
        chain_name=np.zeros(0, "<U4"),
        entity_id=np.zeros(0, "<U8"),
        entity_seq=np.zeros(0, "<U5"),
        entity_seq_start=np.zeros(1, np.int32),
        bonds=np.zeros((0, 2), np.int32),
        order=np.zeros(0, np.int8),
    )
