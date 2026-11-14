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

"""R-groups around each molecule's own most central ring system.

:func:`decompose_central_ring` needs no core: it finds one. A ring's **centrality** is the
mean bond-path distance from its atoms to every heavy atom, and the lowest mean wins; ties
prefer the larger ring, then the ring holding the lowest canonical rank. Rings that share an
atom with the winner — fused, bridged or spiro — join it, so an indole stays whole; a ring
joined only by a bond, such as biphenyl's second ring, does not.

The bonds leaving that system are cut and both ends labelled ``[*:1]``, ``[*:2]``… in the
core's canonical rank order, so the same scaffold numbers its positions the same way. The
substituent labels are those of :mod:`chemsuite.sar.rgroups` (``*Cl``, ``*OC``), so the
columns read the same as R-Group Decomposition's.

Each molecule is decomposed on its own, so two rows share an R number only when they share a
core: ``Core`` says which one each row got.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

from rdkit import Chem, rdBase

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text
from chemsuite.sar.rgroups import substituent_label

#: A structure cell as sent to a worker.
Cell = bytes | str | None
Progress = Callable[[int, int], None]

CHUNK_ROWS = 200
INLINE_CHUNK_ROWS = 20
POOL_MIN_ROWS = 500
#: a ring system smaller than this is not worth calling a core
MIN_CORE_ATOMS = 3
_FAILURES = (ValueError, RuntimeError)


@dataclass
class CentralRings:
    """Per row: the core and its substituent labels. Blank where there was no ring."""

    cores: list[str] = field(default_factory=list)
    #: position name ("R1") → label per row ("" where the row has nothing there)
    labels: dict[str, list[str]] = field(default_factory=dict)
    failed: int = 0

    def __len__(self) -> int:
        return len(self.cores)

    def positions(self) -> tuple[str, ...]:
        """The position names in order (R1, R2, …)."""
        return tuple(sorted(self.labels, key=lambda name: int(name[1:])))

    def extend(self, other: CentralRings) -> None:
        here = len(self.cores)
        self.cores.extend(other.cores)
        for name, values in other.labels.items():
            column = self.labels.setdefault(name, [""] * here)
            column.extend([""] * (here - len(column)))
            column.extend(values)
        for name, column in self.labels.items():
            if name not in other.labels:
                column.extend([""] * len(other.cores))
        self.failed += other.failed

    def add(self, core: str, groups: Sequence[str]) -> None:
        """One row's result."""
        here = len(self.cores)
        self.cores.append(core)
        for index, label in enumerate(groups, start=1):
            name = f"R{index}"
            column = self.labels.setdefault(name, [""] * here)
            column.extend([""] * (here - len(column)))
            column.append(label)
        for column in self.labels.values():
            column.extend([""] * (here + 1 - len(column)))

    def decomposed(self) -> int:
        return sum(1 for core in self.cores if core)


# -- finding the ring system -------------------------------------------------------------------


def _rings(mol: Chem.Mol) -> list[tuple[int, ...]]:
    try:
        found = Chem.GetSymmSSSR(mol)
    except _FAILURES:
        return []
    return [tuple(int(a) for a in ring) for ring in found if len(ring) >= 3]


def _heavy_atoms(mol: Chem.Mol) -> list[int]:
    return [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]


def _mean_distances(mol: Chem.Mol, heavy: Sequence[int]) -> list[float]:
    """Mean topological distance from each atom to every heavy atom."""
    matrix = Chem.GetDistanceMatrix(mol)
    means = [float("inf")] * mol.GetNumAtoms()
    for atom in range(mol.GetNumAtoms()):
        row = matrix[atom]
        taken = [float(row[target]) for target in heavy if float(row[target]) >= 0]
        if taken:
            means[atom] = sum(taken) / len(taken)
    return means


def _system(rings: Sequence[tuple[int, ...]], chosen: int) -> set[int]:
    """Every ring joined to *chosen* through a shared atom (fused, bridged or spiro)."""
    parent = list(range(len(rings)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    sets = [set(ring) for ring in rings]
    for left in range(len(rings)):
        for right in range(left + 1, len(rings)):
            if sets[left] & sets[right]:
                a, b = find(left), find(right)
                if a != b:
                    parent[b] = a
    root = find(chosen)
    atoms: set[int] = set()
    for index, ring in enumerate(sets):
        if find(index) == root:
            atoms |= ring
    return atoms


def most_central_ring_atoms(mol: Chem.Mol) -> tuple[int, ...] | None:
    """The heavy atoms of the ring system holding the most central ring, or None."""
    rings = _rings(mol)
    heavy = _heavy_atoms(mol)
    if not rings or not heavy:
        return None
    means = _mean_distances(mol, heavy)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True, includeChirality=True))
    chosen = min(
        range(len(rings)),
        key=lambda i: (
            sum(means[atom] for atom in rings[i]) / len(rings[i]),
            -len(rings[i]),
            min(ranks[atom] for atom in rings[i]),
        ),
    )
    atoms = {a for a in _system(rings, chosen) if mol.GetAtomWithIdx(a).GetAtomicNum() > 1}
    if len(atoms) < MIN_CORE_ATOMS:
        return None
    return tuple(sorted(atoms))


# -- cutting it out ------------------------------------------------------------------------------


def _prepared(mol: Chem.Mol) -> Chem.Mol | None:
    """A sanitized copy without explicit hydrogens; the caller's molecule is untouched."""
    try:
        work = Chem.Mol(mol)
        Chem.SanitizeMol(work)
    except _FAILURES:
        return None
    try:
        return Chem.RemoveHs(work) or work
    except _FAILURES:
        return work


def _cut_bonds(mol: Chem.Mol, core: set[int]) -> list[int]:
    """Bonds leaving the core for a heavy atom; hydrogens stay implicit."""
    cuts = []
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (begin in core) == (end in core):
            continue
        outside = end if begin in core else begin
        if mol.GetAtomWithIdx(outside).GetAtomicNum() > 1:
            cuts.append(bond.GetIdx())
    return cuts


def _links(mol: Chem.Mol) -> dict[int, tuple[int, int]]:
    """Each dummy atom → (its neighbour, its isotope).

    ``FragmentOnBonds`` labels a dummy with the atom index on the far side of the cut, so the
    two dummies of one bond point at each other.
    """
    links = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        neighbours = atom.GetNeighbors()
        if len(neighbours) == 1:
            links[atom.GetIdx()] = (neighbours[0].GetIdx(), int(atom.GetIsotope()))
    return links


def _partner(links: dict[int, tuple[int, int]], index: int) -> int | None:
    neighbour, far = links[index]
    for other, (its_neighbour, its_far) in links.items():
        if other != index and its_neighbour == far and its_far == neighbour:
            return other
    return None


def _subset(mol: Chem.Mol, atoms: Sequence[int]) -> tuple[Chem.Mol, dict[int, int]]:
    """The induced subgraph of *atoms*, with dummy isotopes cleared so rank is unaffected."""
    out = Chem.RWMol()
    old_to_new: dict[int, int] = {}
    keep = set(atoms)
    for old in atoms:
        atom = Chem.Atom(mol.GetAtomWithIdx(old))
        atom.SetAtomMapNum(0)
        if atom.GetAtomicNum() == 0:
            atom.SetIsotope(0)
        old_to_new[old] = out.AddAtom(atom)
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if begin in keep and end in keep:
            out.AddBond(old_to_new[begin], old_to_new[end], bond.GetBondType())
            made = out.GetBondBetweenAtoms(old_to_new[begin], old_to_new[end])
            if made is not None:
                made.SetIsAromatic(bond.GetIsAromatic())
    return out.GetMol(), old_to_new


def _rankable(mol: Chem.Mol) -> None:
    try:
        Chem.SanitizeMol(mol)
    except _FAILURES:
        mol.UpdatePropertyCache(strict=False)


def _label_cuts(mol: Chem.Mol, core: set[int], links: dict[int, tuple[int, int]]) -> bool:
    """Number both ends of each cut ``[*:1]``, ``[*:2]``… in the core's rank order."""
    core_dummies = [i for i, (neighbour, _far) in links.items() if neighbour in core]
    if not links or len(core_dummies) != len(links) // 2:
        return False
    ranked, old_to_new = _subset(mol, sorted(core) + sorted(core_dummies))
    _rankable(ranked)
    ranks = list(Chem.CanonicalRankAtoms(ranked, breakTies=True, includeChirality=True))
    used: set[int] = set()
    for number, index in enumerate(sorted(core_dummies, key=lambda i: ranks[old_to_new[i]]), 1):
        partner = _partner(links, index)
        if partner is None or index in used or partner in used:
            return False
        used |= {index, partner}
        mol.GetAtomWithIdx(index).SetAtomMapNum(number)
        mol.GetAtomWithIdx(partner).SetAtomMapNum(number)
    return len(used) == len(links)


def _split(mol: Chem.Mol, core: set[int]) -> tuple[Chem.Mol, list[Chem.Mol]] | None:
    """The core fragment and the substituent fragments that carry a dummy.

    A fragment holding both a core atom and another heavy atom means a substituent was not
    cut, and the decomposition is refused.
    """
    groups = Chem.GetMolFrags(mol)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    found: Chem.Mol | None = None
    substituents: list[Chem.Mol] = []
    for indices, frag in zip(groups, frags, strict=True):
        heavy = [i for i in indices if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
        if any(i in core for i in heavy):
            if any(i not in core for i in heavy):
                return None
            found = frag
        elif any(mol.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in indices):
            substituents.append(frag)
    return None if found is None else (found, substituents)


def _loose_dummy(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() == 0 and a.GetAtomMapNum() <= 0 for a in mol.GetAtoms())


def _order(frag: Chem.Mol) -> tuple[int, int]:
    maps = [
        a.GetAtomMapNum()
        for a in frag.GetAtoms()
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() > 0
    ]  # noqa: E501
    return (min(maps), len(maps)) if maps else (10**9, frag.GetNumAtoms())


def _smiles(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(copy)
    except _FAILURES:
        kekulize = Chem.SanitizeFlags.SANITIZE_KEKULIZE
        try:
            Chem.SanitizeMol(copy, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ kekulize)
        except _FAILURES:
            copy.UpdatePropertyCache(strict=False)
    try:
        return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True) or ""
    except _FAILURES:
        return ""


@dataclass(frozen=True)
class Split:
    """One molecule's core and the substituents on it, in position order."""

    core: str
    groups: tuple[str, ...]


def decompose_central_ring(mol: Chem.Mol | None) -> Split | None:
    """*mol* split into its most central ring system and the substituents on it, or None.

    The core keeps its numbered attachment points; each substituent is labelled as R-Group
    Decomposition labels it (``*Cl``, ``*OC``).
    """
    if mol is None or not mol.GetNumAtoms():
        return None
    with rdBase.BlockLogs():
        work = _prepared(mol)
        if work is None:
            return None
        atoms = most_central_ring_atoms(work)
        if not atoms:
            return None
        core = set(atoms)
        cuts = _cut_bonds(work, core)
        if cuts:
            kinds = [work.GetBondWithIdx(i).GetBondType() for i in cuts]
            cut = Chem.FragmentOnBonds(work, cuts, addDummies=True, bondTypes=kinds)
            links = _links(cut)
            for index in links:
                cut.GetAtomWithIdx(index).SetIsotope(0)
            if not _label_cuts(cut, core, links):
                return None
            split = _split(cut, core)
        else:
            split = _split(work, core)
        if split is None:
            return None
        core_mol, substituents = split
        if substituents and not any(a.GetAtomicNum() == 0 for a in core_mol.GetAtoms()):
            return None
        if any(_loose_dummy(frag) for frag in (core_mol, *substituents)):
            return None
        core_smiles = _smiles(core_mol)
        if not core_smiles:
            return None
        groups = [_smiles(frag) for frag in sorted(substituents, key=_order)]
        if any(not smiles for smiles in groups):
            return None
    return Split(core_smiles, tuple(substituent_label(s) for s in groups))


# -- over a column -------------------------------------------------------------------------------


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def central_ring_chunk(cells: Sequence[Cell]) -> CentralRings:
    """One chunk of cells (a picklable worker entry point)."""
    out = CentralRings()
    with rdBase.BlockLogs():
        for cell in cells:
            mol = _mol(cell)
            if mol is None:
                out.add("", ())
                out.failed += cell is not None
                continue
            split = decompose_central_ring(mol)
            out.add(split.core, split.groups) if split else out.add("", ())
    return out


def _chunked(total: int, executor: Executor | None) -> Iterator[slice]:
    size = CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS
    for start in range(0, total, size):
        yield slice(start, min(start + size, total))


def central_ring_rows(
    cells: Sequence[Cell],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> CentralRings:
    """Decompose every cell. Cancelled: ``partial`` holds the rows already finished."""
    cells = list(cells)
    total = len(cells)
    chunks = (cells[s] for s in _chunked(total, executor))
    out = CentralRings()
    try:
        for part in ordered_map(central_ring_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(len(out), total)
    except Cancelled:
        raise Cancelled(out if len(out) else None) from None
    return out


__all__ = [
    "MIN_CORE_ATOMS",
    "POOL_MIN_ROWS",
    "CentralRings",
    "Split",
    "central_ring_chunk",
    "central_ring_rows",
    "decompose_central_ring",
    "most_central_ring_atoms",
]
