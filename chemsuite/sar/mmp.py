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
"""Matched molecular pairs: compounds that differ by one small change at one place.

Each molecule is cut at every single acyclic bond (RDKit MMPA), and each hydrogen is also
treated as a cut point (so H → Cl pairs are found). The larger piece is the constant part;
the smaller (at most *max_variable* heavy atoms) is the change. Compounds sharing a constant
part pair up; pairs found through several cuts keep the smallest change. Transforms
(A → B) are then summarised over every pair that shows them.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMMPA

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text

Cell = bytes | str | None
HYDROGEN = "[*:1][H]"
MAX_GROUP = 2000  # compounds sharing one constant part; larger groups are trivial cores
CHUNK_ROWS = 200


@dataclass(frozen=True)
class Pair:
    a: int  # input row index (the "from" compound)
    b: int
    constant: str
    from_group: str
    to_group: str
    delta: float  # activity(b) − activity(a)


@dataclass
class Transform:
    from_group: str
    to_group: str
    pairs: list[Pair]

    @property
    def count(self) -> int:
        return len(self.pairs)

    def deltas(self) -> np.ndarray:
        return np.array([p.delta for p in self.pairs], dtype=np.float64)

    @property
    def mean(self) -> float:
        d = self.deltas()
        d = d[np.isfinite(d)]
        return float(d.mean()) if len(d) else math.nan

    @property
    def sd(self) -> float:
        d = self.deltas()
        d = d[np.isfinite(d)]
        return float(d.std(ddof=1)) if len(d) > 1 else math.nan

    @property
    def ups(self) -> int:
        return int(np.sum(self.deltas() > 0))

    @property
    def downs(self) -> int:
        return int(np.sum(self.deltas() < 0))

    def rows(self) -> set[int]:
        return {r for p in self.pairs for r in (p.a, p.b)}


def _heavy(smiles: str) -> int:
    return sum(1 for ch in Chem.MolFromSmiles(smiles).GetAtoms() if ch.GetAtomicNum() > 1)


def _canonical(fragment: str) -> str:
    mol = Chem.MolFromSmiles(fragment)
    return Chem.MolToSmiles(mol) if mol is not None else fragment


def fragment(mol: Chem.Mol, max_variable: int = 10) -> list[tuple[str, str]]:
    """``(constant, change)`` SMILES pairs for one molecule, attachment written ``[*:1]``."""
    out: set[tuple[str, str]] = set()
    total = mol.GetNumHeavyAtoms()
    with rdBase.BlockLogs():
        for _core, chains in rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False):
            pieces = chains.split(".")
            if len(pieces) != 2:
                continue
            sizes = [_heavy(p) for p in pieces]
            for k in (0, 1):
                constant, change = pieces[1 - k], pieces[k]
                if sizes[k] <= max_variable and sizes[k] <= sizes[1 - k]:
                    out.add((_canonical(constant), _canonical(change)))
        # Hydrogen as a cut point: each symmetry-distinct H becomes [*:1].
        seen: set[str] = set()
        for atom in mol.GetAtoms():
            if atom.GetTotalNumHs() == 0 or atom.GetAtomicNum() <= 1:
                continue
            edit = Chem.RWMol(mol)
            dummy = edit.AddAtom(Chem.Atom(0))
            edit.GetAtomWithIdx(dummy).SetAtomMapNum(1)
            edit.AddBond(atom.GetIdx(), dummy, Chem.BondType.SINGLE)
            target = edit.GetAtomWithIdx(atom.GetIdx())
            target.SetNumExplicitHs(atom.GetTotalNumHs() - 1)  # from the unedited molecule
            target.SetNoImplicit(True)
            try:
                Chem.SanitizeMol(edit)
            except (ValueError, RuntimeError):
                continue
            constant = Chem.MolToSmiles(edit)
            if constant not in seen and total >= 1:
                seen.add(constant)
                out.add((constant, HYDROGEN))
    return sorted(out)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def fragment_chunk(args: tuple[int, int, Sequence[Cell]]) -> list[tuple[int, str, str]]:
    """``(row, constant, change)`` records for a chunk (a picklable worker entry point)."""
    start, max_variable, cells = args
    out = []
    for offset, cell in enumerate(cells):
        mol = _mol(cell)
        if mol is None:
            continue
        if mol.GetNumConformers() and mol.GetConformer().Is3D():
            mol = Chem.RemoveHs(mol)
        try:
            for constant, change in fragment(mol, max_variable):
                out.append((start + offset, constant, change))
        except (ValueError, RuntimeError):
            continue
    return out


def _chunks(cells: Sequence[Cell], max_variable: int, size: int) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        yield start, max_variable, list(cells[start : start + size])


def find_pairs(
    cells: Sequence[Cell],
    activity: np.ndarray,
    *,
    max_variable: int = 10,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> list[Pair]:
    """Every matched pair (each unordered pair once, the smallest change kept)."""
    activity = np.asarray(activity, dtype=np.float64)
    groups: dict[str, list[tuple[int, str]]] = defaultdict(list)
    done = 0
    for records in ordered_map(fragment_chunk, _chunks(cells, max_variable, CHUNK_ROWS), executor):
        check(cancel)
        for row, constant, change in records:
            groups[constant].append((row, change))
        done = min(len(cells), done + CHUNK_ROWS)
        if progress is not None:
            progress(done, len(cells))
    size_cache: dict[str, int] = {}

    def size(smiles: str) -> int:
        if smiles not in size_cache:
            size_cache[smiles] = 0 if smiles == HYDROGEN else _heavy(smiles)
        return size_cache[smiles]

    best: dict[tuple[int, int], tuple[int, Pair]] = {}
    for constant, members in groups.items():
        check(cancel)
        if len(members) < 2 or len(members) > MAX_GROUP:
            continue
        for x in range(len(members)):
            row_a, change_a = members[x]
            for y in range(x + 1, len(members)):
                row_b, change_b = members[y]
                if row_a == row_b or change_a == change_b:
                    continue
                a, b, from_g, to_g = (
                    (row_a, row_b, change_a, change_b)
                    if row_a < row_b
                    else (row_b, row_a, change_b, change_a)
                )
                weight = size(from_g) + size(to_g)
                current = best.get((a, b))
                if current is None or weight < current[0]:
                    pair = Pair(a, b, constant, from_g, to_g, float(activity[b] - activity[a]))
                    best[(a, b)] = (weight, pair)
    return [pair for _w, pair in sorted(best.values(), key=lambda wp: (wp[1].a, wp[1].b))]


def transforms(pairs: Sequence[Pair], *, min_pairs: int = 1) -> list[Transform]:
    """Group pairs by change. A → B and B → A are one transform, oriented from hydrogen
    (H → X), otherwise from the smaller group."""
    grouped: dict[tuple[str, str], list[Pair]] = defaultdict(list)
    for p in pairs:
        if _orient_key(p.from_group) <= _orient_key(p.to_group):
            grouped[(p.from_group, p.to_group)].append(p)
        else:
            flipped = Pair(p.b, p.a, p.constant, p.to_group, p.from_group, -p.delta)
            grouped[(p.to_group, p.from_group)].append(flipped)
    out = [Transform(a, b, ps) for (a, b), ps in grouped.items() if len(ps) >= min_pairs]
    out.sort(key=lambda t: (-t.count, t.from_group, t.to_group))
    return out


def _orient_key(group: str) -> tuple:
    return (group != HYDROGEN, len(group), group)


def readable(group: str) -> str:
    """``FC(F)(F)[*:1]`` → ``*C(F)(F)F``; hydrogen → ``H`` (the R-group labels' style)."""
    from .rgroups import substituent_label

    return substituent_label(group)
