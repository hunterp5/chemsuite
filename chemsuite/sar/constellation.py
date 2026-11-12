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
"""Constellation plots: the cores analog series share, placed in chemical space.

The method follows the papers:

- Naveja, J. J., Vogt, M., Stumpfe, D., Medina-Franco, J. L. and Bajorath, J. (2019).
  "Systematic extraction of analogue series from large compound collections using a new
  computational compound–core relationship method." *ACS Omega* 4, 1027–1032 (the
  compound–core relationship, CCR); and Naveja, J. J. et al. (2019). "A general approach for
  retrosynthetic molecular core analysis." *J. Cheminf.* 11, 61.
- Naveja, J. J. and Medina-Franco, J. L. (2019). "Finding constellations in chemical space
  through core analysis." *Front. Chem.* 7, 510 (the plot).

**Cores.** A molecule's putative cores are the fragments retrosynthetic (RECAP) cuts leave that
keep at least two thirds of its heavy atoms, and the molecule itself. Attachment points become
hydrogens, so cores compare by canonical SMILES. The RECAP rules are RDKit's
(``rdkit.Chem.Recap.reactions``), with its rule that no cut may leave a fragment of three
carbons or fewer. The search here expands only fragments that are still large enough to be a
core (every fragment of a fragment is smaller), so a long peptide is quick where a full RECAP
hierarchy is not.

**Analog series** are the connected parts of the graph joining compounds to their cores.

**The plot** keeps the cores at least *min_compounds* compounds share; of cores that hold
exactly the same compounds, only the largest. Each is placed by its Morgan fingerprint (radius
2) on a t-SNE (or UMAP or PCA) map, sized by its compounds and joined to the cores it shares a
compound with.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Recap

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text

Cell = bytes | str | None
FRACTION = 2 / 3
#: fragments looked at for one molecule; past this its cores are the ones found so far
MAX_FRAGMENTS = 5000
#: what is left of a fragment without its attachment points that RECAP refuses to cut off
TINY = ("", "C", "CC", "CCC")
CHUNK_ROWS = 100
SUMMARIES = {
    "mean": "Mean",
    "median": "Median",
    "min": "Minimum",
    "max": "Maximum",
    "count": "Count",
}
MAPS = {"tsne": "t-SNE", "umap": "UMAP", "pca": "PCA"}


# -- the cores of one molecule -----------------------------------------------------------------
def _heavy(mol: Chem.Mol) -> int:
    """Heavy atoms, attachment points not counted."""
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def _capped(fragment: Chem.Mol) -> str | None:
    """Canonical SMILES of *fragment* with its attachment points made hydrogens."""
    edit = Chem.RWMol(fragment)
    for atom in edit.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)
            atom.SetIsotope(0)
            atom.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(edit.GetMol()))
    except (ValueError, RuntimeError):
        return None


def _cuts(mol: Chem.Mol) -> Iterator[list[tuple[str, Chem.Mol]]]:
    """Each single RECAP cut of *mol*: its pieces as ``(SMILES, piece)``."""
    for reaction in Recap.reactions:
        for products in reaction.RunReactants((mol,)):
            pieces = []
            for product in products:
                try:
                    Chem.SanitizeMol(product)
                except (ValueError, RuntimeError):
                    break
                smiles = Chem.MolToSmiles(product)
                if smiles.replace("*", "").replace("()", "") in TINY:
                    break
                pieces.append((smiles, product))
            else:
                yield pieces


def parent(mol: Chem.Mol) -> Chem.Mol:
    """The molecule cores are found for: its largest fragment (no salts), without hydrogens."""
    mol = Chem.RemoveHs(mol)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    return max(fragments, key=_heavy) if len(fragments) > 1 else mol


def putative_cores(mol: Chem.Mol, fraction: float = FRACTION) -> list[str]:
    """The molecule's cores (H-capped canonical SMILES), itself first."""
    mol = parent(mol)
    itself = Chem.MolToSmiles(mol)
    need = fraction * _heavy(mol) - 1e-9
    seen = {itself}
    cores = {itself: None}
    todo = [mol]
    with rdBase.BlockLogs():
        while todo and len(seen) < MAX_FRAGMENTS:
            for pieces in _cuts(todo.pop()):
                for smiles, piece in pieces:
                    if smiles in seen or _heavy(piece) < need:
                        continue
                    seen.add(smiles)
                    todo.append(piece)
                    core = _capped(piece)
                    if core is not None:
                        cores.setdefault(core, None)
    return list(cores)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def cores_chunk(args: tuple[float, Sequence[Cell]]) -> list[tuple[str, ...]]:
    """The cores of each cell of a chunk; () for no structure (a picklable worker)."""
    fraction, cells = args
    out: list[tuple[str, ...]] = []
    for cell in cells:
        mol = _mol(cell)
        try:
            out.append(tuple(putative_cores(mol, fraction)) if mol is not None else ())
        except (ValueError, RuntimeError):
            out.append(())
    return out


def find_cores(
    cells: Sequence[Cell],
    *,
    fraction: float = FRACTION,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> list[tuple[str, ...]]:
    """The cores of every cell, in order (``()`` where there is no structure)."""
    chunks = (
        (fraction, list(cells[start : start + CHUNK_ROWS]))
        for start in range(0, len(cells), CHUNK_ROWS)
    )
    out: list[tuple[str, ...]] = []
    for part in ordered_map(cores_chunk, chunks, executor, cancel=cancel):
        check(cancel)
        out.extend(part)
        if progress is not None:
            progress(len(out), len(cells))
    return out


# -- cores shared across a set -----------------------------------------------------------------
@dataclass
class Constellation:
    """The cores a set of compounds share, numbered 0…k-1."""

    cores: list[str]  # H-capped canonical SMILES
    members: list[np.ndarray]  # each core's compounds (indexes into the input), ascending
    heavy: np.ndarray  # heavy atoms of each core
    series: np.ndarray  # each core's analog series, 1 the largest
    series_size: dict[int, int]  # compounds in each series
    edges: np.ndarray  # k' × 2 pairs of cores sharing a compound
    compounds: int  # compounds with a structure

    def __len__(self) -> int:
        return len(self.cores)

    @property
    def n(self) -> np.ndarray:
        """Compounds holding each core."""
        return np.array([len(m) for m in self.members], dtype=np.int64)


def _components(n: int, groups: Sequence[Sequence[int]]) -> np.ndarray:
    """Connected parts of *n* items, items in a group joined (union-find): a label each."""
    root = list(range(n))

    def find(i: int) -> int:
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    for group in groups:
        first = find(group[0])
        for item in group[1:]:
            other = find(item)
            if other != first:
                root[other] = first
    return np.array([find(i) for i in range(n)], dtype=np.int64)


def constellation(
    core_lists: Sequence[Sequence[str]],
    *,
    min_compounds: int = 2,
    cancel: CancelToken | None = None,
) -> Constellation:
    """The shared cores of compounds whose cores are *core_lists* (from :func:`find_cores`)."""
    holders: dict[str, list[int]] = defaultdict(list)
    for compound, cores in enumerate(core_lists):
        for core in cores:
            holders[core].append(compound)
    check(cancel)
    shared = {core: rows for core, rows in holders.items() if len(rows) >= 2}
    labels = _components(len(core_lists), list(shared.values()))
    # of cores holding the same compounds, the one with the most heavy atoms
    best: dict[tuple[int, ...], tuple[int, str]] = {}
    for core, rows in shared.items():
        if len(rows) < max(min_compounds, 2):
            continue
        heavy = Chem.MolFromSmiles(core).GetNumHeavyAtoms()
        key = tuple(rows)
        if key not in best or (heavy, core) > best[key]:
            best[key] = (heavy, core)
    check(cancel)
    order = sorted(best.items(), key=lambda item: (-len(item[0]), -item[1][0], item[1][1]))
    cores = [core for _rows, (_heavy, core) in order]
    members = [np.array(rows, dtype=np.int64) for rows, _best in order]
    heavy = np.array([h for _rows, (h, _core) in order], dtype=np.int64)
    sizes = np.bincount(labels, minlength=len(core_lists)) if len(core_lists) else np.zeros(0)
    plotted = sorted({int(labels[m[0]]) for m in members}, key=lambda label: (-sizes[label], label))
    number = {label: i + 1 for i, label in enumerate(plotted)}
    series = np.array([number[int(labels[m[0]])] for m in members], dtype=np.int64)
    series_size = {number[label]: int(sizes[label]) for label in plotted}
    by_compound: dict[int, list[int]] = defaultdict(list)
    for index, rows in enumerate(members):
        for row in rows:
            by_compound[int(row)].append(index)
    pairs = {
        (a, b) for held in by_compound.values() for i, a in enumerate(held) for b in held[i + 1 :]
    }
    edges = np.array(sorted(pairs), dtype=np.int64).reshape(-1, 2)
    compounds = sum(1 for cores in core_lists if cores)
    return Constellation(cores, members, heavy, series, series_size, edges, compounds)


def summarise(con: Constellation, values: np.ndarray, how: str = "mean") -> np.ndarray:
    """One number per core from its compounds' *values* (NaN: none of them has one)."""
    if how not in SUMMARIES:
        raise ValueError(f"unknown summary {how!r}")
    values = np.asarray(values, dtype=np.float64)
    out = np.full(len(con), np.nan)
    for i, rows in enumerate(con.members):
        held = values[rows]
        held = held[np.isfinite(held)]
        if how == "count":
            out[i] = len(held)
        elif len(held):
            reduce = {"mean": np.mean, "median": np.median, "min": np.min, "max": np.max}[how]
            out[i] = float(reduce(held))
    return out


def place(
    con: Constellation,
    method: str = "tsne",
    *,
    seed: int = 42,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """k × 2 positions of the cores: a map of their Morgan (r=2) fingerprints. Under five cores
    the map is PCA; one or two cores sit in fixed places."""
    from chemsuite.analysis.embedding import Features, embed
    from chemsuite.analysis.fingerprints import compute_fingerprints

    if method not in MAPS:
        raise ValueError(f"unknown map {method!r}")
    k = len(con)
    if k <= 2:
        return np.array([[-1.0, 0.0], [1.0, 0.0]][:k]).reshape(k, 2) if k == 2 else np.zeros((k, 2))
    fps = compute_fingerprints(con.cores, "morgan2", cancel=cancel)
    features = Features.from_fingerprints(fps)
    if len(features) < 3:
        return _ring(k)
    chosen = method if len(features) >= 5 else "pca"
    coords = embed(features, chosen, seed=seed, progress=progress, cancel=cancel).coords
    if coords.shape[1] < 2:
        coords = np.column_stack([coords, np.zeros(len(coords))])
    missing = ~np.isfinite(coords[:, 0])
    coords[missing] = 0.0
    return np.ascontiguousarray(coords[:, :2], dtype=np.float64)


def _ring(k: int) -> np.ndarray:
    angles = np.linspace(0.0, 2.0 * math.pi, k, endpoint=False)
    return np.column_stack([np.cos(angles), np.sin(angles)])
