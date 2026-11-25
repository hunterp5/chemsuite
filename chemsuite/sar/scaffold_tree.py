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
"""The Scaffold Tree (Schuffenhauer et al., *J. Chem. Inf. Model.* 2007, 47, 47–58).

Each compound's Bemis–Murcko scaffold is pruned one ring at a time down to a single ring; at
each step the paper's prioritization rules choose which ring goes. The scaffolds then form a
tree: which ring systems a set holds and how they nest.

How each rule is computed follows ScaffoldGraph (Scott et al., *Bioinformatics* 2020, 36,
3930; MIT), which implements the paper's rules: ``prioritization/original_rules.py`` and the
ring removal of ``core/fragment.py``. Nothing is imported from it. Rule 7, which ScaffoldGraph
leaves out, is chemsuite's own reading of the paper's rule (see :func:`_rule_7`).

Every scaffold, a compound's own and each pruned one alike, goes through :func:`normalize`
(no stereo, sanitized, canonical SMILES), so the same scaffold is one node however it was
reached, on every platform.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field
from itertools import chain

import numpy as np
from rdkit import Chem, rdBase

from chemsuite.analysis.fingerprints import Cell, cell_mol
from chemsuite.core.jobs import CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.descriptors.scaffold import murcko

Progress = Callable[[int, int], None]
CHUNK_ROWS = 256  # rows per worker task
CHUNK_KEYS = 32  # scaffolds per worker task when finding parents
_ORIGIN = "scaffold_tree_atom"  # each atom's index in the child, kept through the edits
_LINKER = Chem.MolFromSmarts("*!@!=!#*")  # an acyclic single (or aromatic) bond: rule 3


def _hetero(atom: Chem.Atom) -> bool:
    return atom.GetAtomicNum() not in (1, 6)


def normalize(mol: Chem.Mol | None) -> str | None:
    """The canonical SMILES of *mol* without stereo, after a full sanitize (None when it has
    no ring or cannot be sanitized)."""
    if mol is None or not mol.GetNumAtoms():
        return None
    copy = Chem.Mol(mol)
    Chem.RemoveStereochemistry(copy)
    with rdBase.BlockLogs():
        again = Chem.MolFromSmiles(Chem.MolToSmiles(copy))
    if again is None or not again.GetRingInfo().NumRings():
        return None
    return Chem.MolToSmiles(again)


def scaffold_of(cell: Cell) -> str | None:
    """A cell's Murcko scaffold as a node key (None: blank, unreadable or no ring)."""
    mol = cell_mol(cell)
    if mol is None:
        return None
    with rdBase.BlockLogs():
        return normalize(murcko(mol))


@dataclass(frozen=True)
class Candidate:
    """A parent of *child*: *child* with ring *ring* (atom indices of the child) removed."""

    smiles: str
    mol: Chem.Mol  # the parent as sanitized, its atoms tagged with their child index
    ring: tuple[int, ...]
    bonds: tuple[int, ...]  # the removed ring's bonds in the child


def _fix(atom: Chem.Atom, hydrogen: bool) -> None:
    """An atom that lost a neighbor: free its hydrogens and stereo; an aromatic heteroatom
    takes one hydrogen when *hydrogen* (pyrrole-type N)."""
    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    if atom.GetIsAromatic() and atom.GetAtomicNum() != 6:
        atom.SetNumExplicitHs(1 if hydrogen else 0)
        atom.SetNoImplicit(True)
    else:
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)


def _remove(child: Chem.Mol, ring: int, hydrogen: bool) -> Chem.Mol | None:
    """*child* without ring *ring* (SSSR index), sanitized and cut back to its scaffold; None
    when that cannot be done."""
    info = child.GetRingInfo()
    atoms = info.AtomRings()[ring]
    bonds = info.BondRings()[ring]
    edit = Chem.RWMol(child)
    for atom in edit.GetAtoms():
        atom.SetIntProp(_ORIGIN, atom.GetIdx())
    gone = {i for i in atoms if info.NumAtomRings(i) == 1}
    touched: set[int] = set()
    # An atom of the removed ring double-bonded to a ring that stays (C1CCNCC1=C1C=CC=CC1)
    # stays as that ring's exocyclic =CH2, as a Murcko scaffold keeps exocyclic double bonds
    # (ScaffoldGraph does the same).
    for i in list(gone):
        partner = _exocyclic_partner(child, i, gone)
        if partner is None:
            continue
        gone.discard(i)
        touched.add(i)
        for bond in child.GetAtomWithIdx(i).GetBonds():
            other = bond.GetOtherAtomIdx(i)
            if other != partner and other not in gone:
                edit.RemoveBond(i, other)
    # Scheme 4 of the paper: a fused three-membered heterocycle leaves a double bond behind
    # (an epoxide or aziridine back to its alkene).
    shared = [b for b in bonds if info.NumBondRings(b) > 1]
    hetero = sum(_hetero(child.GetAtomWithIdx(i)) for i in atoms)
    if len(atoms) == 3 and hetero == 1 and len(shared) == 1:
        edit.GetBondWithIdx(shared[0]).SetBondType(Chem.BondType.DOUBLE)
        edit.GetBondWithIdx(shared[0]).SetIsAromatic(False)
    for b in bonds:  # a ring bond between two atoms that stay (bridged systems) is broken
        bond = child.GetBondWithIdx(b)
        a, c = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if info.NumBondRings(b) == 1 and a not in gone and c not in gone:
            edit.RemoveBond(a, c)
            touched.update((a, c))
    for i in gone:
        touched.update(n.GetIdx() for n in child.GetAtomWithIdx(i).GetNeighbors())
    for i in touched - gone:
        _fix(edit.GetAtomWithIdx(i), hydrogen)
    for i in sorted(gone, reverse=True):
        edit.RemoveAtom(i)
    mol = edit.GetMol()
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(mol)
    except (ValueError, RuntimeError):
        return None
    with rdBase.BlockLogs():
        pieces = Chem.GetMolFrags(murcko(mol), asMols=True, sanitizeFrags=False)
    for piece in pieces:
        Chem.FastFindRings(piece)
    ringed = [p for p in pieces if p.GetRingInfo().NumRings()]
    if len(ringed) != 1:
        return None
    return ringed[0]


def _exocyclic_partner(mol: Chem.Mol, i: int, gone: set[int]) -> int | None:
    """The ring atom that stays which atom *i* has an acyclic double bond to, if any."""
    for bond in mol.GetAtomWithIdx(i).GetBonds():
        other = bond.GetOtherAtom(mol.GetAtomWithIdx(i))
        if (
            bond.GetBondType() == Chem.BondType.DOUBLE
            and not bond.IsInRing()
            and other.IsInRing()
            and other.GetIdx() not in gone
        ):
            return other.GetIdx()
    return None


def parents(child: Chem.Mol) -> list[Candidate]:
    """Every scaffold *child* (a sanitized scaffold) can become by losing one ring: one per
    ring whose removal leaves one connected scaffold with one ring fewer."""
    info = child.GetRingInfo()
    rings = info.NumRings()
    out = []
    for r in range(rings):
        for hydrogen in (True, False):  # an aromatic heteroatom may need its H, or not
            mol = _remove(child, r, hydrogen)
            if mol is None:
                continue
            smiles = normalize(mol)
            if smiles is None:
                continue
            parent = Chem.MolFromSmiles(smiles)
            if parent.GetRingInfo().NumRings() != rings - 1:
                continue
            out.append(
                Candidate(smiles, mol, tuple(info.AtomRings()[r]), tuple(info.BondRings()[r]))
            )
            break
    return out


# -- the rules -------------------------------------------------------------------------------
# Each takes (child, candidate): a filter rule returns True for the candidates it prefers; a
# min / max rule returns the number to keep the smallest / largest of.


def _removed_atoms(child: Chem.Mol, c: Candidate) -> list[Chem.Atom]:
    return [child.GetAtomWithIdx(i) for i in c.ring]


def _rule_1(child, c) -> bool:
    """Remove heterocycles of size 3 first."""
    return len(c.ring) == 3 and sum(_hetero(a) for a in _removed_atoms(child, c)) == 1


def _rule_2(child, c) -> bool:
    """Do not remove rings of 12 or more atoms while smaller rings can go."""
    return len(c.ring) < 12


def _rule_3(child, c) -> int:
    """The parent with the fewest acyclic linker bonds (min)."""
    return len(c.mol.GetSubstructMatches(_LINKER))


def _delta(c: Candidate) -> int:
    """Δ = nrrb − (rings − 1) of the parent: nrrb counts the bonds rings share (a bond in k
    rings counts k − 1). Fused: 0; spiro or separate rings: below 0; bridged: above 0."""
    parent = Chem.MolFromSmiles(c.smiles)  # held: RingInfo dies with its molecule
    info = parent.GetRingInfo()
    bonds = list(chain(*info.BondRings()))
    nrrb = len(bonds) - len(set(bonds))
    return nrrb - (info.NumRings() - 1)


def _rule_4(child, c) -> int:
    """Keep bridged, spiro and non-linear fused systems: the largest |Δ| (max)."""
    return abs(_delta(c))


def _rule_5(child, c) -> bool:
    """Bridged systems before spiro: Δ ≥ 1."""
    return _delta(c) >= 1


def _rule_6(child, c) -> bool:
    """Remove rings of 3, 5 and 6 atoms first."""
    return len(c.ring) in (3, 5, 6)


def _rule_7(child, c) -> bool:
    """A fully aromatic ring system must not be dissected so that what is left of it is not
    aromatic any more. chemsuite's reading: when the removed ring is part of a ring system
    whose atoms are all aromatic, every atom of that system the parent keeps stays aromatic."""
    info = child.GetRingInfo()
    system = set(c.ring)
    grown = True
    while grown:  # the removed ring's fused system: rings sharing a bond with it, and so on
        grown = False
        for atoms, bonds in zip(info.AtomRings(), info.BondRings(), strict=True):
            if not system.issuperset(atoms) and _shares_bond(child, system, bonds):
                system.update(atoms)
                grown = True
    if not all(child.GetAtomWithIdx(i).GetIsAromatic() for i in system):
        return True
    kept = {
        a.GetIntProp(_ORIGIN): a.GetIsAromatic() for a in c.mol.GetAtoms() if a.HasProp(_ORIGIN)
    }
    if not kept:  # the tags were lost on the way: the rule cannot be judged
        raise RuntimeError("ring removal lost its atom tags")
    return all(kept[i] for i in system if i in kept)


def _shares_bond(mol: Chem.Mol, atoms: set[int], bonds: Sequence[int]) -> bool:
    for b in bonds:
        bond = mol.GetBondWithIdx(b)
        if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms:
            return True
    return False


def _rule_8(child, c) -> int:
    """Remove rings with the fewest heteroatoms first (min)."""
    return sum(_hetero(a) for a in _removed_atoms(child, c))


def _element_rule(number: int, doc: str):
    def rule(child, c) -> int:
        return sum(a.GetAtomicNum() == number for a in _removed_atoms(child, c))

    rule.__doc__ = doc
    return rule


def _rule_10(child, c) -> int:
    """Remove smaller rings first (min)."""
    return len(c.ring)


def _rule_11(child, c) -> bool:
    """Keep non-aromatic rings: remove fully aromatic rings first."""
    return all(child.GetBondWithIdx(b).GetIsAromatic() for b in c.bonds)


def _rule_12(child, c) -> bool:
    """Remove rings first where a linker is attached to a ring heteroatom at either end."""
    ring = set(c.ring)
    ends: set[int] = set()
    for start in ring:
        atom = child.GetAtomWithIdx(start)
        for bond in atom.GetBonds():
            if bond.IsInRing():
                continue
            far = _linker_end(child, start, bond.GetOtherAtomIdx(start))
            if far is not None:
                ends.update((start, far))
    return any(_hetero(child.GetAtomWithIdx(i)) for i in ends)


def _linker_end(mol: Chem.Mol, start: int, first: int) -> int | None:
    """The ring atom a linker leaving *start* through *first* reaches (None: a dead end,
    such as an exocyclic =O)."""
    seen, frontier = {start}, [first]
    while frontier:
        i = frontier.pop()
        if i in seen:
            continue
        seen.add(i)
        atom = mol.GetAtomWithIdx(i)
        if atom.IsInRing():
            return i
        frontier.extend(n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() not in seen)
    return None


#: (name, rule, kind): kind "filter" keeps True, "min" / "max" keep the extreme
RULES: tuple[tuple[str, Callable, str], ...] = (
    ("1", _rule_1, "filter"),
    ("2", _rule_2, "filter"),
    ("3", _rule_3, "min"),
    ("4", _rule_4, "max"),
    ("5", _rule_5, "filter"),
    ("6", _rule_6, "filter"),
    ("7", _rule_7, "filter"),
    ("8", _rule_8, "min"),
    ("9a", _element_rule(7, "Fewest nitrogens in the removed ring (min)."), "min"),
    ("9b", _element_rule(8, "Then fewest oxygens (min)."), "min"),
    ("9c", _element_rule(16, "Then fewest sulfurs (min)."), "min"),
    ("10", _rule_10, "min"),
    ("11", _rule_11, "filter"),
    ("12", _rule_12, "filter"),
)


#: what each rule (and "only") means, for the panel's tooltips
MEANINGS = {
    "only": "every ring removal gives the same parent",
    "1": "a three-membered heterocycle goes first",
    "2": "a ring of 12 or more atoms stays while smaller rings can go",
    "3": "the parent with the fewest linker bonds",
    "4": "a bridged, spiro or non-linear fused system is kept",
    "5": "bridged systems are kept before spiro ones",
    "6": "rings of 3, 5 or 6 atoms go first",
    "7": "the aromatic system left behind stays aromatic",
    "8": "the ring with the fewest heteroatoms goes",
    "9a": "the ring with the fewest nitrogens goes",
    "9b": "the ring with the fewest oxygens goes",
    "9c": "the ring with the fewest sulfurs goes",
    "10": "the smallest ring goes",
    "11": "an aromatic ring goes before a non-aromatic one",
    "12": "a ring whose linker ends on a ring heteroatom goes first",
    "13": "the parent whose SMILES sorts first",
}


def choose(child: Chem.Mol, candidates: Sequence[Candidate]) -> tuple[Candidate, str]:
    """The parent the rules pick, and the rule that decided ("only" for a single candidate,
    "13" for the tie-break: the smallest canonical SMILES). Each rule narrows the candidates
    left; a rule no candidate passes is skipped."""
    if not candidates:
        raise ValueError("no candidate parents")
    remaining = list(candidates)
    if len({c.smiles for c in remaining}) == 1:
        return remaining[0], "only"
    for name, rule, kind in RULES:
        if kind == "filter":
            kept = [c for c in remaining if rule(child, c)]
        else:
            values = [rule(child, c) for c in remaining]
            best = min(values) if kind == "min" else max(values)
            kept = [c for c, v in zip(remaining, values, strict=True) if v == best]
        if kept:
            remaining = kept
        if len({c.smiles for c in remaining}) == 1:
            return remaining[0], name
    return min(remaining, key=lambda c: c.smiles), "13"


def parent_of(smiles: str) -> tuple[str, str] | None:
    """``(parent SMILES, deciding rule)`` of a node, or None for a root."""
    child = Chem.MolFromSmiles(smiles)
    if child is None or child.GetRingInfo().NumRings() < 2:
        return None
    candidates = parents(child)
    if not candidates:
        return None
    chosen, rule = choose(child, candidates)
    return chosen.smiles, rule


# -- the tree --------------------------------------------------------------------------------
@dataclass
class ScaffoldTree:
    """Scaffolds as nodes (``smiles``, ``rings``, ``parent`` index or −1, the ``rule`` that
    pruned each to its parent) and each input row's node (−1: no ring, blank or unreadable).
    """

    smiles: list[str]
    rings: np.ndarray
    parent: np.ndarray
    rule: list[str]
    node: np.ndarray  # per input row
    _children: list[list[int]] = field(default_factory=list, repr=False)
    _members: list[np.ndarray] | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self._children = [[] for _ in self.smiles]
        for i, p in enumerate(self.parent.tolist()):
            if p >= 0:
                self._children[p].append(i)

    def __len__(self) -> int:
        return len(self.smiles)

    def roots(self) -> list[int]:
        return [i for i, p in enumerate(self.parent.tolist()) if p < 0]

    def children(self, node: int) -> list[int]:
        return list(self._children[node])

    def own_counts(self) -> np.ndarray:
        """Rows whose own scaffold is each node."""
        placed = self.node[self.node >= 0]
        return np.bincount(placed, minlength=len(self)) if len(self) else np.zeros(0, int)

    def subtree_counts(self) -> np.ndarray:
        """Rows whose scaffold is each node or built on it."""
        counts = self.own_counts().astype(np.int64)
        for i in self._postorder():
            p = int(self.parent[i])
            if p >= 0:
                counts[p] += counts[i]
        return counts

    def _postorder(self) -> list[int]:
        order, stack = [], list(self.roots())
        while stack:
            i = stack.pop()
            order.append(i)
            stack.extend(self._children[i])
        return order[::-1]  # children before their parents

    def members(self) -> list[np.ndarray]:
        """Per node, the input rows whose scaffold is it or built on it (ascending)."""
        if self._members is None:
            parent = self.parent.tolist()
            lists: list[list[int]] = [[] for _ in self.smiles]
            for row, at in enumerate(self.node.tolist()):
                while at >= 0:  # a row belongs to its node and every node above it
                    lists[at].append(row)
                    at = parent[at]
            self._members = [np.array(rows, dtype=np.int64) for rows in lists]
        return self._members

    def subtree_rows(self, node: int) -> np.ndarray:
        """Input rows whose scaffold is *node* or built on it."""
        return self.members()[node]

    def depth(self, node: int) -> int:
        d = 0
        while self.parent[node] >= 0:
            node = int(self.parent[node])
            d += 1
        return d


def _ring_count(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)  # held: RingInfo dies with its molecule
    return mol.GetRingInfo().NumRings()


def scaffold_chunk(cells: Sequence[Cell]) -> list[str | None]:
    """:func:`scaffold_of` for each cell (a worker's share)."""
    return [scaffold_of(cell) for cell in cells]


def parent_chunk(keys: Sequence[str]) -> list[tuple[str, str, str, int] | None]:
    """For each node key: ``(key, parent, rule, rings)``, or ``(key, "", "", rings)`` for a
    root (a worker's share)."""
    out = []
    for key in keys:
        found = parent_of(key)
        up, why = found if found is not None else ("", "")
        out.append((key, up, why, _ring_count(key)))
    return out


def _chunks(items: Sequence, size: int):
    for start in range(0, len(items), size):
        yield items[start : start + size]


def build(
    cells: Sequence[Cell],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> ScaffoldTree:
    """The scaffold tree of *cells*: each row's scaffold, pruned to its root. With
    *executor*, the scaffolds and then each level of parents are found on its workers."""
    total = len(cells) * 5  # progress units: the scaffolds take the first fifth
    keys: list[str | None] = []
    for part in ordered_map(scaffold_chunk, _chunks(cells, CHUNK_ROWS), executor, cancel=cancel):
        keys.extend(part)
        if progress is not None:  # the scaffolds are about a fifth of the work
            progress(len(keys), total)
    index: dict[str, int] = {}
    smiles: list[str] = []
    parent: list[int] = []
    rule: list[str] = []
    rings: list[int] = []

    def add(key: str) -> int:
        at = index.get(key)
        if at is None:
            at = index[key] = len(smiles)
            smiles.append(key)
            parent.append(-1)
            rule.append("")
            rings.append(0)
        return at

    frontier = list(dict.fromkeys(k for k in keys if k is not None))
    for key in frontier:
        add(key)
    # the first level of parents is most of the work (three fifths); each later level, with
    # fewer and smaller scaffolds, takes half of what is left
    start, span = len(cells), len(cells) * 3
    while frontier:  # one level of parents at a time; new parents make the next level
        found = []
        for part in ordered_map(
            parent_chunk, _chunks(frontier, CHUNK_KEYS), executor, cancel=cancel
        ):
            found.extend(part)
            if progress is not None:
                progress(start + span * len(found) // len(frontier), total)
        start += span
        span = (total - start) // 2
        frontier = []
        for key, up, why, count in found:
            at = index[key]
            rings[at] = count
            if up:
                known = up in index
                parent[at], rule[at] = add(up), why
                if not known:
                    frontier.append(up)
    if progress is not None:
        progress(total, total)
    node = np.array([index[k] if k is not None else -1 for k in keys], dtype=np.int64)
    return ScaffoldTree(
        smiles, np.array(rings, dtype=np.int64), np.array(parent, dtype=np.int64), rule, node
    )
