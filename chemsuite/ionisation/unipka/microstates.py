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

"""A molecule's protonation microstates, by Uni-pKa's iterated template matching.

A port of Uni-pKa's enumerator (``enumerator/main.py``,
https://github.com/dptech-corp/Uni-pKa, Apache-2.0; see ``data/LICENSE-Uni-pKa.txt``). Its
templates (``data/unipka_sites.tsv``, upstream's ``simple_smarts_pattern.tsv``) extend
MolGpKa's (MIT, ``data/LICENSE-MolGpKa.txt``). Each acid template names the atom that loses a
proton, each base template the atom that gains one. From the input, every acid site is
deprotonated (and every base site protonated), the results are filtered, and the pools grow
until they stop changing; one macrostate (total charge) leads to the next.

Upstream works on SMILES; here every heavy atom carries its index + 1 as an atom-map number,
so each microstate maps back to the input's atoms (:func:`ensemble` returns plain per-atom
hydrogen counts and charges). Map numbers do not change what a SMARTS matches, nor which
microstates are equal, so the pools are upstream's.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from importlib import resources

from rdkit import Chem, RDLogger

ACID, BASE = "A", "B"
MAX_ITER = 10  # upstream's command-line default
LOWER, UPPER = -2, 2  # total charges enumerated, as upstream's defaults

#: structures upstream filters out as unreasonable
_FILTERS = (
    "[#6X5]",
    "[#7X5]",
    "[#8X4]",
    "[*r]=[*r]=[*r]",
    "[#1]-[*+1]~[*-1]",
    "[#1]-[*+1]=,:[*]-,:[*-1]",
    "[#1]-[*+1]-,:[*]=,:[*-1]",
    "[*+2]",
    "[*-2]",
    "[#1]-[#8+1].[#8-1,#7-1,#6-1]",
    "[#1]-[#7+1,#8+1].[#7-1,#6-1]",
    "[#1]-[#8+1].[#8-1,#6-1]",
    "[#1]-[#7+1].[#8-1]-[C](-[C,#1])(-[C,#1])",
    "[OX1]=[C]-[OH2+1]",
    "[NX1,NX2H1,NX3H2]=[C]-[O]-[H]",
    "[#6-1]=[*]-[*]",
    "[cX2-1]",
    "[N+1](=O)-[O]-[H]",
)


@dataclass(frozen=True)
class Microstate:
    """Hydrogens and formal charge on each heavy atom of the input, in its order."""

    hydrogens: tuple[int, ...]
    charges: tuple[int, ...]
    smiles: str  # canonical, without map numbers

    @property
    def charge(self) -> int:
        return sum(self.charges)


@cache
def _templates() -> dict[str, tuple[tuple[Chem.Mol, int], ...]]:
    text = resources.files(__package__).joinpath("data/unipka_sites.tsv").read_text("utf-8")
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))[1:]
    out: dict[str, list[tuple[Chem.Mol, int]]] = {ACID: [], BASE: []}
    for _name, smarts, index, kind in rows:
        out[kind].append((Chem.MolFromSmarts(smarts), int(index)))
    return {k: tuple(v) for k, v in out.items()}


@cache
def _filters() -> tuple[Chem.Mol, ...]:
    return tuple(Chem.MolFromSmarts(s) for s in _FILTERS)


def _canon(smiles: str, chiral: bool = True) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=chiral)


def _sites(mol: Chem.Mol, kind: str) -> list[int]:
    """The atoms (of *mol* with explicit hydrogens) the templates of *kind* act on."""
    found = set()
    for pattern, index in _templates()[kind]:
        for match in mol.GetSubstructMatches(pattern):
            found.add(match[index])
    return sorted(found)


def _protonated(mol: Chem.Mol, idx: int, kind: str) -> str | None:
    """*mol* (hydrogens implicit) with a proton taken from (ACID) or given to (BASE) atom *idx*,
    as canonical SMILES; None when RDKit cannot sanitise it. The templates name heavy atoms,
    whose indices are the same with hydrogens added."""
    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(idx)
    if kind == ACID:
        atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
        atom.UpdatePropertyCache()
    else:
        atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
    try:
        out = rw.GetMol()
        Chem.SanitizeMol(out)
        out = Chem.MolFromSmiles(Chem.MolToSmiles(out))
        return None if out is None else _canon(Chem.MolToSmiles(Chem.RemoveHs(out)))
    except (ValueError, RuntimeError):
        return None


def _step(smiles: str, kind: str) -> list[str]:
    """Every microstate one proton away from *smiles* by the templates of *kind*."""
    mol = Chem.MolFromSmiles(smiles)  # as upstream: matched with hydrogens, changed without
    out = {_protonated(mol, idx, kind) for idx in _sites(Chem.AddHs(mol), kind)}
    return sorted(s for s in out if s is not None)


def _reasonable(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if any(mol.HasSubstructMatch(p) for p in _filters()):
        return False
    try:
        Chem.SanitizeMol(mol)
    except (ValueError, RuntimeError):
        return False
    return True


def _stereo_count(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return sum(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())


def _filtered(pool: list[str]) -> list[str]:
    """Upstream's filters: reasonable structures, and of those equal but for stereo, the one
    that keeps the most stereocentres (the first of equals)."""
    best: dict[str, tuple[str, int]] = {}
    for smiles in (s for s in pool if _reasonable(s)):
        flat = _canon(smiles, chiral=False)
        count = _stereo_count(smiles)
        if flat not in best or count > best[flat][1]:
            best[flat] = (smiles, count)
    return [s for s, _n in best.values()]


def _grow(start: list[str], mode: str, max_iter: int = MAX_ITER) -> tuple[list[str], list[str]]:
    """Upstream's ``enumerate_template``: the acid and base pools grown from *start* (the acid
    pool with *mode* ACID, the base pool with BASE), alternately, until neither changes."""
    acids, bases = (list(start), []) if mode == ACID else ([], list(start))
    sizes = (-1, -1)
    i = 0
    while (len(acids), len(bases)) != sizes and i < max_iter:
        sizes = (len(acids), len(bases))
        if (mode == ACID) == (i % 2 == 0):
            for smiles in list(acids):
                bases += _filtered(_step(smiles, ACID))
            acids += _filtered(acids)
        else:
            for smiles in list(bases):
                acids += _filtered(_step(smiles, BASE))
            bases += _filtered(bases)
        acids = sorted(set(_filtered(acids)))
        bases = sorted(set(_filtered(bases)))
        i += 1
    return [_canon(s) for s in acids], [_canon(s) for s in bases]


def _ensemble_smiles(smiles: str, lower: int, upper: int) -> dict[int, list[str]]:
    """Upstream's ``enum_ensemble`` for one molecule: total charge → its microstates."""
    q0 = Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))
    out: dict[int, list[str]] = {q0: [smiles]}
    middle = [smiles]
    if q0 > lower:
        middle, below = _grow(middle, ACID)
        if below:
            out[q0 - 1] = below
        for q in range(q0 - 2, lower, -1):
            if q + 1 in out:
                _a, found = _grow(out[q + 1], ACID)
                if found:
                    out[q] = found
    if q0 < upper:
        above, middle = _grow(middle, BASE)
        if above:
            out[q0 + 1] = above
        for q in range(q0 + 2, upper):
            if q - 1 in out:
                found, _b = _grow(out[q - 1], BASE)
                if found:
                    out[q] = found
    out[q0] = middle
    return out


def _microstate(smiles: str, n: int) -> Microstate:
    mol = Chem.MolFromSmiles(smiles)
    hydrogens, charges = [0] * n, [0] * n
    for atom in mol.GetAtoms():
        k = atom.GetAtomMapNum() - 1
        if 0 <= k < n:
            hydrogens[k] = atom.GetTotalNumHs()
            charges[k] = atom.GetFormalCharge()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Microstate(tuple(hydrogens), tuple(charges), Chem.MolToSmiles(mol))


def ensemble(mol: Chem.Mol, lower: int = LOWER, upper: int = UPPER) -> dict[int, list[Microstate]]:
    """The protonation microstates of *mol* (heavy atoms, its hydrogens implicit or counted),
    by total charge from *lower* to *upper*, each as hydrogens and charges per heavy atom of
    *mol* (those of *mol* itself among them)."""
    RDLogger.DisableLog("rdApp.*")
    try:
        heavy = Chem.RemoveHs(mol)
        mapped = Chem.Mol(heavy)
        for atom in mapped.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        found = _ensemble_smiles(_canon(Chem.MolToSmiles(mapped)), lower, upper)
    finally:
        RDLogger.EnableLog("rdApp.*")
    n = heavy.GetNumAtoms()
    return {q: [_microstate(s, n) for s in pool] for q, pool in sorted(found.items())
            if lower <= q <= upper}  # fmt: skip
