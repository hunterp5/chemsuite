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
"""R-group decomposition on a shared core (RDKit), with readable substituent labels.

A core is SMILES or SMARTS; label attachment points ``[*:1]``, ``[*:2]``… to fix which
position is R1, R2… (otherwise RDKit numbers every open position it finds). A substituent is
labelled by its canonical SMILES with ``*`` at the attachment (``*Cl``, ``*OC``), and ``H``
for hydrogen, so the same group always gets the same label.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdRGroupDecomposition as rgd

from chemsuite.core.jobs import CancelToken, check

HYDROGEN = "H"
_POSITION = re.compile(r"^R(\d+)$")
_MAP = re.compile(r"\[\*:(\d+)\]")
CHUNK = 500


def parse_core(text: str) -> Chem.Mol:
    """The core as a query: SMARTS when it parses as SMARTS with query features, else SMILES."""
    text = text.strip()
    if not text:
        raise ValueError("Enter a core (SMILES or SMARTS), or suggest one from the selection.")
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            mol = Chem.MolFromSmarts(text)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError(f"Not a valid core: {text}")
    return mol


def has_labelled_positions(core: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() == 0 and a.GetAtomMapNum() > 0 for a in core.GetAtoms())


def suggest_core(mols: Sequence[Chem.Mol], *, timeout: int = 10) -> str:
    """The maximum common substructure of *mols* (rings kept whole), as SMARTS."""
    mols = [m for m in mols if m is not None]
    if len(mols) < 2:
        raise ValueError("Select at least two structures to suggest a core.")
    params = rdFMCS.MCSParameters()
    params.Timeout = timeout
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    result = rdFMCS.FindMCS(mols, params)
    if not result.numAtoms:
        raise ValueError("These structures share no common core.")
    return result.smartsString


def substituent_label(smiles: str) -> str:
    """``C[*:1]`` → ``*C``; ``[H][*:1]`` → ``H``. Several attachments keep their numbers."""
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return smiles
    dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if all(a.GetAtomicNum() in (0, 1) for a in mol.GetAtoms()) and len(dummies) == 1:
        return HYDROGEN
    if len(dummies) == 1:
        dummies[0].SetAtomMapNum(0)
    with rdBase.BlockLogs():
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
        except (ValueError, RuntimeError):
            return smiles
        return Chem.MolToSmiles(mol)


@dataclass
class Decomposition:
    """Per input row: whether it matched, and each position's substituent label."""

    positions: tuple[str, ...]
    labels: dict[str, list[str]]  # position → label per row ("" where unmatched)
    matched: np.ndarray  # bool per row
    core: str  # the core with numbered attachment points, as RDKit reported it
    fragments: dict[tuple[str, str], str] = field(default_factory=dict)  # (pos, label) → SMILES

    def __len__(self) -> int:
        return len(self.matched)

    def varying(self) -> tuple[str, ...]:
        """Positions with more than one substituent among matched rows."""
        out = []
        for pos in self.positions:
            values = {v for v, ok in zip(self.labels[pos], self.matched, strict=True) if ok}
            if len(values) > 1:
                out.append(pos)
        return tuple(out)

    def counts(self, position: str) -> Counter[str]:
        return Counter(v for v, ok in zip(self.labels[position], self.matched, strict=True) if ok)


def _params(core: Chem.Mol, exhaustive: bool) -> rgd.RGroupDecompositionParameters:
    params = rgd.RGroupDecompositionParameters()
    params.removeHydrogensPostMatch = True
    params.onlyMatchAtRGroups = has_labelled_positions(core)
    params.matchingStrategy = (
        rgd.RGroupMatching.Exhaustive if exhaustive else rgd.RGroupMatching.GreedyChunks
    )
    params.timeout = 60.0
    return params


def _position_key(name: str) -> int:
    match = _POSITION.match(name)
    return int(match.group(1)) if match else 10**6


def decompose(
    mols: Sequence[Chem.Mol | None],
    core: Chem.Mol,
    *,
    exhaustive: bool = False,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> Decomposition:
    """Decompose every molecule onto *core*.

    With numbered attachment points the work runs in chunks (each position keeps its number);
    without, RDKit must see every molecule at once to number positions consistently.
    """
    n = len(mols)
    labelled = has_labelled_positions(core)
    size = CHUNK if labelled else max(n, 1)
    rows: list[dict | None] = [None] * n
    cores: Counter[str] = Counter()
    for start in range(0, n, size):
        check(cancel)
        index = [i for i in range(start, min(n, start + size)) if mols[i] is not None]
        batch = []
        for i in index:
            mol = mols[i]
            if mol.GetNumConformers() and mol.GetConformer().Is3D():
                mol = Chem.RemoveHs(mol)
            batch.append(mol)
        if batch:
            with rdBase.BlockLogs():
                found, unmatched = rgd.RGroupDecompose(
                    [core], batch, asSmiles=True, asRows=True, options=_params(core, exhaustive)
                )
            missing = set(unmatched)
            matched_index = [i for k, i in enumerate(index) if k not in missing]
            for i, row in zip(matched_index, found, strict=True):
                rows[i] = row
                cores[row.get("Core", "")] += 1
        if progress is not None:
            progress(min(n, start + size), n)
    positions = sorted(
        {k for row in rows if row for k in row if _POSITION.match(k)}, key=_position_key
    )
    labels: dict[str, list[str]] = {p: [""] * n for p in positions}
    fragments: dict[tuple[str, str], str] = {}
    cache: dict[str, str] = {}
    for i, row in enumerate(rows):
        if row is None:
            continue
        for pos in positions:
            smiles = row.get(pos, "")
            if not smiles:
                label = HYDROGEN  # RDKit leaves out a position that is not substituted
                smiles = f"[H][*:{_position_key(pos)}]"
            else:
                label = cache.get(smiles) or cache.setdefault(smiles, substituent_label(smiles))
            labels[pos][i] = label
            fragments.setdefault((pos, label), smiles)
    matched = np.array([row is not None for row in rows], dtype=np.bool_)
    core_text = cores.most_common(1)[0][0] if cores else ""
    return Decomposition(tuple(positions), labels, matched, core_text, fragments)


def assemble(core: str, fragments: dict[str, str]) -> Chem.Mol | None:
    """Join a numbered core with one fragment SMILES per position (``[*:n]`` ↔ ``[*:n]``)."""
    with rdBase.BlockLogs():
        parts = [core, *fragments.values()]
        mol = Chem.MolFromSmiles(".".join(parts))
        if mol is None:
            return None
        try:
            mol = Chem.molzip(mol)
            mol = Chem.RemoveHs(mol)
            Chem.SanitizeMol(mol)
        except (ValueError, RuntimeError):
            return None
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return None
    return mol


def core_positions(core: str) -> list[int]:
    return sorted({int(m) for m in _MAP.findall(core)})
