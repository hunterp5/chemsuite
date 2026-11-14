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

"""Building a library from a core and a list of substituents per position.

The core carries numbered attachment points (``[*:1]``, ``[*:2]``…), as
:mod:`chemsuite.sar.rgroups` and :mod:`chemsuite.sar.central_ring` write them, so a
decomposition feeds straight back in: take the substituents a series already uses and make
every combination of them, including the ones nobody has made.

A substituent is written as R-Group Decomposition labels it (``*Cl``, ``*OC``), or ``H`` to
leave the position unsubstituted. Each combination is joined with RDKit's ``molzip``, which
bonds the dummy atoms that share an attachment number.

The result is :class:`chemsuite.analysis.reactions.Library`, the same one reaction enumeration
gives, so both fill a new window the same way.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from concurrent.futures import Executor

import numpy as np
from rdkit import Chem, rdBase

from chemsuite.analysis.reactions import Library, combinations_count
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map

Progress = Callable[[int, int], None]

CHUNK = 2000
INLINE_CHUNK = 200
POOL_MIN_COMBINATIONS = 5000
#: a position left unsubstituted
HYDROGEN = "H"
_POSITION = re.compile(r"\[\*:(\d+)\]")
_FAILURES = (ValueError, RuntimeError)


class CoreError(ValueError):
    """The core cannot be used to build a library."""


def core_positions(core: str) -> tuple[int, ...]:
    """The attachment numbers of *core*, in order: ``[*:1]`` and ``[*:2]`` give (1, 2)."""
    return tuple(sorted({int(n) for n in _POSITION.findall(core)}))


def parse_core(core: str) -> Chem.Mol:
    """*core* as a molecule with numbered attachment points, checked."""
    text = core.strip()
    if not text:
        raise CoreError("Enter a core with numbered attachment points, such as c1ccc([*:1])cc1.")
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise CoreError(f"Not a valid core: {core}")
    if not core_positions(text):
        raise CoreError("The core has no numbered attachment point ([*:1], [*:2] …).")
    return mol


def is_hydrogen(text: str) -> bool:
    """Whether *text* means "leave this position unsubstituted"."""
    stripped = text.strip()
    return not stripped or stripped.upper() == HYDROGEN


def parse_group(text: str, position: int) -> Chem.Mol | None:
    """One substituent, with its attachment numbered *position*.

    None means it cannot be used — either because it is :func:`is_hydrogen` (no substituent)
    or because it is not one attachment point on a valid structure. The two are told apart by
    asking :func:`is_hydrogen` first.
    """
    if is_hydrogen(text):
        return None
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(text.strip())
    if mol is None:
        return None
    dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 1:
        return None  # a substituent attaches at exactly one point
    dummies[0].SetAtomMapNum(position)
    return mol


def _without(core: Chem.Mol, positions: Sequence[int]) -> Chem.Mol | None:
    """*core* with the attachment points of *positions* removed, leaving hydrogens."""
    work = Chem.RWMol(core)
    drop = [
        a.GetIdx()
        for a in work.GetAtoms()
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() in set(positions)
    ]
    for index in sorted(drop, reverse=True):
        work.RemoveAtom(index)
    out = work.GetMol()
    try:
        Chem.SanitizeMol(out)
    except _FAILURES:
        return None
    return out


def build(core: Chem.Mol, groups: Sequence[Chem.Mol | None], positions: Sequence[int]) -> str:
    """One product: *core* with ``groups[i]`` at ``positions[i]``. "" when it cannot be made.

    A None group leaves its position unsubstituted.
    """
    with rdBase.BlockLogs():
        blank = [p for p, g in zip(positions, groups, strict=True) if g is None]
        work = _without(core, blank) if blank else Chem.Mol(core)
        if work is None:
            return ""
        for group in groups:
            if group is not None:
                work = Chem.CombineMols(work, group)
        try:
            out = Chem.molzip(work)
            Chem.SanitizeMol(out)
        except _FAILURES:
            return ""
        if any(a.GetAtomicNum() == 0 for a in out.GetAtoms()):
            return ""  # an attachment point was left open
        try:
            return Chem.MolToSmiles(out)
        except _FAILURES:
            return ""


def enumerate_chunk(args: tuple[str, list[list[str]], int, int]) -> Library:
    """Combinations ``start`` to ``stop`` (flat indices, last list fastest)."""
    core_text, groups, start, stop = args
    core = parse_core(core_text)
    positions = core_positions(core_text)
    # each entry is (usable, molecule): a molecule, or None for an unsubstituted position;
    # a substituent that is neither makes its combinations fail rather than passing for H
    parsed: list[list[tuple[bool, Chem.Mol | None]]] = [
        [
            (True, None)
            if is_hydrogen(text)
            else ((g := parse_group(text, positions[i])) is not None, g)
            for text in column
        ]
        for i, column in enumerate(groups)
    ]
    sizes = [len(column) for column in groups]
    out = Library([[] for _ in groups])
    for flat in range(start, stop):
        index = np.unravel_index(flat, sizes)
        out.combinations += 1
        chosen = [parsed[i][int(k)] for i, k in enumerate(index)]
        if not all(usable for usable, _mol in chosen):
            out.failed += 1
            continue
        product = build(core, [mol for _usable, mol in chosen], positions)
        if not product:
            out.failed += 1
            continue
        for i, k in enumerate(index):
            out.blocks[i].append(int(k))
        out.products.append(product)
    return out


def enumerate_rgroups(
    core: str,
    groups: Sequence[Sequence[str]],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Library:
    """Every combination of one substituent per position.

    ``groups[i]`` holds the substituents for the *i*-th attachment point of *core*, in the
    order :func:`core_positions` gives. Cancelled: ``partial`` holds the products already made.
    """
    positions = core_positions(core)
    parse_core(core)
    if len(groups) != len(positions):
        raise CoreError(
            f"The core has {len(positions)} attachment point(s); choose substituents for each."
        )
    lists = [list(column) for column in groups]
    if any(not column for column in lists):
        raise CoreError("Every position needs at least one substituent.")
    total = combinations_count([len(column) for column in lists])
    size = CHUNK if executor is not None else INLINE_CHUNK
    chunks = ((core, lists, start, min(start + size, total)) for start in range(0, total, size))
    out = Library([[] for _ in lists])
    try:
        for part in ordered_map(enumerate_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(out.combinations, total)
    except Cancelled:
        raise Cancelled(out if out.combinations else None) from None
    return out


__all__ = [
    "is_hydrogen",
    "HYDROGEN",
    "POOL_MIN_COMBINATIONS",
    "CoreError",
    "build",
    "core_positions",
    "enumerate_chunk",
    "enumerate_rgroups",
    "parse_core",
    "parse_group",
]
