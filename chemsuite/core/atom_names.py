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

"""Atoms named in a cell: ``C7``, or ``C7 0.81, N3 0.05, C12 0.04`` (Qt-free).

Columns that point at particular atoms — the sites of metabolism a P450 isoform attacks, and
anything written the same way — name each atom by its element and its index in the molecule,
best first, optionally with a score. This module reads that text back into atom indices so a
structure column's drawings can show the atoms, and nothing else in the app has to know which
tool wrote the column.

The element is checked against the molecule's atom at that index before anything is drawn, so a
column of indices for *another* molecule, or a column that merely looks like atom names,
highlights nothing rather than the wrong atoms.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

from rdkit import Chem

#: One named atom: an element symbol, an index, and an optional score after a space.
_ATOM = re.compile(r"^([A-Z][a-z]?)(\d+)(?:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))?$")
#: How many of a column's values are read when deciding whether it names atoms.
JUDGED = 200


@dataclass(frozen=True)
class AtomName:
    """``C7 0.81``: the atom at index 7, which is carbon, scored 0.81."""

    element: str
    atom: int
    score: float | None = None


def read_atom_names(text: str) -> list[AtomName] | None:
    """The atoms *text* names, in the order it names them; None when it does not read as atom
    names at all.

    Blank text names no atoms, which is not a failure: a column may be filled for some rows
    only.
    """
    items = [item.strip() for item in text.split(",") if item.strip()]
    names = []
    for item in items:
        found = _ATOM.match(item)
        if found is None:
            return None
        element, index, score = found.groups()
        names.append(AtomName(element, int(index), None if score is None else float(score)))
    return names


def names_atoms(texts: Iterable[str]) -> bool:
    """Whether *texts* are atom names: at least one value names atoms and none fails to read.

    What makes a column offerable as the source of a highlight. Only :data:`JUDGED` values are
    read, which is plenty to tell a column of sites from a column of prose.
    """
    found = False
    for position, text in enumerate(texts):
        if position >= JUDGED:
            break
        names = read_atom_names(text)
        if names is None:
            return False
        found = found or bool(names)
    return found


def atoms_named(mol: Chem.Mol, text: str) -> tuple[tuple[int, float | None], ...]:
    """The atoms of *mol* that *text* names, in its order, each with its score.

    An atom named twice is kept once, at its first place. A name whose index is outside *mol*
    or whose element is not the element of that atom is left out, so a column written for other
    structures draws nothing.
    """
    names = read_atom_names(text)
    if not names:
        return ()
    atoms = mol.GetNumAtoms()
    found: dict[int, float | None] = {}
    for name in names:
        if name.atom >= atoms or name.atom in found:
            continue
        if mol.GetAtomWithIdx(name.atom).GetSymbol() != name.element:
            continue
        found[name.atom] = name.score
    return tuple(found.items())


__all__ = ["JUDGED", "AtomName", "atoms_named", "names_atoms", "read_atom_names"]
