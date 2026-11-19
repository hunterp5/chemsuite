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

"""The atoms a structure column's drawings mark, and what marks them.

Two kinds, one per column at a time (see :meth:`.TableView.set_highlight`), both display only:

* a :class:`Highlight` — the atoms matching a substructure query, the same query in every row,
  which Structure Search sets;
* an :class:`AtomSource` — the atoms another column's values *name* (``C7``, or
  ``C7 0.81, N3 0.05``; :mod:`chemsuite.core.atom_names` reads them), which differ row by row
  and are resolved while painting. The sites of metabolism are written that way, so Highlight
  Atoms ▸ in a structure column's header menu marks them on the structures.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem

from chemsuite.analysis.structure_search import match_atoms
from chemsuite.core.atom_names import JUDGED, atoms_named, names_atoms
from chemsuite.core.depict import depiction_mol, drawn_highlight
from chemsuite.core.table import TEXT, Table
from chemsuite.core.table.columns import Column, ParsedColumn

if TYPE_CHECKING:
    from .table_view import TableView


@dataclass(frozen=True)
class Highlight:
    """Atoms to highlight in a structure column: those matching *query* (a substructure)."""

    key: str  # identifies the query in the image cache
    query: Chem.Mol
    stereo: bool = False


@dataclass(frozen=True)
class AtomSource:
    """Atoms to highlight in a structure column: those the column *source* names, row by row."""

    source: str


@dataclass(frozen=True)
class AtomHighlight:
    """One row's share of an :class:`AtomSource`: the atoms it names, best first."""

    key: str  # identifies the row's atoms in the image cache
    atoms: tuple[tuple[int, float | None], ...]  # (atom index, score)


def resolve(table: Table, column: str, source: AtomSource, slot: int) -> AtomHighlight | None:
    """The atoms named for one row of structure *column*; None when nothing names any.

    Resolved while painting, so the drawings follow *source*: write that column again and the
    marks move with it.
    """
    if source.source not in table:
        return None
    text = table.column(source.source).display(slot)
    structures = table.column(column)
    if not text or not isinstance(structures, ParsedColumn):
        return None
    mol = structures.obj(slot)
    if mol is None:
        return None
    atoms = atoms_named(mol, text)
    return AtomHighlight(f"{source.source}|{text}", atoms) if atoms else None


def marked(
    mol: Chem.Mol, highlight: Highlight | AtomHighlight, dark: bool
) -> tuple[Chem.Mol, list[int], list[int], dict[int, tuple[float, float, float]]]:
    """The molecule as it is drawn, and the atoms, bonds and atom colors to draw marked."""
    shown = depiction_mol(mol)
    if isinstance(highlight, AtomHighlight):
        atoms, colors = drawn_highlight(shown, highlight.atoms, dark)
        return shown, atoms, [], colors
    atoms, bonds = match_atoms(shown, highlight.query, highlight.stereo)
    return shown, atoms, bonds, {}


# -- the columns that can mark atoms, and the view's current choice ---------------------------
def set_source(view: TableView, column: str, source: str | None) -> None:
    """Mark the atoms *source*'s values name in *column*'s drawings (None: stop)."""
    view.set_highlight(column, AtomSource(source) if source else None)


def sources(view: TableView) -> dict[str, str]:
    """Structure column -> the column whose values name the atoms its drawings mark."""
    return {
        name: highlight.source
        for name, highlight in view.highlights().items()
        if isinstance(highlight, AtomSource)
    }


def atom_name_columns(table: Table, structures: str) -> list[str]:
    """The text columns whose values name atoms of *structures*, in table order.

    A value must read as atom names *and* name atoms the molecule beside it really has, so a
    column of compound ids, or of indices written for other structures, is not offered.
    """
    slots = np.arange(min(len(table), JUDGED))
    drawn = table.column(structures)
    found = []
    for name in table.column_names:
        if table.kind(name) != TEXT:
            continue
        texts = table.column(name).displays(slots)
        if names_atoms(texts) and _names_these_atoms(drawn, texts, slots):
            found.append(name)
    return found


def _names_these_atoms(drawn: Column, texts: Sequence[str], slots: np.ndarray) -> bool:
    """Whether *texts* name atoms of the molecules in *drawn*: every row that has both a value
    and a molecule must name at least one atom of it, and some row must."""
    if not isinstance(drawn, ParsedColumn):
        return False
    found = False
    for text, slot in zip(texts, slots.tolist(), strict=True):
        mol = drawn.obj(slot) if text else None
        if mol is None:
            continue
        if not atoms_named(mol, text):
            return False
        found = True
    return found
