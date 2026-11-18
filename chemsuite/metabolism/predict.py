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

"""Sites of metabolism for many structures: chunked, cancellable, in the job's own thread.

The pKa driver (``ionisation/predict.py``) farms its rows out to worker processes because
one molecule takes seconds; this model takes about a millisecond, so a thousand rows are
a second of work and the processes would cost more than they save. The shape is otherwise
the same: chunks, progress after each, and a :class:`Cancelled` carrying the rows already
finished.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

from rdkit import Chem, rdBase

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import mol_from_blob, mol_from_text

from . import ISOFORMS
from .sites import MetaboliteSite, predict_sites

#: A structure cell: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
#: Rows per chunk: enough that progress and a cancel stay prompt without being chatty.
CHUNK_ROWS = 25
#: How many ranked atoms the ``sites`` text lists.
RANKED = 3

_FAILURES = (ValueError, RuntimeError, ArithmeticError, KeyError, IndexError)


@dataclass
class SiteResults:
    """One entry per row in each list, in order; "" or NaN where a row has no value."""

    #: per isoform, the highest-scoring atom as its element and index, e.g. ``"C7"``
    top_site: dict[str, list[str]] = field(default_factory=dict)
    #: per isoform, that atom's score
    top_score: dict[str, list[float]] = field(default_factory=dict)
    #: per isoform, the :data:`RANKED` best atoms with their scores
    sites: dict[str, list[str]] = field(default_factory=dict)
    #: per isoform, the share of the top atom's environments seen in training
    domain: dict[str, list[float]] = field(default_factory=dict)
    #: per isoform, the ranked atom indices themselves, for highlighting the structure
    ranked_atoms: dict[str, list[list[int]]] = field(default_factory=dict)
    failed: int = 0

    @classmethod
    def empty(cls, isoforms: Sequence[str]) -> SiteResults:
        return cls(
            top_site={name: [] for name in isoforms},
            top_score={name: [] for name in isoforms},
            sites={name: [] for name in isoforms},
            domain={name: [] for name in isoforms},
            ranked_atoms={name: [] for name in isoforms},
        )

    def __len__(self) -> int:
        return len(next(iter(self.top_score.values()))) if self.top_score else 0

    def extend(self, other: SiteResults) -> None:
        for isoform, values in other.top_site.items():
            self.top_site[isoform].extend(values)
            self.top_score[isoform].extend(other.top_score[isoform])
            self.sites[isoform].extend(other.sites[isoform])
            self.domain[isoform].extend(other.domain[isoform])
            self.ranked_atoms[isoform].extend(other.ranked_atoms[isoform])
        self.failed += other.failed


def atom_name(mol: Chem.Mol, atom: int) -> str:
    """How a site is named in a column: its element and its index in the molecule."""
    return f"{mol.GetAtomWithIdx(atom).GetSymbol()}{atom}"


def sites_text(mol: Chem.Mol, ranked: Sequence[MetaboliteSite], count: int = RANKED) -> str:
    """The *count* best atoms, best first, each with its score."""
    return ", ".join(f"{atom_name(mol, site.atom)} {site.score:.2f}" for site in ranked[:count])


def is_blank(cell: Cell) -> bool:
    """A cell with nothing in it: no structure, and so not a failure to count."""
    return cell is None or (isinstance(cell, (bytes, str)) and not cell.strip())


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell.strip() else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def _blank(results: SiteResults, isoforms: Sequence[str]) -> None:
    for isoform in isoforms:
        results.top_site[isoform].append("")
        results.top_score[isoform].append(math.nan)
        results.sites[isoform].append("")
        results.domain[isoform].append(math.nan)
        results.ranked_atoms[isoform].append([])


def predict_chunk(cells: Sequence[Cell], isoforms: Sequence[str]) -> SiteResults:
    """Results for one chunk of cells."""
    out = SiteResults.empty(isoforms)
    with rdBase.BlockLogs():
        for cell in cells:
            mol = _mol(cell)
            ranked: dict[str, list[MetaboliteSite]] = {}
            if mol is not None:
                try:
                    ranked = predict_sites(mol, tuple(isoforms))
                except _FAILURES:
                    ranked = {}
            if not any(ranked.values()):
                out.failed += not is_blank(cell)
                _blank(out, isoforms)
                continue
            for isoform in isoforms:
                best = ranked[isoform][0]
                out.top_site[isoform].append(atom_name(mol, best.atom))
                out.top_score[isoform].append(best.score)
                out.sites[isoform].append(sites_text(mol, ranked[isoform]))
                out.domain[isoform].append(best.domain)
                out.ranked_atoms[isoform].append([site.atom for site in ranked[isoform][:RANKED]])
    return out


def predict_metabolism(
    cells: Sequence[Cell],
    isoforms: Sequence[str] = ISOFORMS,
    *,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> SiteResults:
    """The ranked sites of metabolism for each cell, for each of *isoforms*.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` holds the results for
    the rows finished so far (the first rows, in order), or None if none were.
    """
    isoforms = tuple(isoforms)
    results = SiteResults.empty(isoforms)
    total = len(cells)
    for start in range(0, total, CHUNK_ROWS):
        if cancel is not None and cancel.cancelled:
            raise Cancelled(results if len(results) else None)
        results.extend(predict_chunk(cells[start : start + CHUNK_ROWS], isoforms))
        if progress is not None:
            progress(len(results), total)
    return results
