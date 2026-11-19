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

"""How many charges a molecule's microstates span: start where the model was evaluated, widen
while the populated forms reach an end.

Upstream's enumerator loop stops one charge short of each end it is given: its default −2..+2
(:data:`START_RANGE`) gives charges −1..+1, which is what docs/UNIPKA_PLAN.md U2 evaluated.
:func:`widen` steps out one charge at a time, up to :data:`WIDEST_RANGE` (charges −4..+3:
citrate, ATP), while the caller says an end is populated. Each microstate costs the model about
0.75 s, so a molecule with more than :data:`MAX_MICROSTATES` raises
:class:`TooManyMicrostates`, and a wider range that would exceed it is not taken.

:func:`plausible` leaves out the microstates that protonate an oxygen or sulfur (a carbonyl or
an alcohol: oxonium ions). Upstream's templates make them, but their pKa values lie below
about −2, outside what the model was trained on, where it puts them at 1–6 (acetic acid's at
1.85). On the external tests, leaving them out changes nothing (RMSE 0.808 either way).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

from rdkit import Chem

from .microstates import Microstate

START_RANGE, WIDEST_RANGE = (-2, 2), (-5, 4)
MAX_MICROSTATES = 48

T = TypeVar("T")
Pools = dict[int, list[Microstate]]


class TooManyMicrostates(Exception):
    """A molecule has more microstates than the model is given (:data:`MAX_MICROSTATES`)."""

    def __init__(self, count: int) -> None:
        super().__init__(f"{count} protonation states")
        self.count = count


def plausible(pools: Pools, mol: Chem.Mol) -> Pools:
    """*pools* (microstates of *mol*) less those with an O or S more positive than in *mol*."""
    atoms = [(k, a.GetFormalCharge()) for k, a in enumerate(mol.GetAtoms())
             if a.GetAtomicNum() in (8, 16)]  # fmt: skip
    out = {}
    for q, pool in pools.items():
        kept = [m for m in pool if all(m.charges[k] <= q0 for k, q0 in atoms)]
        if kept:
            out[q] = kept
    return out


def widen(
    pools_of: Callable[[int, int], Pools],
    accept: Callable[[Pools], list[T]],
    ends: Callable[[list[T], int, int], tuple[bool, bool]],
) -> list[T]:
    """``accept(pools)`` for the widest charge range worth taking.

    *pools_of* enumerates for a range (as :func:`.microstates.ensemble` takes it); *accept*
    scores the pools; *ends*, given that result and the lowest and highest charge found, says
    whether to widen below and above. Raises :class:`TooManyMicrostates` when even
    :data:`START_RANGE` has too many.
    """
    (lower, upper), charges, found = START_RANGE, set(), []
    while True:
        pools = pools_of(lower, upper)
        count = sum(len(pool) for pool in pools.values())
        if count > MAX_MICROSTATES:
            if not charges:
                raise TooManyMicrostates(count)
            return found
        if set(pools) == charges:
            return found  # widening found no new charge
        charges, found = set(pools), accept(pools)
        if not found:
            return found
        below, above = ends(found, min(pools), max(pools))
        low = lower - 1 if below and lower > WIDEST_RANGE[0] else lower
        high = upper + 1 if above and upper < WIDEST_RANGE[1] else upper
        if (low, high) == (lower, upper):
            return found
        lower, upper = low, high
