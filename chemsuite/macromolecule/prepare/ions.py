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

"""The formal charge of a monatomic ion.

Files often give an ion no charge (a PDB file's charge column is usually blank), and a force
field needs one. A charge the file states is kept; otherwise it comes from:

- the component's name, where the wwPDB Chemical Component Dictionary has one component per
  oxidation state (``FE`` is Fe³⁺ and ``FE2`` Fe²⁺, ``CU`` Cu²⁺ and ``CU1`` Cu⁺);
- else the element, where it has one common ion (Na⁺, Mg²⁺, Zn²⁺, Cl⁻).

An ion with neither (:func:`charge` returns None) is reported, never guessed.
"""

from __future__ import annotations

#: components whose name fixes the oxidation state (the CCD's formal charges)
BY_NAME = {
    "FE": 3, "FE2": 2, "CU": 2, "CU1": 1, "MN": 2, "MN3": 3, "CO": 2, "3CO": 3, "NI": 2,
    "3NI": 3, "CR": 3, "AU": 1, "AU3": 3, "TL": 1, "HG": 2, "PB": 2, "PT": 2, "PD": 2,
}  # fmt: skip
#: elements with one common ion, by atomic number
BY_ELEMENT = {
    3: 1, 11: 1, 19: 1, 37: 1, 55: 1,  # Li Na K Rb Cs
    4: 2, 12: 2, 20: 2, 38: 2, 56: 2,  # Be Mg Ca Sr Ba
    30: 2, 48: 2,  # Zn Cd
    13: 3, 31: 3,  # Al Ga
    47: 1,  # Ag
    9: -1, 17: -1, 35: -1, 53: -1,  # F Cl Br I
}  # fmt: skip


def charge(name: str, element: int, stated: int = 0) -> int | None:
    """The charge of the ion in component *name* (atomic number *element*); *stated* is the
    file's charge for it, kept when it is not zero. None when it cannot be told."""
    if stated:
        return int(stated)
    found = BY_NAME.get(name.strip().upper())
    if found is not None:
        return found
    return BY_ELEMENT.get(int(element))
