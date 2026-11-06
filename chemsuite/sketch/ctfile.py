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

#
# The valence rules below are a line-by-line port of ``Atom.calcValence`` (default mode) in
# Ketcher, ketcher-core/src/domain/entities/atom.ts at commit 2b031cc2, and its group numbers
# from domain/constants/elements.ts. Copyright 2021 EPAM Systems, Licensed under the Apache
# License, Version 2.0 (chemsuite/sketch/data/LICENSE-Ketcher.txt, NOTICE-Ketcher.txt).
# Changed from the original: translated to Python, restricted to the elements below.

"""The CTfile default valence model, as other programs read a molfile.

The sketcher's own hydrogens come from RDKit (``valence.py``). This model is used only when
writing a molfile: an atom whose hydrogen count a CTfile reader would work out differently
gets an explicit valence, so ChemDraw, Marvin or Ketcher read the same structure
(docs/sketcher/valence-review.md). Elements outside the table always get one.
"""

from __future__ import annotations

#: Ketcher's group numbers (domain/constants/elements.ts at the same commit)
ELEMENTS = {
    "H": 1, "Li": 1, "Na": 1, "K": 1, "Cu": 1, "Ag": 1,
    "Be": 2, "Mg": 2, "Ca": 2, "Zn": 2, "Hg": 2,
    "B": 3, "Al": 3, "Ga": 3, "Tl": 3,
    "C": 4, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4,
    "N": 5, "P": 5, "As": 5, "Sb": 5, "Bi": 5,
    "O": 6, "S": 6, "Se": 6, "Te": 6,
    "F": 7, "Cl": 7, "Br": 7, "I": 7,
    "Fe": 8, "Pt": 8,
}  # fmt: skip
INVALID = None  # the model reports a valence error


def default_hydrogens(label: str, charge: int, conn: int, rad: int) -> int | None:
    """Implicit H count by the CTfile default model; INVALID (None) for a valence error."""
    group = ELEMENTS[label]
    a = abs(charge)

    def result(value: int) -> int | None:
        return INVALID if value < 0 else value

    if group == 1:
        if label in ("H", "Li", "Na", "K", "Rb", "Cs", "Fr"):
            return result(1 - rad - conn - a)
        return 0
    if group == 2:
        return 0 if conn + rad + a in (0, 2) else INVALID
    if group == 3:
        if label in ("B", "Al", "Ga", "In"):
            return result(4 - rad - conn) if charge == -1 else result(3 - rad - conn - a)
        if label == "Tl":
            if charge == -1:
                return result(2 - rad - conn) if rad + conn <= 2 else result(4 - rad - conn)
            if charge == -2:
                return result(3 - rad - conn) if rad + conn <= 3 else result(5 - rad - conn)
            if rad + conn + a <= 1:
                return result(1 - rad - conn - a)
            return result(3 - rad - conn - a)
        return 0
    if group == 4:
        if label in ("C", "Si", "Ge"):
            return result(4 - rad - conn - a)
        if label in ("Sn", "Pb"):
            return result((2 if conn + rad + a <= 2 else 4) - rad - conn - a)
        return 0
    if group == 5:
        if label in ("N", "P"):
            if charge == 1:
                return result(4 - rad - conn)
            if charge == 2:
                return result(3 - rad - conn)
            return result((3 if rad + conn + a <= 3 else 5) - rad - conn - a)
        if label in ("Bi", "Sb", "As"):
            if charge == 1:
                if label != "As" and rad + conn <= 2:
                    return result(2 - rad - conn)
                return result(4 - rad - conn)
            if charge == 2:
                return result(3 - rad - conn)
            return result((3 if rad + conn <= 3 else 5) - rad - conn - a)
        return 0
    if group == 6:
        if label == "O":
            return result(3 - rad - conn) if charge >= 1 else result(2 - rad - conn - a)
        if label in ("S", "Se", "Po"):
            if charge == 1:
                return result((3 if conn <= 3 else 5) - rad - conn)
            total = conn + rad + a
            return result((2 if total <= 2 else 4 if total <= 4 else 6) - rad - conn - a)
        if label == "Te":
            if charge in (-1, 0, 2) and conn <= 2:
                return result(2 - rad - conn - a)
            if charge in (0, 2):
                if conn <= 4:
                    return result(4 - rad - conn - a)
                if charge == 0 and conn <= 6:
                    return result(6 - rad - conn - a)
                return INVALID
            return 0
        return 0
    if group == 7:
        if label == "F":
            return result(1 - rad - conn - a)
        if label in ("Cl", "Br", "I", "At"):
            if charge == 1:
                if conn <= 2:
                    return result(2 - rad - conn)
                if conn in (3, 5) or conn >= 7:
                    return INVALID
            elif charge == 0:
                if conn <= 1:
                    return result(1 - rad - conn)
                if conn in (2, 4, 6):
                    return 0 if rad == 1 else INVALID
                if conn > 7:
                    return INVALID
        return 0
    if group == 8:
        if label == "Pt":
            total = conn + rad + a
            if total <= 2:
                return result(2 - rad - conn - a)
            if total <= 4:
                return result(4 - rad - conn - a)
            return INVALID
        return 0 if conn + rad + a == 0 else INVALID
    raise AssertionError(label)
