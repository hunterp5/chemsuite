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
"""Concentrations in one unit (Tools ▸ Standardize Values): potencies given as a value and a unit
(nM, µM, mg/mL…) row by row, converted to the unit chosen.

Molar units convert among themselves by powers of ten, and mass per volume units likewise; a
molar value converts to a mass per volume one (or back) through the row's molecular weight:
``g/L = mol/L × g/mol``.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import Cancelled, CancelToken, check

#: molar units, in mol/L
MOLAR = {"M": 1.0, "mM": 1e-3, "µM": 1e-6, "nM": 1e-9, "pM": 1e-12, "fM": 1e-15}
#: mass per volume units, in g/L
MASS = {
    "g/L": 1.0,
    "mg/mL": 1.0,
    "mg/L": 1e-3,
    "µg/mL": 1e-3,
    "µg/L": 1e-6,
    "ng/mL": 1e-6,
    "ng/L": 1e-9,
    "pg/mL": 1e-9,
}
UNITS = (*MOLAR, *MASS)

_PREFIX = {"": "", "m": "m", "u": "µ", "µ": "µ", "μ": "µ", "mc": "µ", "n": "n", "p": "p",
           "f": "f"}  # fmt: skip
_WORDS = {"molar": "M", "millimolar": "mM", "micromolar": "µM", "nanomolar": "nM",
          "picomolar": "pM", "femtomolar": "fM"}  # fmt: skip
_MOLAR_RE = re.compile(r"^(mc|[munpfµμ]?)(m|mol/l|mol l-1|mol\.l-1)$")
_MASS_RE = re.compile(r"^(mc|[munpµμ]?)g\s*/\s*(m|d)?l$")
_NUMBER_RE = re.compile(r"^\s*(<=|>=|≤|≥|<|>|~|=)?\s*([-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?)\s*$")


def parse_unit(text: str | None) -> str | None:
    """*text* as one of :data:`UNITS` (None if it is not a concentration read here). Case does
    not matter (``uM``, ``UM`` and ``µM`` are all micromolar); ``mol/L`` and the words
    (``nanomolar``) are read too."""
    if text is None:
        return None
    key = str(text).strip().replace(" ", "").lower()
    if not key:
        return None
    if key in _WORDS:
        return _WORDS[key]
    found = _MOLAR_RE.match(key)
    if found:
        return _PREFIX[found.group(1)] + "M"
    found = _MASS_RE.match(key)
    if found:
        prefix, per = _PREFIX[found.group(1)], found.group(2) or ""
        if per == "d":  # per deciliter: not read (mg/dL is rare for potencies)
            return None
        unit = f"{prefix}g/{per}L" if per else f"{prefix}g/L"
        return unit if unit in MASS else None
    return None


def parse_value(value: object) -> tuple[str, float]:
    """``(qualifier, number)`` from a cell: a number, or text such as ``>10`` or ``≤ 0.5``
    (``nan`` if there is no number)."""
    if value is None:
        return "", np.nan
    if isinstance(value, (int, float, np.floating)):
        return "", float(value)
    found = _NUMBER_RE.match(str(value))
    if found is None:
        return "", np.nan
    qualifier = {"≤": "<=", "≥": ">=", "=": ""}.get(found.group(1) or "", found.group(1) or "")
    return qualifier, float(found.group(2))


def needs_weight(units: Sequence[str | None], target: str) -> bool:
    """Whether some rows cross between molar and mass units (so need a molecular weight)."""
    molar = target in MOLAR
    return any(u is not None and (u in MOLAR) != molar for u in units)


@dataclass
class Converted:
    """Values in the target unit (NaN where they could not be), the qualifiers read with them
    (``""`` for none), and how many rows had no unit read or no molecular weight to use."""

    values: np.ndarray
    qualifiers: list[str]
    unknown_units: int
    no_weight: int


def convert(
    values: Sequence[object],
    units: Sequence[str | None],
    target: str,
    weights: np.ndarray | None = None,
) -> Converted:
    """*values* (numbers, or text with qualifiers) in *units* (as typed: ``uM``…) converted to
    *target*. *weights* (g/mol, NaN where unknown) are used where a row crosses between molar
    and mass units."""
    if target not in UNITS:
        raise ValueError(f"unknown unit: {target!r}")
    out = np.full(len(values), np.nan)
    qualifiers = [""] * len(values)
    unknown = no_weight = 0
    to_molar = target in MOLAR
    for i, (value, typed) in enumerate(zip(values, units, strict=True)):
        qualifier, number = parse_value(value)
        if np.isnan(number):
            continue
        unit = parse_unit(typed)
        if unit is None:
            unknown += 1
            continue
        base = number * (MOLAR[unit] if unit in MOLAR else MASS[unit])  # mol/L or g/L
        if (unit in MOLAR) != to_molar:
            weight = np.nan if weights is None else float(weights[i])
            if not weight > 0:
                no_weight += 1
                continue
            base = base / weight if to_molar else base * weight
        out[i] = base / (MOLAR[target] if to_molar else MASS[target])
        qualifiers[i] = qualifier
    return Converted(out, qualifiers, unknown, no_weight)


def molecular_weights(cells: Sequence[bytes | str | None]) -> np.ndarray:
    """Average molecular weight (g/mol) of each structure cell, as drawn (salts included);
    NaN for a blank or unreadable one."""
    from rdkit.Chem import Descriptors

    from chemsuite.core.mols import mol_from_blob, mol_from_text

    out = np.full(len(cells), np.nan)
    for i, cell in enumerate(cells):
        mol = None
        if isinstance(cell, bytes) and cell:
            mol = mol_from_blob(cell)
        elif isinstance(cell, str) and cell:
            mol = mol_from_text(cell)
        if mol is not None:
            out[i] = Descriptors.MolWt(mol)
    return out


#: rows weighed and converted between checks for a cancel
CHUNK_ROWS = 500


def convert_rows(
    values: Sequence[object],
    units: Sequence[str | None],
    target: str,
    cells: Sequence[bytes | str | None] | None = None,
    *,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> Converted:
    """:func:`convert` with molecular weights from structure *cells* (None: no row needs
    one), a chunk at a time. If *cancel* is set, raises :class:`Cancelled` whose ``partial``
    holds the first rows done (None if none were)."""
    total = len(values)
    parts: list[Converted] = []

    def joined() -> Converted:
        return Converted(
            np.concatenate([p.values for p in parts]) if parts else np.empty(0),
            [q for p in parts for q in p.qualifiers],
            sum(p.unknown_units for p in parts),
            sum(p.no_weight for p in parts),
        )

    for start in range(0, total, CHUNK_ROWS):
        try:
            check(cancel)
        except Cancelled:
            raise Cancelled(joined() if parts else None) from None
        stop = min(total, start + CHUNK_ROWS)
        weights = molecular_weights(cells[start:stop]) if cells is not None else None
        parts.append(convert(values[start:stop], units[start:stop], target, weights))
        if progress is not None:
            progress(stop, total)
    return joined() if parts else convert([], [], target)
