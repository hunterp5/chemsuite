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
"""Hansch analysis: activity from substituent constants on a shared core.

    activity = a·π + b·π² + ρ·σ + c·(steric) + d

Constants come from a CSV the user keeps (published tables are not redistributed here):

    smiles,pi,sigma_m,sigma_p,f,r,es,mr,upsilon,sigma_star

``smiles`` is ``H`` or an attachment SMILES such as ``*C`` or ``*Cl``. A blank cell means the
constant is unknown; rows needing it are left out. Terms either sum over the positions
(classical Hansch: Σπ, Σσ) or get one coefficient per position.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .regression import Fit, least_squares
from .rgroups import HYDROGEN, Decomposition, substituent_label

FIELDS = ("pi", "sigma_m", "sigma_p", "f", "r", "es", "mr", "upsilon", "sigma_star")
HEADER = ",".join(("smiles", *FIELDS))
STERIC = {"": "None", "es": "Taft Es", "mr": "Molar refractivity (MR)", "upsilon": "Charton υ"}
ELECTRONIC = {"sigma": "Hammett σ", "field_resonance": "Field and resonance (F, R)", "": "None"}
SYMBOLS = {
    "pi": "π",
    "pi2": "π²",
    "sigma": "σ",
    "f": "F",
    "r": "R",
    "es": "Es",
    "mr": "MR",
    "upsilon": "υ",
}
#: Scales defined as zero at hydrogen (MR, Es and υ are not: hydrogen needs its own row).
ZERO_AT_HYDROGEN = frozenset({"pi", "pi2", "sigma", "f", "r"})
SKIP_UNMATCHED = "no core match"
SKIP_ACTIVITY = "no activity"
SKIP_CONSTANT = "substituent missing from the constants"


def table_key(text: str) -> str:
    """The label a constants-file entry matches: ``*C``, ``[*]C`` or ``C[*:1]`` → ``*C``."""
    text = text.strip()
    if text.upper() == HYDROGEN or text in ("[H]", "*[H]", "[*][H]"):
        return HYDROGEN
    return substituent_label(text)


def load_constants(path: str | Path) -> dict[str, dict[str, float]]:
    """``{label: {field: value}}``; blank cells are left out. Raises ValueError with the line."""
    table: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(line for line in handle if not line.lstrip().startswith("#"))
        names = {(n or "").strip().lower() for n in reader.fieldnames or ()}
        if "smiles" not in names:
            raise ValueError(f"The constants file needs a header like: {HEADER}")
        for number, row in enumerate(reader, start=2):
            row = {(k or "").strip().lower(): (v or "").strip() for k, v in row.items()}
            if not row.get("smiles"):
                continue
            values = {}
            for name in FIELDS:
                if row.get(name):
                    try:
                        values[name] = float(row[name])
                    except ValueError:
                        raise ValueError(
                            f"Line {number}: {name} is not a number: {row[name]!r}"
                        ) from None
            table[table_key(row["smiles"])] = values
    if not table:
        raise ValueError("The constants file has no substituents.")
    return table


@dataclass(frozen=True)
class HanschSettings:
    pi: bool = True
    pi_squared: bool = False
    electronic: str = "sigma"  # "sigma", "field_resonance" or ""
    sigma: dict[str, str] = field(default_factory=dict)  # position → "meta" / "para"
    default_sigma: str = "para"
    steric: str = ""  # "", "es", "mr", "upsilon"
    per_position: bool = False  # one coefficient per position instead of sums

    def terms(self) -> list[str]:
        out = ["pi"] if self.pi else []
        if self.pi and self.pi_squared:
            out.append("pi2")
        if self.electronic == "sigma":
            out.append("sigma")
        elif self.electronic == "field_resonance":
            out += ["f", "r"]
        if self.steric:
            out.append(self.steric)
        return out


@dataclass
class Hansch:
    fit: Fit
    terms: list[str]
    positions: tuple[str, ...]
    used: np.ndarray
    predicted: np.ndarray  # per input row (NaN where it cannot predict)
    skipped: Counter[str]
    missing: Counter[str]  # substituent labels with no constant needed by the model
    activity_name: str = "activity"

    def equation(self) -> str:
        parts = []
        for name, value in zip(self.fit.names, self.fit.coefficients, strict=True):
            if math.isfinite(value):
                sign = "−" if value < 0 else "+"
                parts.append(f" {sign} {abs(value):.3g}·{name}")
        text = f"{self.activity_name} = {self.fit.intercept:.3g}" + "".join(parts)
        return text


def _value(constants: dict[str, dict[str, float]], label: str, term: str, sigma: str) -> float:
    entry = constants.get(label)
    if entry is None:
        return math.nan
    if term == "pi2":
        v = entry.get("pi", math.nan)
        return v * v
    if term == "sigma":
        return entry.get("sigma_m" if sigma == "meta" else "sigma_p", math.nan)
    return entry.get(term, math.nan)


def hansch(
    decomposition: Decomposition,
    activity: np.ndarray,
    constants: dict[str, dict[str, float]],
    settings: HanschSettings,
    *,
    positions: Sequence[str] | None = None,
    activity_name: str = "activity",
) -> Hansch:
    terms = settings.terms()
    if not terms:
        raise ValueError("Choose at least one term (π, σ or a steric constant).")
    positions = tuple(positions or decomposition.varying())
    if not positions:
        raise ValueError("No position varies: Hansch analysis needs different substituents.")
    activity = np.asarray(activity, dtype=np.float64)
    n = len(decomposition)
    names = (
        [f"{SYMBOLS[t]}({p})" for p in positions for t in terms]
        if settings.per_position
        else [f"Σ{SYMBOLS[t]}" for t in terms]
    )
    x = np.full((n, len(names)), math.nan)
    missing: Counter[str] = Counter()
    for i in range(n):
        if not decomposition.matched[i]:
            continue
        per_term = {t: [] for t in terms}
        for p in positions:
            label = decomposition.labels[p][i]
            sigma = settings.sigma.get(p, settings.default_sigma)
            for t in terms:
                value = _value(constants, label, t, sigma)
                if label == HYDROGEN and not math.isfinite(value) and t in ZERO_AT_HYDROGEN:
                    value = 0.0  # hydrogen is the zero of π, σ, F and R by definition
                if not math.isfinite(value):
                    missing[label] += 1
                per_term[t].append(value)
        if settings.per_position:
            x[i] = [per_term[t][k] for k in range(len(positions)) for t in terms]
        else:
            x[i] = [sum(per_term[t]) for t in terms]
    skipped: Counter[str] = Counter()
    used = []
    for i in range(n):
        if not decomposition.matched[i]:
            skipped[SKIP_UNMATCHED] += 1
        elif not math.isfinite(activity[i]):
            skipped[SKIP_ACTIVITY] += 1
        elif not np.isfinite(x[i]).all():
            skipped[SKIP_CONSTANT] += 1
        else:
            used.append(i)
    if len(used) < len(names) + 2:
        raise ValueError(
            f"Only {len(used)} compound(s) have every constant needed for {len(names)} term(s); "
            "add constants or choose fewer terms."
        )
    fit = least_squares(x[used], activity[used], names)
    predicted = np.full(n, math.nan)
    ok = np.isfinite(x).all(axis=1)
    coefficients = np.nan_to_num(fit.coefficients, nan=0.0)
    predicted[ok] = fit.intercept + x[ok] @ coefficients
    return Hansch(fit, terms, positions, np.array(used), predicted, skipped, missing, activity_name)
