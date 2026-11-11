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
"""Free-Wilson / Fujita-Ban: activity as a sum of one contribution per substituent.

One least-squares fit, reported two ways. Fujita-Ban gives each substituent's effect against a
reference at its position (hydrogen when present there, else the most common group).
Free-Wilson recentres each position so the occurrence-weighted contributions sum to zero, plus
an overall mean. Predictions are identical. Unmade combinations of seen substituents are
ranked best-first, so even a vast combination space is never enumerated in full.
"""

from __future__ import annotations

import heapq
import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .regression import Fit, least_squares
from .rgroups import HYDROGEN, Decomposition

SKIP_UNMATCHED = "no core match"
SKIP_ACTIVITY = "no activity"
SKIP_RARE = "rare substituent"


@dataclass(frozen=True)
class Term:
    position: str
    label: str
    count: int
    reference: bool
    fujita_ban: float  # effect vs the reference at this position (0 for the reference)
    error: float
    free_wilson: float  # recentred contribution

    @property
    def estimable(self) -> bool:
        return math.isfinite(self.fujita_ban)


@dataclass
class FreeWilson:
    positions: tuple[str, ...]
    fit: Fit
    terms: list[Term]
    reference: dict[str, str]
    mean: float  # the Free-Wilson overall mean (μ)
    used: np.ndarray  # input row indices in the fit
    predicted: np.ndarray  # per input row: the model's value (NaN where it cannot predict)
    skipped: Counter[str] = field(default_factory=Counter)
    observed: set[tuple[str, ...]] = field(default_factory=set)

    def term(self, position: str, label: str) -> Term | None:
        return next((t for t in self.terms if t.position == position and t.label == label), None)

    def predict(self, choice: Sequence[str]) -> float:
        """The model's activity for one label per position (NaN if any is unknown)."""
        total = self.fit.intercept
        for position, label in zip(self.positions, choice, strict=True):
            term = self.term(position, label)
            if term is None or not term.estimable:
                return math.nan
            total += term.fujita_ban
        return total

    def best_unmade(
        self, limit: int, *, higher_is_better: bool = True
    ) -> list[tuple[tuple[str, ...], float]]:
        """The *limit* best predicted combinations not in the data, best first."""
        sign = 1.0 if higher_is_better else -1.0
        options = []
        for position in self.positions:
            terms = [t for t in self.terms if t.position == position and t.estimable]
            terms.sort(key=lambda t: -sign * t.fujita_ban)
            options.append(terms)
        if not all(options):
            return []

        def score(index: tuple[int, ...]) -> float:
            return sum(options[p][i].fujita_ban for p, i in enumerate(index))

        start = (0,) * len(options)
        heap = [(-sign * score(start), start)]
        seen = {start}
        out: list[tuple[tuple[str, ...], float]] = []
        while heap and len(out) < limit:
            _key, index = heapq.heappop(heap)
            labels = tuple(options[p][i].label for p, i in enumerate(index))
            if labels not in self.observed:
                out.append((labels, self.fit.intercept + score(index)))
            for p in range(len(index)):
                if index[p] + 1 < len(options[p]):
                    nxt = (*index[:p], index[p] + 1, *index[p + 1 :])
                    if nxt not in seen:
                        seen.add(nxt)
                        heapq.heappush(heap, (-sign * score(nxt), nxt))
        return out

    def combinations(self) -> int:
        sizes = [
            sum(1 for t in self.terms if t.position == p and t.estimable) for p in self.positions
        ]
        return math.prod(sizes) if sizes else 0


def _reference(counts: Counter[str]) -> str:
    if HYDROGEN in counts:
        return HYDROGEN
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def free_wilson(  # noqa: C901 - the published method, step by step
    decomposition: Decomposition,
    activity: np.ndarray,
    *,
    positions: Sequence[str] | None = None,
    min_count: int = 2,
) -> FreeWilson:
    """Fit the additive model. *positions* default to those whose substituent varies.

    Rows with a substituent seen fewer than *min_count* times are left out of the fit (one
    occurrence fits exactly and says nothing); they are counted in ``skipped``.
    """
    activity = np.asarray(activity, dtype=np.float64)
    positions = tuple(positions or decomposition.varying())
    if not positions:
        raise ValueError("No position varies: Free-Wilson needs different substituents.")
    n = len(decomposition)
    skipped: Counter[str] = Counter()
    candidates = []
    for i in range(n):
        if not decomposition.matched[i]:
            skipped[SKIP_UNMATCHED] += 1
        elif not math.isfinite(activity[i]):
            skipped[SKIP_ACTIVITY] += 1
        else:
            candidates.append(i)
    rows = list(candidates)
    while True:  # dropping a rare group can make another rare: repeat until stable
        counts = {p: Counter(decomposition.labels[p][i] for i in rows) for p in positions}
        rare = [
            i
            for i in rows
            if any(counts[p][decomposition.labels[p][i]] < min_count for p in positions)
        ]
        if not rare:
            break
        skipped[SKIP_RARE] += len(rare)
        dropped = set(rare)
        rows = [i for i in rows if i not in dropped]
    if len(rows) < 3:
        raise ValueError(
            f"Only {len(rows)} compound(s) remain after leaving out rare substituents; "
            "lower the minimum count or add data."
        )
    reference = {p: _reference(counts[p]) for p in positions}
    columns: list[tuple[str, str]] = []
    for p in positions:
        for label, _count in sorted(counts[p].items(), key=lambda kv: (-kv[1], kv[0])):
            if label != reference[p]:
                columns.append((p, label))
    x = np.zeros((len(rows), len(columns)))
    for j, (p, label) in enumerate(columns):
        x[:, j] = [decomposition.labels[p][i] == label for i in rows]
    fit = least_squares(x, activity[rows], [f"{p}={label}" for p, label in columns])

    terms: list[Term] = []
    mean = fit.intercept
    for p in positions:
        effects = {reference[p]: (0.0, math.nan)}
        for j, (q, label) in enumerate(columns):
            if q == p:
                effects[label] = (float(fit.coefficients[j]), float(fit.errors[j]))
        total = sum(counts[p].values())
        shift = sum(counts[p][lab] * eff for lab, (eff, _e) in effects.items()) / total
        mean += shift
        for label, (effect, error) in sorted(
            effects.items(), key=lambda kv: (-counts[p][kv[0]], kv[0])
        ):
            terms.append(
                Term(
                    p, label, counts[p][label], label == reference[p], effect, error, effect - shift
                )
            )
    model = FreeWilson(positions, fit, terms, reference, mean, np.array(rows), np.full(n, np.nan))
    model.skipped = skipped
    model.observed = {
        tuple(decomposition.labels[p][i] for p in positions)
        for i in range(n)
        if decomposition.matched[i]
    }
    for i in range(n):  # predict every matched row whose groups the model knows
        if decomposition.matched[i]:
            model.predicted[i] = model.predict([decomposition.labels[p][i] for p in positions])
    return model
