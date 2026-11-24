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
"""Everything a fit needs to be repeated: the feature pool, the split, the budget, the seed.

A :class:`Recipe` is plain data with a JSON round-trip. A fit stores its recipe in the model
archive; applying and refining a model replay it, so nothing about a model is implicit and the
same recipe, data and seed give the same model back.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

from chemsuite.analysis import fingerprints
from chemsuite.descriptors.catalog import CATALOG, PKA_KEYS, SHAPE_GROUP

#: What a model predicts. The task is detected from the target column, not asked for.
TASKS: dict[str, str] = {"regression": "Regression", "classification": "Classification"}

#: How the rows are divided. Scaffold is the default: a random split on a congeneric series
#: flatters every model, because the test rows' near neighbors sit in the training set.
SPLITS: dict[str, str] = {
    "scaffold": "By scaffold",
    "cluster": "By cluster",
    "random": "Random",
    "column": "By a column",
}

#: The fewest rows with a value that a fit will accept. Below this nothing can be estimated
#: honestly: five outer folds would hold out three rows each.
MIN_ROWS = 20
#: At most this many distinct values for a text or integer target to be classes, not numbers.
MAX_CLASSES = 10


@dataclass(frozen=True)
class Budget:
    """How much search a fit may spend.

    *configurations* is per candidate, before halving; *outer* and *inner* are fold counts.
    *candidates* limits the zoo to its cheapest entries, or is None for all of them.
    """

    key: str
    title: str
    configurations: int
    outer: int
    inner: int
    candidates: int | None = None

    @property
    def fits(self) -> int:
        """The order of magnitude of estimator fits, for the dialog's tooltip."""
        return self.configurations * self.inner * (self.outer + 1)


BUDGETS: dict[str, Budget] = {
    b.key: b
    for b in (
        Budget("quick", "Quick", configurations=8, outer=3, inner=3, candidates=4),
        Budget("standard", "Standard", configurations=24, outer=5, inner=5),
        Budget("thorough", "Thorough", configurations=80, outer=5, inner=5),
    )
}
DEFAULT_BUDGET = "standard"


def default_descriptors() -> tuple[str, ...]:
    """Every numeric descriptor that needs no conformer and no model of its own.

    The 3D shape group needs conformers the rows may not have, and the pKa descriptors cost a
    model run per row, so both are offered but neither is on by default.
    """
    return tuple(
        key
        for key, d in CATALOG.items()
        if not d.is_text and d.group != SHAPE_GROUP and key not in PKA_KEYS
    )


def default_fingerprints() -> tuple[str, ...]:
    """Morgan counts: how often each feature occurs, not just whether it is there."""
    return (fingerprints.variant("morgan2", counts=True),)


@dataclass(frozen=True)
class Recipe:
    """The feature pool and the fitting choices, as plain data."""

    task: str
    target: str
    structures: str = ""
    #: catalog keys (:data:`chemsuite.descriptors.catalog.CATALOG`)
    descriptors: tuple[str, ...] = ()
    #: fingerprint variant keys (:func:`chemsuite.analysis.fingerprints.variant`)
    fingerprints: tuple[str, ...] = ()
    #: numeric columns the table already has
    columns: tuple[str, ...] = ()
    split: str = "scaffold"
    #: the column the "column" split orders by
    split_column: str = ""
    #: the similarity at which the "cluster" split's Butina clusters are cut
    cluster_similarity: float = 0.6
    budget: str = DEFAULT_BUDGET
    #: the confidence a prediction interval or class set is built for
    confidence: float = 0.9
    seed: int = 0
    #: the selection cascade: drop a column correlated above this with one already kept, and
    #: stop at this many features
    correlation: float = 0.95
    max_features: int = 200
    #: candidates to leave out of the zoo, by key
    without: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.task not in TASKS:
            raise ValueError(f"unknown task {self.task!r}")
        if not self.target:
            raise ValueError("a recipe needs a target column")
        if self.split not in SPLITS:
            raise ValueError(f"unknown split {self.split!r}")
        if self.split == "column" and not self.split_column:
            raise ValueError("the column split needs a column to order by")
        if self.budget not in BUDGETS:
            raise ValueError(f"unknown budget {self.budget!r}")
        if not 0.5 <= self.confidence < 1.0:
            raise ValueError(f"a confidence is 0.5 up to 1, not {self.confidence}")
        if not 0.0 < self.correlation <= 1.0:
            raise ValueError(f"a correlation threshold is above 0 up to 1, not {self.correlation}")
        if self.max_features < 1:
            raise ValueError("a fit keeps at least one feature")
        unknown = [key for key in self.descriptors if key not in CATALOG]
        if unknown:
            raise ValueError(f"unknown descriptor(s): {', '.join(sorted(unknown))}")
        for key in self.fingerprints:
            fingerprints.parse(key)  # raises ValueError on anything it cannot read
        if (self.descriptors or self.fingerprints) and not self.structures:
            raise ValueError("descriptors and fingerprints need a structure column")
        if not (self.descriptors or self.fingerprints or self.columns):
            raise ValueError("a recipe needs at least one kind of feature")

    @property
    def budget_spec(self) -> Budget:
        return BUDGETS[self.budget]

    @property
    def needs_structures(self) -> bool:
        return bool(self.descriptors or self.fingerprints)

    def with_pool(
        self,
        descriptors: tuple[str, ...] | None = None,
        fingerprints: tuple[str, ...] | None = None,
        columns: tuple[str, ...] | None = None,
    ) -> Recipe:
        """The same recipe over another feature pool."""
        return replace(
            self,
            descriptors=self.descriptors if descriptors is None else descriptors,
            fingerprints=self.fingerprints if fingerprints is None else fingerprints,
            columns=self.columns if columns is None else columns,
        )

    def as_json(self) -> dict[str, Any]:
        """Plain types only, in a stable order, for the model archive's manifest."""
        return {
            "task": self.task,
            "target": self.target,
            "structures": self.structures,
            "descriptors": list(self.descriptors),
            "fingerprints": list(self.fingerprints),
            "columns": list(self.columns),
            "split": self.split,
            "split_column": self.split_column,
            "cluster_similarity": self.cluster_similarity,
            "budget": self.budget,
            "confidence": self.confidence,
            "seed": self.seed,
            "correlation": self.correlation,
            "max_features": self.max_features,
            "without": list(self.without),
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Recipe:
        """A recipe from :meth:`as_json`. Missing keys take their defaults; unknown keys are
        an error, so a recipe written by a newer chemsuite is refused rather than guessed at."""
        known = set(cls.__dataclass_fields__)
        unknown = sorted(set(data) - known)
        if unknown:
            raise ValueError(f"the recipe has setting(s) this chemsuite does not know: {unknown}")
        values = dict(data)
        for key in ("descriptors", "fingerprints", "columns", "without"):
            if key in values:
                values[key] = tuple(values[key])
        return cls(**values)


def detect_task(numeric: bool, distinct: int) -> str:
    """Whether a target column holds a regression or a classification target.

    Text is always classes. A number column is classes only when few values repeat over the
    whole column — a 0/1 flag, or a graded score — and a regression target otherwise. The
    dialog shows what this decided, so a graded score can be modelled either way.
    """
    if not numeric:
        return "classification"
    return "classification" if 2 <= distinct <= MAX_CLASSES else "regression"
