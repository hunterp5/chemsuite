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
"""What a similarity search takes and gives: the references, the methods and the scores."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from chemsuite.analysis.fingerprints import (
    ASYMMETRIC_METRICS,
    DEFAULT_KIND,
    KINDS,
    METRICS,
    Fingerprints,
)

#: the methods
FINGERPRINT = "fingerprint"  # any fingerprint kind, the 2D pharmacophore ("pharm2d") among them
MCS = "mcs"
SCAFFOLD = "scaffold"
GENERIC_SCAFFOLD = "generic_scaffold"
SHAPE = "shape"
COLOR = "shape_color"
FIELD = "field"
METHODS = (FINGERPRINT, MCS, SCAFFOLD, GENERIC_SCAFFOLD, SHAPE, COLOR, FIELD)
#: methods too slow to run on every row of a big table or file (see ``Search.prefilter``)
EXPENSIVE = frozenset({MCS, SHAPE, COLOR, FIELD})
THREE_D = frozenset({SHAPE, COLOR, FIELD})

#: how one row's scores against several references become one
FUSIONS = {"max": "Best", "mean": "Mean"}
#: how several methods' scores become one
CONSENSUS = {"mean": "Mean", "rank": "Mean Rank", "normalized": "Normalized Mean"}

CONSENSUS_COLUMN = "Consensus"
NEAREST_COLUMN = "Nearest Reference"
_TITLES = {
    MCS: "MCS Tanimoto",
    SCAFFOLD: "Scaffold Match",
    GENERIC_SCAFFOLD: "Generic Scaffold Match",
    SHAPE: "Shape Tanimoto",
    COLOR: "Shape + Features",
    FIELD: "Field Similarity",
}


@dataclass(frozen=True)
class Reference:
    """A structure to compare with. *row*, when the reference is one of the molecules being
    scored, is its position among them, so it is not compared with itself (-1: none)."""

    mol: Chem.Mol
    name: str
    row: int = -1


@dataclass(frozen=True)
class Method:
    """One way of scoring. *fingerprint* and *metric* are for the fingerprint method; *alpha*
    weighs the bits only the reference has and *beta* those only the row has (Tversky)."""

    kind: str = FINGERPRINT
    fingerprint: str = DEFAULT_KIND
    metric: str = "tanimoto"
    alpha: float = 0.5
    beta: float = 0.5

    def __post_init__(self) -> None:
        if self.kind not in METHODS:
            raise ValueError(f"unknown method: {self.kind!r}")
        if self.kind == FINGERPRINT:
            if self.fingerprint not in KINDS:
                raise ValueError(f"unknown fingerprint: {self.fingerprint!r}")
            if self.metric not in METRICS and self.metric not in ASYMMETRIC_METRICS:
                raise ValueError(f"unknown similarity metric: {self.metric!r}")
        if not 0.0 <= self.alpha <= 1.0 or not 0.0 <= self.beta <= 1.0:
            raise ValueError("Tversky weights must be between 0 and 1")
        if self.alpha == self.beta == 0.0:
            raise ValueError("Tversky weights cannot both be 0")

    @property
    def column(self) -> str:
        """The name of its score column."""
        if self.kind != FINGERPRINT:
            return _TITLES[self.kind]
        metric = METRICS.get(self.metric) or ASYMMETRIC_METRICS[self.metric]
        return f"{KINDS[self.fingerprint].column} {metric}"

    @property
    def expensive(self) -> bool:
        return self.kind in EXPENSIVE


@dataclass(frozen=True)
class Search:
    """The methods (the first one is primary: it names each row's nearest reference), how to
    fuse the scores over references and combine them over methods, and *prefilter*: how many of
    the rows the primary fingerprint ranks best go on to the expensive methods (0: all)."""

    methods: tuple[Method, ...]
    fuse: str = "max"
    consensus: str = "mean"
    prefilter: int = 500

    def __post_init__(self) -> None:
        if not self.methods:
            raise ValueError("choose at least one method")
        if len({m.column for m in self.methods}) != len(self.methods):
            raise ValueError("a method is chosen twice")
        if self.fuse not in FUSIONS:
            raise ValueError(f"unknown fusion: {self.fuse!r}")
        if self.consensus not in CONSENSUS:
            raise ValueError(f"unknown consensus: {self.consensus!r}")
        if self.prefilter < 0:
            raise ValueError("the prefilter cannot be negative")

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(m.column for m in self.methods)


@dataclass
class Scores:
    """Results row for row: each method's score (NaN: not scored), the index of each row's
    nearest reference by the primary method (-1: none), the consensus when there are several
    methods, and the fingerprints computed on the way (for a cache)."""

    columns: dict[str, np.ndarray]
    nearest: np.ndarray
    consensus: np.ndarray | None = None
    fingerprints: dict[str, Fingerprints] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.nearest)

    def head(self, rows: int) -> Scores:
        """The first *rows* rows."""
        consensus = None if self.consensus is None else self.consensus[:rows]
        columns = {name: values[:rows] for name, values in self.columns.items()}
        return Scores(columns, self.nearest[:rows], consensus, self.fingerprints)
