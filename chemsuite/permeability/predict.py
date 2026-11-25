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

"""Permeability for many structures: chunked, cancellable, in the job's own thread.

The shape is ``metabolism/predict.py``'s: chunks, progress after each one, and a
:class:`Cancelled` carrying the rows already finished. The work per row is a feature pool and
a numpy forward pass, so a thousand rows are a second or two and worker processes would cost
more than they save — except for the descriptor pool itself, which takes an *executor* when the
caller has one.

**Why it predicts per chunk rather than building one pool.** The pool is a few thousand columns
wide, so a hundred thousand rows at once would be gigabytes. Each chunk is featurized,
predicted and reduced to the handful of numbers a column needs; the per-endpoint results are
stacked at the end into one :class:`~chemsuite.qsar.model.Prediction`, which is what
``tools/qsar/columns.py`` already knows how to write.

A row whose structure RDKit cannot read gets no prediction rather than a filled-in guess: the
pool marks it unusable and :meth:`Endpoint Results.blank` puts NaN and an empty class set
there. A blank cell is not counted as a failure; an unreadable one is.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

import numpy as np

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.qsar import features
from chemsuite.qsar.model import Model, Prediction

from . import ENDPOINTS, Endpoint
from .models import model as shipped_model

#: A structure cell: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
#: Rows per chunk: enough that progress and a cancel stay prompt without being chatty, and
#: small enough that one chunk's feature pool stays a few megabytes.
CHUNK_ROWS = 200

Progress = Callable[[int, int], None]


def is_blank(cell: Cell) -> bool:
    """A cell with nothing in it: no structure, and so not a failure to count."""
    return cell is None or (isinstance(cell, (bytes, str)) and not cell.strip())


@dataclass
class EndpointResults:
    """One endpoint's predictions, growing a chunk at a time."""

    endpoint: Endpoint
    predicted: list[float] = field(default_factory=list)
    low: list[float] = field(default_factory=list)
    high: list[float] = field(default_factory=list)
    distance: list[float] = field(default_factory=list)
    inside: list[bool] = field(default_factory=list)
    #: classification only: one score row per row, and the class set as a boolean row
    scores: list[np.ndarray] = field(default_factory=list)
    sets: list[np.ndarray] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.predicted)

    @property
    def classifying(self) -> bool:
        return self.endpoint.task == "classification"

    def blank(self, classes: int) -> None:
        """A row the model has nothing to say about."""
        self.predicted.append(np.nan)
        self.low.append(np.nan)
        self.high.append(np.nan)
        self.distance.append(np.nan)
        self.inside.append(False)
        if self.classifying:
            self.scores.append(np.full(classes, np.nan))
            self.sets.append(np.zeros(classes, dtype=np.bool_))

    def add(self, said: Prediction, row: int) -> None:
        """Row *row* of a chunk's prediction."""
        self.predicted.append(float(said.predicted[row]))
        self.low.append(np.nan if said.low is None else float(said.low[row]))
        self.high.append(np.nan if said.high is None else float(said.high[row]))
        self.distance.append(np.nan if said.distance is None else float(said.distance[row]))
        self.inside.append(bool(said.inside[row]) if said.inside is not None else False)
        if self.classifying:
            self.scores.append(
                np.array([]) if said.scores is None else np.asarray(said.scores[row])
            )
            self.sets.append(
                np.array([], dtype=np.bool_)
                if said.sets is None
                else np.asarray(said.sets[row], dtype=np.bool_)
            )

    def extend(self, other: EndpointResults) -> None:
        self.predicted.extend(other.predicted)
        self.low.extend(other.low)
        self.high.extend(other.high)
        self.distance.extend(other.distance)
        self.inside.extend(other.inside)
        self.scores.extend(other.scores)
        self.sets.extend(other.sets)

    def as_prediction(self) -> Prediction:
        """The rows so far as one :class:`Prediction`, which writes columns like any other."""
        numbers = np.asarray(self.predicted, dtype=np.float64)
        if self.classifying:
            return Prediction(
                predicted=np.nan_to_num(numbers, nan=-1.0).astype(np.int64),
                scores=np.vstack(self.scores) if self.scores else None,
                sets=np.vstack(self.sets) if self.sets else None,
                distance=np.asarray(self.distance, dtype=np.float64),
                inside=np.asarray(self.inside, dtype=np.bool_),
            )
        return Prediction(
            predicted=numbers,
            low=np.asarray(self.low, dtype=np.float64),
            high=np.asarray(self.high, dtype=np.float64),
            distance=np.asarray(self.distance, dtype=np.float64),
            inside=np.asarray(self.inside, dtype=np.bool_),
        )


@dataclass
class Results:
    """What a run produced: one :class:`EndpointResults` per endpoint, in the order asked for."""

    per_endpoint: dict[str, EndpointResults] = field(default_factory=dict)
    failed: int = 0

    @classmethod
    def empty(cls, endpoints: Sequence[Endpoint]) -> Results:
        return cls({e.id: EndpointResults(e) for e in endpoints})

    def __len__(self) -> int:
        return len(next(iter(self.per_endpoint.values()))) if self.per_endpoint else 0

    def extend(self, other: Results) -> None:
        for endpoint_id, results in other.per_endpoint.items():
            self.per_endpoint[endpoint_id].extend(results)
        self.failed += other.failed


def _classes(model: Model) -> int:
    return 0 if model.classes is None else int(model.classes.size)


def predict_chunk(
    cells: Sequence[Cell],
    endpoints: Sequence[Endpoint],
    *,
    executor: Executor | None = None,
) -> Results:
    """One chunk of cells, for each endpoint."""
    out = Results.empty(endpoints)
    counted = False
    for endpoint in endpoints:
        model = shipped_model(endpoint)
        results = out.per_endpoint[endpoint.id]
        pool = features.build(cells, model.recipe, executor=executor)
        said = model.predict(pool.x)
        classes = _classes(model)
        for row, cell in enumerate(cells):
            if pool.usable[row]:
                results.add(said, row)
                continue
            results.blank(classes)
            # one failure per row, however many endpoints could not read it
            if not counted and not is_blank(cell):
                out.failed += 1
        counted = True
    return out


def predict_permeability(
    cells: Sequence[Cell],
    endpoints: Sequence[Endpoint] = ENDPOINTS,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Results:
    """Every chosen endpoint's prediction for each cell, in order.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` holds the results for the
    rows finished so far (the first rows, in order), or None if none were.
    """
    endpoints = tuple(endpoints)
    results = Results.empty(endpoints)
    total = len(cells)
    for start in range(0, total, CHUNK_ROWS):
        if cancel is not None and cancel.cancelled:
            raise Cancelled(results if len(results) else None)
        results.extend(
            predict_chunk(cells[start : start + CHUNK_ROWS], endpoints, executor=executor)
        )
        if progress is not None:
            progress(len(results), total)
    return results
