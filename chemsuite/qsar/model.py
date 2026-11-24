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
"""The saved model: a ``.csqm`` archive, and what it predicts.

A ``.csqm`` file is an uncompressed ``.npz`` — a format number, the manifest as JSON, and
plain numeric or fixed-width text arrays beside it — read with ``allow_pickle=False``, the
same shape as :mod:`chemsuite.macromolecule.structure.store`. **Opening a model never runs
code.** A model file can be mailed to a colleague, and the worst a bad one can do is fail to
load.

What the archive carries, and why each part has to be there:

- **the recipe**, so the feature pool can be rebuilt column for column on rows the model has
  never seen;
- **the fill values, the kept columns and the scaling**, because a feature means nothing
  without the arithmetic that put it on a scale;
- **the estimator as arrays** (:mod:`chemsuite.qsar.export`), checked against scikit-learn
  when it was written;
- **the calibration** (:mod:`chemsuite.qsar.conformal`), so an interval means the same thing
  in a later session as it did at the fit;
- **the training rows' keys and values**, because refining a model means refitting its recipe
  on the old rows plus the new ones — a frozen estimator cannot learn, and a model that
  cannot say what it was trained on cannot be refined honestly either;
- **the lineage**, one line per fit or refit, so a model that has been round the loop three
  times says so.

The training keys and the calibration rows are the model's data. Sharing a model shares them,
which the user guide says in as many words.
"""

from __future__ import annotations

import io
import json
import zipfile
from dataclasses import dataclass, field, replace
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np

import chemsuite

from .apply import ModelError, Weights
from .apply import predicted as apply_predicted
from .apply import values as apply_values
from .conformal import Reliability
from .prepare import Prepare
from .recipe import Recipe

#: What a saved model is called. chemsuite QSAR model.
SUFFIX = ".csqm"
FORMAT = 1
_FORMAT_KEY = "_format"
_META_KEY = "_meta"
#: Arrays that are not part of the estimator's own weights.
_OWN = ("kept", "names", "classes", "train_keys", "train_y")


@dataclass(frozen=True)
class Training:
    """Which rows a model was fitted on, so it can be fitted again with more.

    *keys* identify a compound — an InChIKey, where the table could compute one — and are how
    a refit tells an updated measurement from a new compound.
    """

    keys: np.ndarray  # fixed-width text, one per training row
    y: np.ndarray  # the value each row was fitted on (a class index for a classification)

    def __post_init__(self) -> None:
        if len(self.keys) != len(self.y):
            raise ValueError("every training row needs a key and a value")

    def __len__(self) -> int:
        return len(self.y)


@dataclass(frozen=True)
class Prediction:
    """What a model says about some rows."""

    #: a number per row, or the index of the chosen class
    predicted: np.ndarray
    #: regression: the ends of the interval; otherwise None
    low: np.ndarray | None = None
    high: np.ndarray | None = None
    #: classification: one score column per class, and the prediction set
    scores: np.ndarray | None = None
    sets: np.ndarray | None = None
    #: how far each row sits from the training set, and whether that is too far
    distance: np.ndarray | None = None
    inside: np.ndarray | None = None

    def set_names(self, classes: np.ndarray) -> list[str]:
        """Each row's prediction set as text, for a column: "high", "low, high" or "" when the
        model cannot place the row at this confidence."""
        if self.sets is None:
            raise ValueError("this prediction has no class sets")
        return [", ".join(classes[row].tolist()) for row in self.sets]


@dataclass(frozen=True)
class Model:
    """A fitted model that predicts without scikit-learn."""

    recipe: Recipe
    weights: Weights
    prepare: Prepare
    kept: np.ndarray
    names: tuple[str, ...]
    candidate: str
    params: dict[str, Any]
    scores: dict[str, float]
    reliability: Reliability
    training: Training
    classes: np.ndarray | None = None
    leaderboard: tuple[dict[str, Any], ...] = ()
    lineage: tuple[dict[str, Any], ...] = ()
    created: str = ""
    chemsuite_version: str = ""
    objective: str = ""
    label: str = ""

    @property
    def task(self) -> str:
        return self.recipe.task

    @property
    def features(self) -> int:
        return len(self.kept)

    @property
    def title(self) -> str:
        """What to call this model in a panel or a file dialog."""
        return self.label or f"{self.recipe.target} ({self.candidate})"

    def kept_names(self) -> tuple[str, ...]:
        """The pool names of the columns this model actually uses."""
        return tuple(self.names[i] for i in self.kept.tolist())

    # -- predicting ------------------------------------------------------------------------
    def design(self, pool: np.ndarray) -> np.ndarray:
        """The matrix the estimator sees, from a feature pool built on this model's recipe."""
        pool = np.asarray(pool, dtype=np.float64)
        if pool.shape[1] != len(self.names):
            raise ModelError(
                f"the pool has {pool.shape[1]} columns; this model was fitted on "
                f"{len(self.names)}. Rebuild it from the model's own recipe."
            )
        return self.prepare.design(pool, self.kept)

    def predict(self, pool: np.ndarray) -> Prediction:
        """Everything the model has to say about these rows: the prediction, how sure it is,
        and whether the rows are the kind of thing it knows about."""
        design = self.design(pool)
        said = apply_values(self.weights, design)
        chosen = apply_predicted(self.weights, design)
        distance = self.reliability.distance(design)
        inside = distance <= self.reliability.domain_edge
        if self.task == "regression":
            low, high = self.reliability.interval(said, design)
            return Prediction(chosen, low=low, high=high, distance=distance, inside=inside)
        return Prediction(
            chosen,
            scores=said,
            sets=self.reliability.class_sets(said),
            distance=distance,
            inside=inside,
        )

    # -- the archive -----------------------------------------------------------------------
    def with_lineage(self, note: str, scores: dict[str, float], rows: int) -> Model:
        """This model with one more line of history."""
        entry = {"date": date.today().isoformat(), "note": note, "rows": rows, "scores": scores}
        return replace(self, lineage=(*self.lineage, entry))

    def as_meta(self) -> dict[str, Any]:
        return {
            "format": FORMAT,
            "created": self.created or date.today().isoformat(),
            "chemsuite": self.chemsuite_version or chemsuite.__version__,
            "label": self.label,
            "recipe": self.recipe.as_json(),
            "candidate": self.candidate,
            "params": _plain(self.params),
            "family": self.weights.family,
            "link": self.weights.link,
            "weights": _plain(self.weights.meta),
            "scores": _plain(self.scores),
            "objective": self.objective,
            "leaderboard": [_plain(entry) for entry in self.leaderboard],
            "reliability": self.reliability.as_json(),
            "lineage": [_plain(entry) for entry in self.lineage],
            "training": {"rows": len(self.training)},
        }

    def pack(self) -> bytes:
        """The whole model as the bytes of a ``.csqm`` file."""
        arrays: dict[str, np.ndarray] = {
            _FORMAT_KEY: np.array(FORMAT, dtype=np.int32),
            _META_KEY: np.frombuffer(json.dumps(self.as_meta()).encode(), dtype=np.uint8),
            "kept": np.asarray(self.kept, dtype=np.int64),
            "names": np.asarray(self.names, dtype=np.str_),
            "train_keys": np.asarray(self.training.keys, dtype=np.str_),
            "train_y": np.asarray(self.training.y, dtype=np.float64),
        }
        if self.classes is not None:
            arrays["classes"] = np.asarray(self.classes, dtype=np.str_)
        arrays.update(self.prepare.as_arrays())
        arrays.update(self.reliability.as_arrays())
        for name, value in self.weights.arrays.items():
            arrays[f"w_{name}"] = np.asarray(value)
        buffer = io.BytesIO()
        np.savez(buffer, **arrays)
        return buffer.getvalue()

    def to_bytes(self) -> bytes:
        """:meth:`pack`, under the name :class:`chemsuite.core.assets.Packable` asks for, so a
        document can keep a model as an asset and hand back the decoded model rather than the
        bytes it was stored as."""
        return self.pack()

    @property
    def nbytes(self) -> int:
        """Roughly what the model holds in memory, for the document's asset accounting."""
        own = sum(
            value.nbytes
            for value in (*self.weights.arrays.values(), *self.reliability.as_arrays().values())
        )
        return own + self.kept.nbytes + self.training.y.nbytes + self.prepare.fill.nbytes

    def save(self, path: Path) -> Path:
        """Write the model to *path*, giving it the ``.csqm`` suffix if it has none."""
        path = Path(path)
        if not path.suffix:
            path = path.with_suffix(SUFFIX)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(self.pack())
        return path


def _plain(value: Any) -> Any:
    """Numpy scalars and tuples as the plain types JSON can hold."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_plain(item) for item in value.tolist()]
    if isinstance(value, float) and (value != value or value in (float("inf"), float("-inf"))):
        return str(value)  # JSON has no NaN or infinity; they come back as text
    return value


def _number(value: Any) -> float:
    """A score read back from the manifest: a number, or NaN where it was not one."""
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def unpack(data: bytes) -> Model:
    """A model from the bytes of a ``.csqm`` file. Never unpickles."""
    try:
        with np.load(io.BytesIO(data), allow_pickle=False) as npz:
            arrays = {name: npz[name] for name in npz.files}
    except (ValueError, OSError, zipfile.BadZipFile) as error:
        raise ModelError(f"not a chemsuite model: {error}") from error
    version = int(np.ravel(arrays.pop(_FORMAT_KEY, np.array(0)))[0])
    if version < 1:
        raise ModelError("not a chemsuite model: it has no format number")
    if version > FORMAT:
        raise ModelError(f"saved by a newer chemsuite (format {version}; this one reads {FORMAT})")
    if _META_KEY not in arrays:
        raise ModelError("the model has no manifest")
    meta = json.loads(arrays[_META_KEY].tobytes().decode())
    missing = [name for name in ("kept", "names", "train_keys", "train_y") if name not in arrays]
    if missing:
        raise ModelError(f"the model is missing {', '.join(missing)}")
    weights = Weights(
        meta["family"],
        meta["link"],
        {name[2:]: value for name, value in arrays.items() if name.startswith("w_")},
        meta.get("weights", {}),
    )
    scores = {key: _number(value) for key, value in meta.get("scores", {}).items()}
    return Model(
        recipe=Recipe.from_json(meta["recipe"]),
        weights=weights,
        prepare=Prepare.from_arrays(arrays),
        kept=np.asarray(arrays["kept"], dtype=np.int64),
        names=tuple(arrays["names"].tolist()),
        candidate=meta["candidate"],
        params=meta.get("params", {}),
        scores=scores,
        reliability=Reliability.from_saved(arrays, meta["reliability"]),
        training=Training(arrays["train_keys"], np.asarray(arrays["train_y"], dtype=np.float64)),
        classes=arrays.get("classes"),
        leaderboard=tuple(meta.get("leaderboard", ())),
        lineage=tuple(meta.get("lineage", ())),
        created=meta.get("created", ""),
        chemsuite_version=meta.get("chemsuite", ""),
        objective=meta.get("objective", ""),
        label=meta.get("label", ""),
    )


def load(path: Path) -> Model:
    """The model saved at *path*."""
    return unpack(Path(path).read_bytes())


@dataclass(frozen=True)
class Lineage:
    """One line of a model's history, as the panel shows it."""

    date: str
    note: str
    rows: int
    scores: dict[str, float] = field(default_factory=dict)
