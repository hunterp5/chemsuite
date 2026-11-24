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
"""Export for Python: a pickle chemsuite writes and never reads.

A chemist with a model worth keeping often wants it in a notebook next to pandas, and in that
world a model is a pickled scikit-learn estimator. So chemsuite writes one — and reads none.
The direction matters: unpickling executes whatever the file says, so a model file passed
between people is a way to run code on the recipient's machine. chemsuite's own format
(:mod:`chemsuite.qsar.model`) is arrays and JSON for exactly that reason, and Apply Model…
accepts only that.

Two things go in the folder beside the pickle. A **manifest** naming the Python, numpy and
scikit-learn versions it was written under, because a pickle is only readable by a close
enough combination and the error when it is not is unhelpful. And a short **script** showing
how to featurize a SMILES the same way the model expects, since a bare estimator is useless
without the recipe that produced its columns — the pickle holds the estimator, not the pool.

The estimator itself is only in memory while the document that fitted it is open, so this is
offered for a model fitted in this session; a model reopened from a saved document is arrays,
which is all the archive keeps, and the honest answer is to say so rather than to refit
something slightly different behind the user's back.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import chemsuite
from chemsuite.qsar.model import SUFFIX

from .part import QsarModels

#: fitted estimators from this session, by asset id — never saved, never reloaded
_FITTED: dict[tuple[int, str], object] = {}

SCRIPT = '''"""Apply the exported model to new structures.

The pickle holds the fitted estimator only. These are the columns it expects, in this order,
built exactly as chemsuite built them; feed it anything else and the numbers are meaningless.
"""

import joblib
import numpy as np

model = joblib.load({pickle!r})
# The features, in order. See features.json beside this script for the full list.
FEATURES = {features}


def design(rows):
    """rows: one dict of {{feature name: value}} per compound."""
    x = np.array([[float(row.get(name, np.nan)) for name in FEATURES] for row in rows])
    fill = np.array({fill})
    holes = ~np.isfinite(x)
    x[holes] = np.take(fill, np.where(holes)[1])
    {scaling}
    return x


# print(model.predict(design([...])))
'''


def remember(part: QsarModels, asset: str, estimator: object) -> None:
    """Keep the fitted estimator for this session, so it can be exported."""
    _FITTED[(id(part), asset)] = estimator


def fitted_estimator(part: QsarModels, asset: str) -> object | None:
    return _FITTED.get((id(part), asset))


def exportable(part: QsarModels, asset: str) -> str:
    """ "" when the model can be exported, or why it cannot."""
    if fitted_estimator(part, asset) is not None:
        return ""
    return (
        "This model was reopened from a saved document, which keeps it as arrays rather than "
        "as a scikit-learn object, so there is nothing to pickle. Fit it again in this "
        f"session to export it, or use Save Model… for the portable {SUFFIX} file."
    )


def export_for_python(part: QsarModels, asset: str, path: Path) -> Path:
    """Write the pickle, its manifest, its feature list and a usage script. Returns the pickle."""
    import joblib
    import sklearn

    estimator = fitted_estimator(part, asset)
    if estimator is None:
        raise ValueError(exportable(part, asset))
    model = part.model(asset)
    path = Path(path)
    if path.suffix != ".joblib":
        path = path.with_suffix(".joblib")
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(estimator, path)

    features = list(model.kept_names())
    manifest = {
        "written_by": f"chemsuite {chemsuite.__version__}",
        "pickle": path.name,
        "warning": (
            "Loading a pickle runs code. Only load this file if you trust where it came "
            "from. chemsuite itself never reads it back."
        ),
        "needs": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scikit-learn": sklearn.__version__,
            "joblib": joblib.__version__,
        },
        "model": {
            "target": model.recipe.target,
            "task": model.recipe.task,
            "candidate": model.candidate,
            "params": model.params,
            "classes": None if model.classes is None else model.classes.tolist(),
            "scores": model.scores,
        },
        "features": features,
    }
    path.with_suffix(".json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    _write_script(path, model, features)
    return path


def _write_script(path: Path, model, features: list[str]) -> None:
    prepare = model.prepare
    kept = model.kept
    fill = np.asarray(prepare.fill, dtype=np.float64)[kept]
    scaling = "# this model's features were used unscaled"
    if prepare.scale is not None and prepare.center is not None:
        center = np.asarray(prepare.center, dtype=np.float64)
        scale = np.asarray(prepare.scale, dtype=np.float64)
        scaling = (
            f"x = (x - np.array({_as_list(center)})) / np.array({_as_list(scale)})"
            "  # the same scaling the fit used"
        )
    path.with_name(f"{path.stem}_apply.py").write_text(
        SCRIPT.format(
            pickle=path.name,
            features=_as_list(features),
            fill=_as_list(fill),
            scaling=scaling,
        ),
        encoding="utf-8",
    )


def _as_list(values) -> str:
    """A short python literal for the script, on one line."""
    if isinstance(values, np.ndarray):
        return repr([round(float(v), 10) for v in values.tolist()])
    return repr(list(values))
