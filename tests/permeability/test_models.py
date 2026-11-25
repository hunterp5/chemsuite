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

"""The shipped models: the file is the pinned one, and opening it never runs code.

A model refitted without updating its pin fails here rather than shipping quietly — the
arrangement ``tests/ionisation/test_weights.py`` uses for the pKa weights. Where an archive is
not built in this checkout the pin check skips, because the training data is not committed;
what must hold either way is that the package and the pins agree about which endpoints exist.
"""

from __future__ import annotations

import pickle
import zipfile

import numpy as np
import pytest

from chemsuite.permeability import ENDPOINTS, Endpoint, by_id
from chemsuite.permeability import models as shipped


def ids() -> list[str]:
    return [endpoint.id for endpoint in ENDPOINTS]


# -- what holds whether or not the archives are built ----------------------------------------
def test_every_endpoint_has_a_pin_recorded_for_it() -> None:
    assert sorted(shipped.PINS) == sorted(ids())


def test_every_endpoint_names_a_distinct_archive() -> None:
    files = [endpoint.file for endpoint in ENDPOINTS]
    assert len(set(files)) == len(files)
    assert all(file.endswith(".csqm") for file in files)


def test_an_unknown_endpoint_is_an_error() -> None:
    with pytest.raises(KeyError):
        by_id("no_such_endpoint")


def test_a_model_that_is_not_installed_says_so_rather_than_guessing(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(shipped, "DATA", tmp_path)
    endpoint = ENDPOINTS[0]
    assert not shipped.installed(endpoint)
    with pytest.raises(FileNotFoundError):
        shipped.model(endpoint)


# -- what holds for an archive that is built -------------------------------------------------
@pytest.mark.parametrize("endpoint_id", ids())
def test_the_shipped_archive_is_the_pinned_one(endpoint_id: str) -> None:
    endpoint = by_id(endpoint_id)
    if not shipped.installed(endpoint):
        pytest.skip(f"{endpoint.file} is not built in this checkout")
    assert shipped.mismatch(endpoint) == ""


@pytest.mark.parametrize("endpoint_id", ids())
def test_the_shipped_model_loads_and_matches_its_endpoint(endpoint_id: str) -> None:
    endpoint = by_id(endpoint_id)
    if not shipped.installed(endpoint):
        pytest.skip(f"{endpoint.file} is not built in this checkout")
    model = shipped.model(endpoint)

    assert model.recipe.task == endpoint.task
    assert model.recipe.target == endpoint.id
    assert model.features > 0
    if endpoint.task == "classification":
        assert model.classes is not None and model.classes.size >= 2
    else:
        assert model.classes is None


@pytest.mark.parametrize("endpoint_id", ids())
def test_the_same_model_is_read_once_per_process(endpoint_id: str) -> None:
    endpoint = by_id(endpoint_id)
    if not shipped.installed(endpoint):
        pytest.skip(f"{endpoint.file} is not built in this checkout")
    assert shipped.model(endpoint) is shipped.model(endpoint)


def test_reading_a_model_never_unpickles(tmp_path, monkeypatch) -> None:
    """The property that makes a shipped model safe: a file with a pickle in it is refused,
    not executed. Written here as well as in tests/qsar, because this is the package that
    ships archives to every user."""
    endpoint = Endpoint("x", "X", "X", "regression", "", "x.csqm", "")
    monkeypatch.setattr(shipped, "DATA", tmp_path)
    path = tmp_path / "x.csqm"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("evil.npy", pickle.dumps({"owned": True}))

    assert shipped.installed(endpoint)
    with pytest.raises((ValueError, OSError, KeyError, zipfile.BadZipFile)):
        shipped.model(endpoint)


@pytest.mark.parametrize("endpoint_id", ids())
def test_a_pinned_archive_is_a_plain_npz_of_arrays(endpoint_id: str) -> None:
    endpoint = by_id(endpoint_id)
    if not shipped.installed(endpoint):
        pytest.skip(f"{endpoint.file} is not built in this checkout")
    with np.load(shipped.path_for(endpoint), allow_pickle=False) as archive:
        assert archive.files, "an archive with no arrays in it is not a model"
