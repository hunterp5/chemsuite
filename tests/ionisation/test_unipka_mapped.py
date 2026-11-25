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

"""The pKa parameters mapped from one cache file that every worker process shares."""

from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pytest

from chemsuite.ionisation.unipka import mapped


def _archive(path: Path, scale: float = 1.0) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(0)
    params = {
        "f0/a.weight": rng.normal(size=(3, 4)) * scale,
        "f0/a.bias": rng.normal(size=4) * scale,
        "f1/a.weight": rng.normal(size=(2, 2)) * scale,
    }
    np.savez_compressed(path, mean=np.float64(5.0), **params)
    return params


def test_the_parameters_are_the_archives_as_float32(tmp_path) -> None:
    source = tmp_path / "w.npz"
    params = _archive(source)
    found = mapped.load(source, tmp_path / "cache")
    assert found is not None and set(found) == {"f0", "f1"}
    for key, values in params.items():
        fold, name = key.split("/", 1)
        got = found[fold][name]
        assert got.dtype == np.float32 and not got.flags.writeable
        np.testing.assert_array_equal(got, values.astype(np.float32))


def test_the_cache_is_written_once_and_then_mapped(tmp_path, monkeypatch) -> None:
    source, cache = tmp_path / "w.npz", tmp_path / "cache"
    _archive(source)
    assert mapped.load(source, cache) is not None
    monkeypatch.setattr(mapped, "_write", lambda *a: pytest.fail("written again"))
    assert mapped.load(source, cache) is not None
    assert not list(cache.glob("*.lock")) and not list(cache.glob("*.tmp"))


def test_a_changed_archive_gets_a_new_cache(tmp_path) -> None:
    source, cache = tmp_path / "w.npz", tmp_path / "cache"
    _archive(source)
    mapped.load(source, cache)
    params = _archive(source, scale=2.0)
    os.utime(source, ns=(time.time_ns(), time.time_ns() + 10**9))  # a different mtime
    found = mapped.load(source, cache)
    assert found is not None
    np.testing.assert_array_equal(found["f0"]["a.bias"], params["f0/a.bias"].astype(np.float32))


def test_a_crashed_writers_lock_is_taken_over(tmp_path) -> None:
    source, cache = tmp_path / "w.npz", tmp_path / "cache"
    _archive(source)
    cache.mkdir()
    lock = cache / f"{mapped._stem(source)}.lock"
    lock.touch()
    old = time.time() - mapped.STALE_LOCK_S - 10
    os.utime(lock, (old, old))
    assert mapped.load(source, cache) is not None
    assert not lock.exists()


def test_a_cache_that_cannot_be_written_falls_back(tmp_path, monkeypatch) -> None:
    source = tmp_path / "w.npz"
    _archive(source)
    blocked = tmp_path / "file"
    blocked.write_text("")  # a file where the folder would be
    assert mapped.load(source, blocked / "cache") is None


def test_an_incomplete_cache_is_rewritten(tmp_path) -> None:
    source, cache = tmp_path / "w.npz", tmp_path / "cache"
    _archive(source)
    mapped.load(source, cache)
    (data,) = cache.glob("*.f32")
    data.write_bytes(data.read_bytes()[:8])  # truncated, as by a crash
    found = mapped.load(source, cache)
    assert found is not None and found["f0"]["a.weight"].shape == (3, 4)


def test_the_predictor_maps_the_shipped_weights(tmp_path, monkeypatch) -> None:
    from chemsuite.ionisation.unipka import Predictor, weights

    monkeypatch.setenv("CHEMSUITE_CACHE_DIR", str(tmp_path))
    shared = Predictor(weights.PATH, n_3d=1)
    own = Predictor(weights.PATH, n_3d=1, shared=False)
    assert not shared.nets[0].p["embed_tokens.weight"].flags.writeable
    assert shared.free_energy("CCC(=O)O") == own.free_energy("CCC(=O)O")
