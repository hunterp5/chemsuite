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
"""The start-up OpenGL check: software OpenGL on Windows, 3D disabled with a reason."""

from __future__ import annotations

import json

import pytest

from chemsuite.app import graphics
from chemsuite.core.settings import valid

from .conftest import gui

pytestmark = gui


@pytest.fixture(autouse=True)
def no_problem():
    graphics.set_problem(None)
    yield
    graphics.set_problem(None)


class Child:
    """A stand-in for the check's child process."""

    def __init__(self, exitcode: int | None) -> None:
        self.exitcode = exitcode
        self.killed = False

    def join(self, _timeout: float) -> None:
        pass

    def kill(self) -> None:
        self.killed = True


def test_only_windows_ever_draws_in_software(tmp_path) -> None:
    for mode in ("auto", "software"):
        assert not graphics.use_software(mode, tmp_path, platform="darwin")
        assert not graphics.use_software(mode, tmp_path, platform="linux")
    assert graphics.use_software("software", tmp_path, platform="win32")
    assert not graphics.use_software("off", tmp_path, platform="win32")


@pytest.mark.parametrize("exitcode", [1, -11])  # an old OpenGL, or a driver that crashed
def test_a_driver_without_3d_switches_windows_to_software(tmp_path, monkeypatch, exitcode) -> None:
    monkeypatch.setattr(graphics, "_start", lambda cache: Child(exitcode))
    assert graphics.use_software("auto", tmp_path, platform="win32")
    cache = json.loads((tmp_path / graphics.CACHE_FILE).read_text(encoding="utf-8"))
    assert cache == {"driver_ok": False}


def test_a_hung_check_counts_as_no_driver(tmp_path, monkeypatch) -> None:
    child = Child(None)
    monkeypatch.setattr(graphics, "_start", lambda cache: child)
    assert graphics.use_software("auto", tmp_path, platform="win32")
    assert child.killed


def test_a_good_driver_keeps_the_gpu(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(graphics, "_start", lambda cache: Child(0))
    assert not graphics.use_software("auto", tmp_path, platform="win32")


def test_a_known_answer_starts_at_once_and_is_refreshed(tmp_path, monkeypatch) -> None:
    cache = tmp_path / graphics.CACHE_FILE
    cache.write_text('{"driver_ok": false}', encoding="utf-8")
    started = []
    monkeypatch.setattr(graphics, "_start", lambda path: started.append(path) or Child(0))
    assert graphics.use_software("auto", tmp_path, platform="win32")  # not waiting for Child(0)
    assert started == [cache]


@pytest.mark.parametrize(
    ("version", "reason"),
    [
        (None, graphics.NO_OPENGL),
        ((1, 1, False), "OpenGL 1.1; 3D needs 3.3"),  # Windows' own fallback
        ((3, 0, False), "OpenGL 3.0; 3D needs 3.3"),
        ((3, 3, False), None),
        ((4, 1, False), None),  # macOS
        ((3, 0, True), None),  # OpenGL ES 3
        ((2, 0, True), "OpenGL 2.0; 3D needs 3.0"),
    ],
)
def test_the_check_says_why_3d_cannot_draw(monkeypatch, version, reason) -> None:
    monkeypatch.setattr(graphics, "opengl_version", lambda: version)
    found = graphics.check()
    assert found is None if reason is None else reason in found
    assert graphics.problem() == found


def test_turned_off_in_settings(tmp_path) -> None:
    assert not graphics.prepare("off", tmp_path)
    assert graphics.problem() == graphics.TURNED_OFF
    assert graphics.plots_problem() == graphics.TURNED_OFF
    assert valid("performance.graphics", "software") and not valid("performance.graphics", "gpu")


def test_the_check_runs_as_its_own_process(tmp_path) -> None:
    """The real child process: whatever this machine's answer, it is saved for the next start."""
    cache = tmp_path / graphics.CACHE_FILE
    child = graphics._start(cache)
    child.join(60)
    saved = json.loads(cache.read_text(encoding="utf-8"))["driver_ok"]
    assert child.exitcode == (0 if saved else 1)
