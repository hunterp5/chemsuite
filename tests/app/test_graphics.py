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
"""The start-up graphics check: software drawing on request, 3D greyed out with a reason."""

from __future__ import annotations

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


def test_automatic_leaves_the_graphics_to_the_system() -> None:
    environ: dict[str, str] = {}
    graphics.prepare("auto", environ)
    assert environ == {} and graphics.problem() is None


def test_software_asks_direct3d_and_mesa_for_their_software_renderers() -> None:
    environ: dict[str, str] = {}
    graphics.prepare("software", environ)
    assert environ == graphics.SOFTWARE_ENV
    assert graphics.problem() is None


def test_turned_off_in_settings() -> None:
    graphics.prepare("off", {})
    assert graphics.problem() == graphics.TURNED_OFF
    assert graphics.check() == graphics.TURNED_OFF  # the check never overrides it
    assert valid("performance.graphics", "software") and not valid("performance.graphics", "gpu")


def test_graphics_that_do_not_start_grey_3d_out(monkeypatch) -> None:
    monkeypatch.setattr(graphics, "_starts", lambda: False)
    assert graphics.check() == graphics.NO_GPU == graphics.problem()


def test_graphics_that_start_leave_3d_on(monkeypatch) -> None:
    monkeypatch.setattr(graphics, "_starts", lambda: True)
    assert graphics.check() is None


def test_the_check_asks_qt() -> None:
    """The real probe answers without raising (offscreen here: the answer may be no)."""
    assert graphics._starts() in (True, False)
    assert graphics._starts("win32") is True  # Direct3D falls back to WARP by itself
