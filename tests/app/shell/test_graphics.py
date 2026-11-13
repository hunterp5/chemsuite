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

from chemsuite.app.shell import graphics
from chemsuite.core.settings import valid
from tests.gui import gui

pytestmark = gui


@pytest.fixture(autouse=True)
def no_problem():
    graphics.set_problem(None)
    yield
    graphics.set_problem(None)


def test_automatic_leaves_the_graphics_to_the_system() -> None:
    environ: dict[str, str] = {}
    graphics.prepare("auto", environ, platform="linux")
    assert environ == {} and graphics.problem() is None


def test_windows_and_macos_compose_every_window_through_the_rhi() -> None:
    # so docking the first 3D view does not make Qt close and reopen the window
    for platform in ("win32", "darwin"):
        environ: dict[str, str] = {}
        graphics.prepare("auto", environ, platform=platform)
        assert environ == graphics.RHI_WINDOWS
    kept = {"QT_WIDGETS_RHI": "0"}  # a choice made outside the app stands
    graphics.prepare("auto", kept, platform="win32")
    assert kept == {"QT_WIDGETS_RHI": "0"}


def test_software_asks_direct3d_and_mesa_for_their_software_renderers() -> None:
    environ: dict[str, str] = {}
    graphics.prepare("software", environ, platform="linux")
    assert environ == graphics.SOFTWARE_ENV
    assert graphics.problem() is None


def test_turned_off_in_settings() -> None:
    environ: dict[str, str] = {}
    graphics.prepare("off", environ, platform="win32")
    assert environ == {}  # no 3D: the windows draw as they always have
    assert graphics.problem() == graphics.TURNED_OFF
    assert graphics.check() == graphics.TURNED_OFF  # the check never overrides it
    assert valid("performance.graphics", "software") and not valid("performance.graphics", "gpu")


def test_graphics_that_do_not_start_grey_3d_out(monkeypatch) -> None:
    monkeypatch.setattr(graphics, "_starts", lambda: False)
    assert graphics.check() == graphics.NO_GPU == graphics.problem()


def test_graphics_that_start_leave_3d_on(monkeypatch) -> None:
    monkeypatch.setattr(graphics, "_starts", lambda: True)
    assert graphics.check() is None


def test_the_check_asks_qt(qapp) -> None:
    """The real probe answers without raising (offscreen here: the answer may be no)."""
    assert graphics._starts() in (True, False)
    assert graphics._starts("win32") is True  # Direct3D falls back to WARP by itself


def test_the_check_needs_no_metal_outside_macos(qapp, monkeypatch) -> None:
    """PySide6 has QRhiMetalInitParams only on macOS: the Windows and Linux checks must not
    touch it (they did once, and chemsuite would not start there)."""
    import sys
    import types

    from PySide6 import QtGui

    class WithoutMetal(types.ModuleType):  # QtGui as Windows and Linux builds have it
        def __getattr__(self, name: str):
            if name == "QRhiMetalInitParams":
                raise AttributeError(name)
            return getattr(QtGui, name)

    monkeypatch.setitem(sys.modules, "PySide6.QtGui", WithoutMetal("PySide6.QtGui"))
    assert graphics._starts("win32") is True
    assert graphics._starts("linux") in (True, False)
