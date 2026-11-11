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
"""Whether this computer can draw 3D, decided once at start-up.

The Structure Viewer's 3D and the 3D plots draw through Qt's RHI, in the system's own graphics
API: Direct3D 11 on Windows, Metal on macOS, OpenGL on Linux. Remote desktops, virtual
machines, VDI and locked-down or missing drivers are where that can fail:

- **Windows** needs nothing here: when Direct3D cannot use the graphics card, Qt falls back to
  WARP, Windows' built-in software renderer, by itself.
- :func:`prepare` (before the QApplication) applies the ``performance.graphics`` setting, which
  IT can force through the policy file: ``software`` makes Windows use WARP and Linux Mesa's
  software OpenGL from the start; ``off`` disables 3D.
- :func:`check` (after the QApplication) asks Qt whether the system's graphics API can be
  started (Metal or OpenGL; PySide6 cannot start Direct3D from Python). When it cannot,
  :func:`problem` says why, and 3D views and plots are greyed out with that reason instead of
  failing when opened. A failure only found when drawing (Windows) is reported by the view
  itself, through :func:`set_problem`.
"""

from __future__ import annotations

import logging
import os
import sys

from chemsuite.core.settings import GRAPHICS_OFF, GRAPHICS_SOFTWARE

log = logging.getLogger(__name__)

NO_GPU = "3D drawing is not available on this computer's graphics."
TURNED_OFF = "3D drawing is turned off in Settings."
#: environment that asks for software drawing: Qt's widgets on Direct3D, Mesa's OpenGL
SOFTWARE_ENV = {"QSG_RHI_PREFER_SOFTWARE_RENDERER": "1", "LIBGL_ALWAYS_SOFTWARE": "1"}

_problem: str | None = None


def problem() -> str | None:
    """Why 3D cannot be drawn here, or None when it can (as far as is known)."""
    return _problem


def set_problem(text: str | None) -> None:
    global _problem
    _problem = text


def prepare(mode: str, environ=os.environ) -> None:
    """Apply the 3D graphics setting (call before the QApplication exists)."""
    if mode == GRAPHICS_OFF:
        set_problem(TURNED_OFF)
    elif mode == GRAPHICS_SOFTWARE:
        log.info("3D draws in software (Settings)")
        environ.update(SOFTWARE_ENV)


def check() -> str | None:
    """Ask Qt whether the system's graphics API starts; returns (and remembers) why 3D cannot
    be drawn, if so."""
    if _problem is None and not _starts():
        set_problem(NO_GPU)
    log.info("3D graphics: %s", _problem or "available")
    return _problem


def _starts(platform: str = sys.platform) -> bool:
    if platform == "win32":
        return True  # Direct3D 11, with WARP behind it
    from PySide6.QtGui import QRhi

    if platform == "darwin":
        from PySide6.QtGui import QRhiMetalInitParams  # only macOS builds of PySide6 have it

        return QRhi.probe(QRhi.Implementation.Metal, QRhiMetalInitParams())
    from PySide6.QtGui import QOpenGLContext, QRhiGles2InitParams

    if not QOpenGLContext().create():  # no OpenGL at all: the surface below would crash
        return False
    params = QRhiGles2InitParams()
    surface = QRhiGles2InitParams.newFallbackSurface()
    params.fallbackSurface = surface
    try:
        return QRhi.probe(QRhi.Implementation.OpenGLES2, params)
    finally:
        surface.destroy()
