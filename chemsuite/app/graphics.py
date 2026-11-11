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

The Structure Viewer's 3D and the 3D plots draw with OpenGL. Remote desktops, virtual
machines, VDI and locked-down or missing drivers often give no usable OpenGL (Windows'
built-in fallback is OpenGL 1.1). So:

1. :func:`prepare` (before the QApplication) chooses the OpenGL. On Windows it asks a short
   child process whether the driver gives :data:`MIN_VERSION`; if not, Qt's bundled software
   OpenGL (``opengl32sw.dll``, slower but complete) is used instead. The answer is cached and
   refreshed in the background for the next start, so only the first start waits for it.
2. :func:`check` (after the QApplication) tests the OpenGL actually in use. When it cannot
   draw 3D, :func:`problem` says why, and 3D views and plots are disabled with that reason
   instead of failing when opened.

The ``performance.graphics`` setting (which IT can force through the policy file) overrides
the automatic choice: ``software`` always uses software OpenGL (Windows), ``off`` disables 3D.
"""

from __future__ import annotations

import json
import logging
import multiprocessing
import sys
from importlib.util import find_spec
from multiprocessing.process import BaseProcess
from pathlib import Path

from chemsuite.core.settings import GRAPHICS_OFF, GRAPHICS_SOFTWARE

log = logging.getLogger(__name__)

MIN_VERSION = (3, 3)  # instanced drawing and the shaders of both 3D views
MIN_ES_VERSION = (3, 0)  # the same on OpenGL ES (some Linux systems)
PROBE_TIMEOUT_S = 20
CACHE_FILE = "graphics.json"

NO_OPENGL = "3D drawing is not available: this computer's graphics give no OpenGL."
OLD_OPENGL = (
    "3D drawing is not available: this computer's graphics give OpenGL {0}.{1}; 3D needs {2}.{3}."
)
TURNED_OFF = "3D drawing is turned off in Settings."
NO_GPU = "3D drawing is not available on this computer's graphics."
NO_PYOPENGL = "3D plots need the PyOpenGL package."

_problem: str | None = None
_refresh: BaseProcess | None = None  # the background probe for the next start


def problem() -> str | None:
    """Why 3D cannot be drawn here, or None when it can (or before :func:`check`)."""
    return _problem


def plots_problem() -> str | None:
    """Why the 3D plots cannot be drawn here (no PyOpenGL, or the graphics), or None."""
    return NO_PYOPENGL if find_spec("OpenGL") is None else _problem


def set_problem(text: str | None) -> None:
    global _problem
    _problem = text


# -- before the QApplication --------------------------------------------------------------------
def prepare(mode: str, cache_folder: Path) -> bool:
    """Choose the OpenGL for this run; True when software OpenGL is used."""
    from PySide6.QtCore import QCoreApplication, Qt

    from .molview.renderer import prepare_opengl

    prepare_opengl()
    if mode == GRAPHICS_OFF:
        set_problem(TURNED_OFF)
    software = use_software(mode, cache_folder)
    if software:
        log.info("using software OpenGL")
        QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_UseSoftwareOpenGL)
    return software


def use_software(mode: str, cache_folder: Path, platform: str = sys.platform) -> bool:
    """Whether to draw with Qt's software OpenGL (only shipped for Windows)."""
    if platform != "win32" or mode == GRAPHICS_OFF:
        return False
    return mode == GRAPHICS_SOFTWARE or not _driver_ok(cache_folder)


def _driver_ok(folder: Path) -> bool:
    """Whether the graphics driver gives MIN_VERSION: cached, else asked of a child process."""
    global _refresh
    cache = folder / CACHE_FILE
    try:
        known = json.loads(cache.read_text(encoding="utf-8"))["driver_ok"]
    except (OSError, ValueError, KeyError, TypeError):
        known = None
    if isinstance(known, bool):
        _refresh = _start(cache)  # for the next start: drivers change
        return known
    child = _start(cache)
    child.join(PROBE_TIMEOUT_S)
    if child.exitcode is None:
        child.kill()
    ok = child.exitcode == 0
    log.info("OpenGL check: the driver %s", "can draw 3D" if ok else f"cannot ({child.exitcode})")
    _write(cache, ok)  # also when the check crashed or hung
    return ok


def _start(cache: Path) -> BaseProcess:
    """The check in a child process (spawned, so it also works in frozen builds)."""
    child = multiprocessing.get_context("spawn").Process(
        target=run_probe, args=(str(cache),), name="opengl-check", daemon=True
    )
    child.start()
    return child


def _write(cache: Path, ok: bool) -> None:
    try:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps({"driver_ok": ok}), encoding="utf-8")
    except OSError:
        log.warning("could not save the OpenGL check", exc_info=True)


def run_probe(cache: str) -> None:
    """The child process: save whether the driver's OpenGL can draw 3D; exit 0 if so, else 1.

    A driver that crashes or hangs here takes only the child down.
    """
    from PySide6.QtGui import QGuiApplication

    app = QGuiApplication([])  # noqa: F841 - OpenGL needs the application to exist
    ok = _usable(opengl_version())
    _write(Path(cache), ok)
    sys.exit(0 if ok else 1)


# -- after the QApplication ---------------------------------------------------------------------
def opengl_version() -> tuple[int, int, bool] | None:
    """The OpenGL Qt gives this process (major, minor, is OpenGL ES), or None when none."""
    from PySide6.QtGui import QOffscreenSurface, QOpenGLContext

    context = QOpenGLContext()
    if not context.create():
        return None
    surface = QOffscreenSurface()
    surface.setFormat(context.format())
    surface.create()
    if not surface.isValid() or not context.makeCurrent(surface):
        return None
    fmt = context.format()
    context.doneCurrent()
    return fmt.majorVersion(), fmt.minorVersion(), context.isOpenGLES()


def _needed(version: tuple[int, int, bool]) -> tuple[int, int]:
    return MIN_ES_VERSION if version[2] else MIN_VERSION


def _usable(version: tuple[int, int, bool] | None) -> bool:
    return version is not None and version[:2] >= _needed(version)


def check() -> str | None:
    """Test the OpenGL in use; returns (and remembers) why 3D cannot be drawn, if so."""
    if _problem is not None:
        return _problem
    version = opengl_version()
    if version is None:
        set_problem(NO_OPENGL)
    elif not _usable(version):
        set_problem(OLD_OPENGL.format(*version[:2], *_needed(version)))
    log.info("OpenGL %s%s", version, f": {_problem}" if _problem else "")
    return _problem
