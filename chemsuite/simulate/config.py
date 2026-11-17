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

"""Where OpenMM's interpreter is.

Two paths come from ``Settings ▸ External Tools``: the Python of a native OpenMM environment
(``openmm_python``) and, on Windows, the Python of one inside WSL (``openmm_wsl_python``, a Linux
path). The native one is used when it is set; the WSL one otherwise. Neither is looked up on
``PATH``: the ``python`` there is almost never the OpenMM environment.
"""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from chemsuite.core.jobs.lines import is_app_python
from chemsuite.core.settings import Settings

PYTHON_TOOL = "openmm_python"
WSL_TOOL = "openmm_wsl_python"
TOOLS = frozenset({PYTHON_TOOL, WSL_TOOL})
#: Settings ▸ External Tools ▸ WSL: the wsl.exe that runs a WSL interpreter
WSL_EXE_TOOL = "wsl"


@dataclass(frozen=True)
class Interpreter:
    python: str
    wsl: bool = False  # a Linux path, run through wsl.exe
    launcher: str = "wsl.exe"  # the wsl.exe that runs it (a WSL interpreter only)


def wsl_launcher(settings: Settings, text: str = "") -> str:
    """The wsl.exe set in Settings (*text* wins while a dialog is open), else ``wsl.exe``
    found on PATH."""
    return (text or "").strip() or settings.tool_paths.get(WSL_EXE_TOOL, "").strip() or "wsl.exe"


def native_python(settings: Settings, text: str = "") -> str:
    """The configured native interpreter (*text* wins while a dialog is open), or ``""``."""
    raw = (text or "").strip() or settings.tool_paths.get(PYTHON_TOOL, "").strip()
    if not raw:
        return ""
    path = Path(raw).expanduser()
    if path.is_file():
        return str(path)
    return shutil.which(raw) or ""


def wsl_python(settings: Settings, text: str = "") -> str:
    """The configured WSL interpreter (a Linux path, not checked here), or ``""``; never off
    Windows."""
    if sys.platform != "win32":
        return ""
    return linux_path((text or "").strip() or settings.tool_paths.get(WSL_TOOL, "").strip())


def linux_path(path: str) -> str:
    """*path* as WSL sees it: Browse… gives a WSL file as ``\\\\wsl.localhost\\Ubuntu\\home\\…``
    (or ``\\\\wsl$\\…``), which ``wsl.exe -e`` cannot run, so the share and distribution are
    dropped (``/home/…``). Anything else is returned as it is."""
    unc = path.replace("\\", "/")
    for share in ("//wsl.localhost/", "//wsl$/"):
        if unc.lower().startswith(share):
            rest = unc[len(share) :].split("/", 1)
            return "/" + rest[1] if len(rest) == 2 else path
    return path


def interpreter(
    settings: Settings, native_text: str = "", wsl_text: str = ""
) -> Interpreter | None:
    """The interpreter to run OpenMM in: the native one if set, else the WSL one."""
    native = native_python(settings, native_text)
    if native:
        return Interpreter(native)
    wsl = wsl_python(settings, wsl_text)
    return Interpreter(wsl, wsl=True, launcher=wsl_launcher(settings)) if wsl else None


def ready_problem(settings: Settings, native_text: str = "", wsl_text: str = "") -> str:
    """Why OpenMM cannot run, or ``""`` when it can be tried. Starts no process."""
    found = interpreter(settings, native_text, wsl_text)
    if found is None:
        return "OpenMM is not set up: set its Python under Settings ▸ External Tools"
    if not found.wsl and is_app_python(found.python):
        return "That interpreter is the one running chemsuite: OpenMM needs its own environment"
    return ""
