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

"""External programs chemsuite can use, and where to find them.

A path set in Settings (or forced by the IT policy) wins; otherwise the program is looked up on
``PATH`` under its usual name. Tools come from MolManager's Configure menu; the features that use
most of them arrive later, but their paths can be set now.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass

from .settings import Settings

TEST_TIMEOUT_S = 15


@dataclass(frozen=True)
class ExternalTool:
    id: str
    title: str
    help: str
    names: tuple[str, ...] = ()  # executable names looked up on PATH
    kind: str = "executable"  # "executable", "file" (e.g. a .jar or config) or "folder"
    version_args: tuple[str, ...] = ("--version",)
    windows_only: bool = False


TOOLS: tuple[ExternalTool, ...] = (
    ExternalTool(
        "obabel",
        "Open Babel",
        "obabel, for Confab systematic conformers.",
        ("obabel",),
        version_args=("-V",),
    ),
    ExternalTool(
        "confgen",
        "CONFORGE",
        "CDPKit's confgen, for CONFORGE conformers.",
        ("confgen",),
        version_args=("-V",),
    ),
    ExternalTool(
        "gnina", "Gnina", "Gnina docking; on Windows the Linux binary run through WSL.", ("gnina",)
    ),
    ExternalTool(
        "ambertools",
        "AmberTools",
        "antechamber, for GAFF parameters.",
        ("antechamber",),
        version_args=("-h",),
    ),
    ExternalTool(
        "mafft", "MAFFT", "MAFFT sequence alignment (executable or install folder).", ("mafft",)
    ),
    ExternalTool(
        "aizynthfinder_python",
        "AiZynthFinder Python",
        "Python interpreter with AiZynthFinder installed, for retrosynthesis.",
        ("python3", "python"),
    ),
    ExternalTool(
        "aizynthfinder_config",
        "AiZynthFinder config",
        "AiZynthFinder's config.yml.",
        kind="file",
        version_args=(),
    ),
    ExternalTool(
        "wsl",
        "WSL",
        "wsl.exe, so Windows can run Linux-only tools.",
        ("wsl",),
        version_args=("--status",),
        windows_only=True,
    ),
)
BY_ID = {tool.id: tool for tool in TOOLS}


def visible_tools() -> list[ExternalTool]:
    return [t for t in TOOLS if not t.windows_only or sys.platform == "win32"]


def resolve(tool_id: str, settings: Settings) -> str | None:
    """The configured path, else the first of the tool's names found on ``PATH``."""
    configured = settings.tool_paths.get(tool_id, "").strip()
    if configured:
        return configured
    for name in BY_ID[tool_id].names:
        found = shutil.which(name)
        if found:
            return found
    return None


def test_tool(tool_id: str, path: str) -> tuple[bool, str]:
    """Run the tool's version command (or check the file exists): ``(ok, message)``."""
    from pathlib import Path

    tool = BY_ID[tool_id]
    target = Path(path)
    if tool.kind == "folder":
        return (target.is_dir(), "Folder found." if target.is_dir() else "No such folder.")
    if not target.exists() and shutil.which(path) is None:
        return False, "No such file."
    if tool.kind == "file" or not tool.version_args:
        return True, "File found."
    try:
        done = subprocess.run(
            [path, *tool.version_args],
            capture_output=True,
            text=True,
            timeout=TEST_TIMEOUT_S,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as err:
        return False, f"Could not run it: {err}"
    output = (done.stdout or done.stderr or "").strip().splitlines()
    first = output[0][:200] if output else f"exit code {done.returncode}"
    return True, first
