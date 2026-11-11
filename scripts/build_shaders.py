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
"""Compile the 3D Viewer's shaders to .qsb files (checked in, so users never need this).

    python scripts/build_shaders.py [--all]

Rebuilds only the shaders whose source changed since the last build (``--all``: every one).
Needs Qt's ``qsb`` tool: set ``QSB`` to its path, or have it on PATH (PySide6-Essentials does
not ship it; Qt's installer, Homebrew's ``qt`` or PySide6-Addons do: on Windows, install
PySide6-Addons of the same version in a scratch venv and point ``QSB`` at its ``qsb.exe``).
Each shader becomes SPIR-V, GLSL 150 and 330 (OpenGL), HLSL 5.0 (Direct3D) and MSL 1.2
(Metal).
"""

from __future__ import annotations

import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHADERS = ROOT / "chemsuite" / "app" / "molview" / "shaders"
TARGETS = ["--glsl", "100es,120,150,330", "--hlsl", "50", "--msl", "12"]
MANIFEST = SHADERS / "manifest.json"  # each source's hash when its .qsb was built


def find_qsb() -> str | None:
    found = os.environ.get("QSB") or shutil.which("qsb")
    if found:
        return found
    candidates = glob.glob("/usr/local/Cellar/qt/*/bin/qsb") + glob.glob(
        "/opt/homebrew/Cellar/qt/*/bin/qsb"
    )
    return sorted(candidates)[-1] if candidates else None


def sources() -> list[Path]:
    return sorted(p for p in SHADERS.iterdir() if p.suffix in (".vert", ".frag"))


def digest(path: Path) -> str:
    """The source's hash, the same whatever line endings a checkout gave it."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def hashes() -> dict[str, str]:
    return {source.name: digest(source) for source in sources()}


def main(argv: list[str]) -> int:
    qsb = find_qsb()
    if qsb is None:
        print("qsb not found: set QSB to Qt's shader baker", file=sys.stderr)
        return 1
    # Only changed sources by default: qsb's output differs in bytes between systems, so a full
    # rebuild on the other machine would rewrite every checked-in .qsb for nothing.
    built = {} if "--all" in argv or not MANIFEST.exists() else json.loads(MANIFEST.read_text())
    for source in sources():
        target = source.with_name(source.name + ".qsb")
        if target.exists() and built.get(source.name) == digest(source):
            continue
        subprocess.run([qsb, *TARGETS, "-o", str(target), str(source)], check=True)
        print(target.relative_to(ROOT))
    text = json.dumps(hashes(), indent=1, sort_keys=True) + "\n"
    MANIFEST.write_text(text, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
