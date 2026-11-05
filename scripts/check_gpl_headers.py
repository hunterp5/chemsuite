#!/usr/bin/env python3
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
"""Check (or, with --fix, add) the GPL header on every first-party Python file."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = ("chemsuite", "tests", "scripts")
MARKER = "This file is part of chemsuite."

HEADER = """\
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
"""


def python_files() -> list[Path]:
    files: list[Path] = []
    for name in SOURCE_DIRS:
        files.extend(sorted((ROOT / name).rglob("*.py")))
    return files


def add_header(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    if text.startswith("#!"):
        shebang, _, rest = text.partition("\n")
        path.write_text(f"{shebang}\n{HEADER}{rest}", encoding="utf-8")
    else:
        sep = "\n" if text.strip() else ""
        path.write_text(f"{HEADER}{sep}{text}", encoding="utf-8")


def main(argv: list[str]) -> int:
    fix = "--fix" in argv
    missing = [p for p in python_files() if MARKER not in p.read_text(encoding="utf-8")[:400]]
    if fix:
        for path in missing:
            add_header(path)
        print(f"added header to {len(missing)} file(s)")
        return 0
    if missing:
        print(
            "Missing GPL header (run: python scripts/check_gpl_headers.py --fix):", file=sys.stderr
        )
        for path in missing:
            print(f"  {path.relative_to(ROOT).as_posix()}", file=sys.stderr)
        return 1
    print(f"ok: {len(python_files())} Python files carry the GPL header")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
