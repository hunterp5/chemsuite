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

"""Architecture rules that keep chemsuite's layers apart.

``chemsuite.core`` (and the ported domain packages) must never depend on Qt or on the
``chemsuite.app`` shell, so the table, files and computations stay testable without a GUI and
reusable from scripts. Tools may use the app only through ``ToolContext``.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "chemsuite"

# Packages that must stay free of Qt and of the GUI shell.
QT_FREE = (
    "core",
    "descriptors",
    "plotting",
    "conformers",
    "molview",
    "analysis",
    "sar",
    "ionisation",
    "sketch",
    "retro",
    "protein",
)
FORBIDDEN_FOR_QT_FREE = ("PySide6", "shiboken6", "pyqtgraph", "chemsuite.app", "chemsuite.tools")


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            found.add(node.module)
    return found


def test_qt_free_packages_do_not_import_qt_or_the_app() -> None:
    offenders: list[str] = []
    for package in QT_FREE:
        for path in sorted((PACKAGE / package).rglob("*.py")):
            for module in _imported_modules(path):
                if any(
                    module == bad or module.startswith(bad + ".") for bad in FORBIDDEN_FOR_QT_FREE
                ):
                    offenders.append(f"{path.relative_to(PACKAGE.parent).as_posix()}: {module}")
    assert not offenders, "Qt-free package imports Qt or the GUI:\n" + "\n".join(offenders)


def test_importing_core_does_not_load_qt() -> None:
    """Catches indirect imports the static scan cannot see (e.g. via a third-party module)."""
    code = (
        "import pkgutil, importlib, sys\n"
        "import chemsuite.core as core\n"
        "for m in pkgutil.walk_packages(core.__path__, core.__name__ + '.'):\n"
        "    importlib.import_module(m.name)\n"
        "bad = sorted(n for n in sys.modules if n.split('.')[0] in ('PySide6', 'pyqtgraph'))\n"
        "print(','.join(bad))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "", f"core loaded Qt modules: {result.stdout.strip()}"
