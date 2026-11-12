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
reusable from scripts. Tools may use the app only through ``ToolContext``. The Protein
Workbench's Qt layer (``chemsuite.protein.ui``) reaches the table only through its bridge,
``chemsuite.tools.protein``.
"""

from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "chemsuite"

# Packages that must stay free of Qt and of the GUI shell.
QT_FREE = (
    "core",
    "io",
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
#: The Qt part of a Qt-free package, scanned by its own rules below.
WORKBENCH_UI = PACKAGE / "protein" / "ui"


def _imported_modules(path: Path) -> set[str]:
    """Every module *path* imports, relative imports resolved."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    package = path.parent.relative_to(PACKAGE.parent).parts
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            base = list(package[: len(package) - node.level + 1]) if node.level else []
            module = ".".join([*base, *([node.module] if node.module else [])])
            found.add(module)
            found.update(  # names that are modules: ``from .. import ui``
                f"{module}.{alias.name}"
                for alias in node.names
                if _is_module(f"{module}.{alias.name}")
            )
    return found


def _is_module(name: str) -> bool:
    *parents, last = name.split(".")
    folder = PACKAGE.parent.joinpath(*parents)
    # exact names: Windows and macOS would find table.py for "Table"
    return folder.is_dir() and bool({last, f"{last}.py"} & set(os.listdir(folder)))


def _qt_free_files(package: str) -> list[Path]:
    return sorted(
        p for p in (PACKAGE / package).rglob("*.py") if not p.is_relative_to(WORKBENCH_UI)
    )


def test_qt_free_packages_do_not_import_qt_or_the_app() -> None:
    offenders: list[str] = []
    for package in QT_FREE:
        for path in _qt_free_files(package):
            for module in _imported_modules(path):
                if any(
                    module == bad or module.startswith(bad + ".") for bad in FORBIDDEN_FOR_QT_FREE
                ):
                    offenders.append(f"{path.relative_to(PACKAGE.parent).as_posix()}: {module}")
    assert not offenders, "Qt-free package imports Qt or the GUI:\n" + "\n".join(offenders)


#: Domain packages sit above core: core never imports them (chemsuite.io, which reads and
#: writes files, sits above core and the sketcher, whose model lays out reactions it writes).
DOMAIN = tuple(p for p in QT_FREE if p != "core")


def _offenders(paths, forbidden) -> list[str]:
    found = []
    for path in paths:
        for module in _imported_modules(path):
            if any(module == bad or module.startswith(bad + ".") for bad in forbidden):
                found.append(f"{path.relative_to(PACKAGE.parent).as_posix()}: {module}")
    return found


def test_core_imports_no_domain_package() -> None:
    offenders = []
    forbidden = [f"chemsuite.{p}" for p in DOMAIN]
    for path in sorted((PACKAGE / "core").rglob("*.py")):
        offenders += _offenders([path], forbidden)
    assert not offenders, "core imports a package above it:\n" + "\n".join(offenders)


def test_the_app_never_imports_a_tool() -> None:
    """Tools plug into the app (the registry in chemsuite.tools lists them); the app reaches a
    tool only through its Tool declaration, never by importing it."""
    found = []
    for path in sorted((PACKAGE / "app").rglob("*.py")):
        for module in _imported_modules(path):
            if module.startswith("chemsuite.tools."):
                found.append(f"{path.relative_to(PACKAGE.parent).as_posix()}: {module}")
    assert not found, "app imports a tool:\n" + "\n".join(found)


def test_importing_core_does_not_load_qt() -> None:
    """Catches indirect imports the static scan cannot see (e.g. via a third-party module).
    Covers core and the file formats above it."""
    code = (
        "import pkgutil, importlib, sys\n"
        "import chemsuite.core, chemsuite.io\n"
        "for package in (chemsuite.core, chemsuite.io):\n"
        "    for m in pkgutil.walk_packages(package.__path__, package.__name__ + '.'):\n"
        "        importlib.import_module(m.name)\n"
        "bad = sorted(n for n in sys.modules if n.split('.')[0] in ('PySide6', 'pyqtgraph'))\n"
        "print(','.join(bad))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "", f"core loaded Qt modules: {result.stdout.strip()}"


#: The workbench's Qt layer uses the app's shell, widgets and 3D base, never the table window
#: or a tool.
WORKBENCH_UI_FORBIDDEN = ("chemsuite.app.document", "chemsuite.app.table", "chemsuite.tools")
#: Workbench UI files that still read the document's Table, until the workspace engine takes
#: their place (refactor Phase C). A ratchet: this set only shrinks.
WORKBENCH_UI_READS_TABLE = {"panel.py", "session.py"}


def test_the_workbench_ui_reaches_the_table_only_through_the_bridge() -> None:
    paths = sorted(WORKBENCH_UI.rglob("*.py"))
    offenders = _offenders(paths, WORKBENCH_UI_FORBIDDEN)
    for path in paths:
        if path.name not in WORKBENCH_UI_READS_TABLE:
            offenders += _offenders([path], ["chemsuite.core.table"])
    assert not offenders, "the workbench UI reaches past its bridge:\n" + "\n".join(offenders)
    stale = sorted(
        name
        for name in WORKBENCH_UI_READS_TABLE
        if not _offenders([WORKBENCH_UI / name], ["chemsuite.core.table"])
    )
    assert not stale, f"no longer read the Table: drop {stale} from WORKBENCH_UI_READS_TABLE"


def test_only_the_bridge_imports_the_workbench_ui() -> None:
    """The engine (the rest of chemsuite.protein) and the app never import the workbench UI;
    the table side opens it through its bridge, chemsuite.tools.protein."""
    bridge = PACKAGE / "tools" / "protein"
    found = []
    for path in sorted(PACKAGE.rglob("*.py")):
        if not (path.is_relative_to(WORKBENCH_UI) or path.is_relative_to(bridge)):
            found += _offenders([path], ["chemsuite.protein.ui"])
    assert not found, "imports the workbench UI from outside its bridge:\n" + "\n".join(found)
