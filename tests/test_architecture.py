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
"""Guards against architecture drift (see docs/ARCHITECTURE.md).

MolManager did not fail in one step: each of its problems (dialogs reaching into the window,
thousands of ``getattr`` fallbacks, a window of 60 mixins) grew one reasonable-looking change
at a time. These tests make each such step visible:

- **Rules** fail on any new breach: tools use the app only through its tool API.
- **Ratchets** hold today's counts as budgets. A change may lower a count (then lower the
  budget here too, so it stays lowered) but never raise it without editing this file, which a
  reviewer sees.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "chemsuite"
TOOLS = PACKAGE / "tools"


def _sources(folder: Path = PACKAGE) -> dict[str, str]:
    return {
        path.relative_to(ROOT).as_posix(): path.read_text(encoding="utf-8")
        for path in sorted(folder.rglob("*.py"))
    }


def _count(pattern: str, folder: Path = PACKAGE) -> dict[str, int]:
    regex = re.compile(pattern)
    counts = {name: len(regex.findall(text)) for name, text in _sources(folder).items()}
    return {name: n for name, n in counts.items() if n}


def _over(found: dict[str, int], budget: dict[str, int], default: int = 0) -> list[str]:
    return [
        f"{name}: {n} (budget {budget.get(name, default)})"
        for name, n in sorted(found.items())
        if n > budget.get(name, default)
    ]


# -- rules ----------------------------------------------------------------------------------------

#: The app modules a tool may import. ToolContext (chemsuite.app.tools) is the API; the rest are
#: shared widgets. Adding one here is an architecture decision: say why in the commit.
TOOL_APP_IMPORTS = {
    "chemsuite.app.tools",  # Tool, ToolContext: the tool API
    "chemsuite.app.jobs",  # lane names
    "chemsuite.app.dialogs",  # the dialog kit every tool dialog is built on
    "chemsuite.app.dialogs.kit",
    "chemsuite.app.flow_layout",  # compact combo boxes for narrow panels
    "chemsuite.app.glyphs",  # toolbar icons
    "chemsuite.app.settings_hub",  # read-only settings (fingerprint defaults)
    "chemsuite.app.structure_copy",  # Copy as SMILES / molfile
    "chemsuite.app.sketcher",  # the structure editor, opened by tools that draw
    "chemsuite.app.plots.panel",  # plots opened by analysis tools
    "chemsuite.app.viewer.panel",  # the Structure Viewer, opened by its tool
    "chemsuite.app.protein.panel",  # the Protein Viewer, opened by its tool
}


def _app_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.level:  # relative: resolve against the file's package
                parts = path.relative_to(ROOT).with_suffix("").parts[: -node.level]
                names = [".".join((*parts, node.module))]
            else:
                names = [node.module]
        else:
            continue
        found.update(n for n in names if n == "chemsuite.app" or n.startswith("chemsuite.app."))
    return found


def test_tools_use_the_app_only_through_its_tool_api() -> None:
    """MolManager's dialogs reached into the window 271 times; tools here get ToolContext and
    a few shared widgets, nothing else of the app."""
    offenders = [
        f"{path.relative_to(ROOT).as_posix()}: {module}"
        for path in sorted(TOOLS.rglob("*.py"))
        for module in sorted(_app_imports(path))
        if module not in TOOL_APP_IMPORTS
    ]
    assert not offenders, (
        "Tools may import only the tool API and shared widgets (TOOL_APP_IMPORTS); use "
        "ToolContext or add a narrow method to it:\n" + "\n".join(offenders)
    )


def test_tools_never_reach_into_the_window() -> None:
    found = _count(r"ctx\._|\._window\b|main_window", TOOLS)
    assert not found, f"tools reach past ToolContext: {found}"


# -- ratchets -------------------------------------------------------------------------------------

#: Tools reading or writing table columns directly instead of through ToolContext /
#: table.read() / apply_results. Lower these as tools move to the context's helpers.
COLUMN_INTERNALS_IN_TOOLS = {
    "chemsuite/tools/retro/dialog.py": 9,
    "chemsuite/tools/utilities/calculator.py": 7,  # formulas read whole columns by design
    "chemsuite/tools/reactions/enumerate.py": 6,
    "chemsuite/tools/reactions/common.py": 4,
    "chemsuite/tools/similarity/common.py": 4,
    "chemsuite/tools/statistics/panel.py": 4,
    "chemsuite/tools/medchem/dialog.py": 3,
    "chemsuite/tools/sar/matrix_panel.py": 3,
    "chemsuite/tools/similarity/search.py": 3,
}
#: Tools changing the table other than through ToolContext.apply_results.
DIRECT_WRITES_IN_TOOLS = {
    "chemsuite/tools/sar/free_wilson_ui.py": 1,
    "chemsuite/tools/reactions/enumerate.py": 1,  # fills a new table, not the document's
}


def test_tools_use_the_table_through_the_context() -> None:
    internals = _count(r"\.column\(|slots_of\(|\.blobs\(|\.mol\(|\.obj\(|\.displays\(", TOOLS)
    writes = _count(
        r"table\.(add_rows|set_cells|set_block|add_column|remove_rows|remove_columns|"
        r"convert_column|rename_column)\(",
        TOOLS,
    )
    over = _over(internals, COLUMN_INTERNALS_IN_TOOLS) + _over(writes, DIRECT_WRITES_IN_TOOLS)
    assert not over, (
        "Tools should read through ToolContext (structure_inputs, table.read) and write "
        "through ctx.apply_results:\n" + "\n".join(over)
    )


#: Whole-package counts of patterns MolManager drowned in (it had ~2,000 getattr/hasattr
#: fallbacks and 1,088 blind excepts). Each use should be rare and deliberate.
BUDGETS = {
    r"\b(?:getattr|hasattr)\(": 73,
    r"except\s+(?:Exception|BaseException)\b": 16,
    r"#\s*type:\s*ignore": 31,
}


def test_fallbacks_and_blind_excepts_do_not_grow() -> None:
    over = []
    for pattern, budget in BUDGETS.items():
        total = sum(_count(pattern).values())
        if total > budget:
            over.append(f"{pattern!r}: {total} (budget {budget})")
    assert not over, (
        "Prefer a typed attribute, a protocol or a narrow except over a fallback:\n"
        + "\n".join(over)
    )


#: Line budgets for the files already over MAX_LINES. MolManager's window grew to 60 mixins
#: one method at a time: split a file (as column_actions.py and list_actions.py were split out
#: of the window) rather than raise its budget.
MAX_LINES = 1000
LONG_FILES = {
    "chemsuite/app/main_window.py": 1796,
    "chemsuite/app/sketcher/canvas.py": 1341,
    "chemsuite/app/table_view.py": 1139,
    "chemsuite/app/viewer/panel.py": 1021,
}


def test_files_stay_within_their_line_budgets() -> None:
    lines = {name: text.count("\n") for name, text in _sources().items()}
    over = _over(lines, LONG_FILES, default=MAX_LINES)
    assert not over, "split these files (see LONG_FILES):\n" + "\n".join(over)
