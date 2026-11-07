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
"""A saved document reopens as it was left (``chemsuite.app.session``)."""

from __future__ import annotations

import ast
import importlib
import pkgutil
from pathlib import Path

import pytest

from .conftest import gui

pytestmark = gui

TOOLS = Path(__file__).parents[2] / "chemsuite" / "tools"


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": [f"C{i}" for i in range(200)],
            "pIC50": [str(5 + (i % 30) / 10) for i in range(200)],
            "Series": ["AB"[i % 2] for i in range(200)],
        },
        undoable=False,
    )
    table.convert_column("pIC50", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    window.resize(1000, 600)
    window.show()
    return window


def reopen(window, qtbot, path: Path) -> None:
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    assert window._save_to(path)
    window._install(Workspace(Table()))
    window.open_path(path)
    qtbot.waitUntil(lambda: window.workspace.path == path, timeout=10_000)
    qtbot.wait(20)


@pytest.mark.parametrize("open_when_saved", [False, True])
def test_filters_pane_reopens_as_it_was(chem, qtbot, tmp_path, open_when_saved) -> None:
    from chemsuite.app.dock_layout import is_open

    chem.filters.add_filter("Series").set_checked({"A"})
    chem.filters.apply_now()
    chem.filter_dock.setVisible(open_when_saved)
    reopen(chem, qtbot, tmp_path / "filters.csd")
    assert chem.filters.cards  # the filter is kept either way
    assert len(chem.workspace.table.visible_oids()) == 100
    assert is_open(chem.filter_dock) == open_when_saved


def test_selection_current_cell_and_scroll_are_kept(chem, qtbot, tmp_path) -> None:
    table, view, model = chem.workspace.table, chem.view, chem.workspace.model
    view.setCurrentIndex(model.index(150, model.column_index("pIC50")))
    view.scrollTo(model.index(120, 0), view.ScrollHint.PositionAtTop)
    chem.workspace.selection.set([int(table.oids[3]), int(table.oids[150])], "test")
    top = view.indexAt(view.viewport().rect().topLeft()).row()
    reopen(chem, qtbot, tmp_path / "view.csd")
    table, view, model = chem.workspace.table, chem.view, chem.workspace.model
    assert chem.workspace.selection.oids == {int(table.oids[3]), int(table.oids[150])}
    current = chem.view.currentIndex()
    assert (current.row(), model.column_name(current.column())) == (150, "pIC50")
    assert view.indexAt(view.viewport().rect().topLeft()).row() == top
    assert not chem.isWindowModified()


def test_damaged_session_parts_are_skipped(chem) -> None:
    from chemsuite.app import session

    session.restore(chem, {part.key: "junk" for part in session.PARTS})
    session.restore(chem, {"view": {"selected": ["x", 10**9], "current": [1], "top": "a"}})
    assert set(session.capture(chem)) <= session.KEYS


def _panel_classes() -> dict[str, str]:
    """Panel classes in chemsuite.tools: tool_id -> class name."""
    found = {}
    for info in pkgutil.walk_packages([str(TOOLS)], "chemsuite.tools."):
        module = importlib.import_module(info.name)
        for name, value in vars(module).items():
            tool_id = getattr(value, "tool_id", None)
            if (
                isinstance(value, type)
                and isinstance(tool_id, str)
                and value.__module__ == info.name
            ):
                found[tool_id] = f"{info.name}.{name}"
    return found


def test_every_tool_panel_is_saved_with_the_document() -> None:
    """A tool that opens a panel must save it: the panel class has ``tool_id`` and
    ``state()``, and the tool with that id has ``restore``. New panels are caught here."""
    from chemsuite.app.main_window import registered_tools

    tools = {tool.id: tool for tool in registered_tools()}
    panels = _panel_classes()
    for tool_id, where in panels.items():
        assert tool_id in tools and tools[tool_id].restore is not None, (
            f"{where} is never reopened: give tool {tool_id!r} a restore"
        )
    # every module that adds a panel defines (or uses) a panel class
    for path in TOOLS.rglob("*.py"):
        tree = ast.parse(path.read_text())
        adds = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "add_panel"
        ]
        if adds:
            text = path.read_text()
            assert "tool_id" in text or "restore=" in text, (
                f"{path} adds a panel that is never saved"
            )
