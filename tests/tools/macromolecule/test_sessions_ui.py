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
# ruff: noqa: F811  (the fixtures below are test_protein_tool's, imported)
"""The workspace as a whole: Save Session, Open (Replace or Add), Clear Session; and a
structure's Rename and Duplicate. Each change to the structures is one undo step."""

from __future__ import annotations

from chemsuite.app.document.tools import ToolContext
from chemsuite.macromolecule.workspace import file as session_file
from chemsuite.macromolecule.workspace.prepared import PreparedRecord
from chemsuite.tools.macromolecule import session_files, tool
from tests.gui import gui

from .test_protein_tool import DATA, network, opened, wait, with_rows  # noqa: F401  (fixtures)

pytestmark = gui


def panel_of(window):
    return ToolContext(window).find_panel(tool.TOOL_ID)


def saved(window, qtbot, tmp_path) -> str:
    path = str(tmp_path / "kinase.cspx")
    wait(qtbot, session_files.save_session(ToolContext(window), path))
    return path


def test_a_saved_session_opens_again_in_place_of_what_is_there(opened, qtbot, tmp_path) -> None:
    ctx = ToolContext(opened)
    session = panel_of(opened).session
    session.select_text("chain E and resi 10-20")
    session.save_selection("site")
    session.show_surface(True)
    looks = session.display("1ATP").state()
    path = saved(opened, qtbot, tmp_path)
    assert session_file.read(path).assets[0].id == "1ATP"
    session_files.clear_session(ctx, tool.KEEP_ROWS)
    assert not session_files.structures(ctx) and not session.named
    wait(qtbot, session_files.open_session(ctx, path))  # nothing open: no question
    table = opened.workspace.table
    assert list(session_files.structures(ctx)) == ["1ATP"]
    assert table.history.undo_label == "Open kinase.cspx"
    assert session.named.keys() == {"site"}
    assert session.display("1ATP").state() == looks


def test_a_session_added_beside_the_same_structure_gets_new_ids(opened, qtbot, tmp_path) -> None:
    ctx = ToolContext(opened)
    session = panel_of(opened).session
    session.select_text("chain E and resi 10-20")
    session.save_selection("site")
    session.measurements.append((("1ATP", 1), ("1ATP", 2)))
    path = saved(opened, qtbot, tmp_path)
    wait(qtbot, session_files.open_session(ctx, path, session_files.ADD))
    assert list(session_files.structures(ctx)) == ["1ATP", "1ATP-2"]
    assert session.named.keys() == {"site", "site_2"}
    assert session.named["site_2"].keys() == {"1ATP-2"}
    assert (("1ATP-2", 1), ("1ATP-2", 2)) in session.measurements
    table = opened.workspace.table
    table.undo()  # the added structures go in one step
    assert list(session_files.structures(ctx)) == ["1ATP"]


def test_replace_asks_first_when_structures_are_open(opened, qtbot, tmp_path, monkeypatch) -> None:
    ctx = ToolContext(opened)
    path = saved(opened, qtbot, tmp_path)
    monkeypatch.setattr(session_files, "_ask_how", lambda ctx, name: None)
    wait(qtbot, session_files.open_session(ctx, path))
    assert list(session_files.structures(ctx)) == ["1ATP"]  # cancelled: as it was
    monkeypatch.setattr(session_files, "_ask_how", lambda ctx, name: session_files.REPLACE)
    wait(qtbot, session_files.open_session(ctx, path))
    assert list(session_files.structures(ctx)) == ["1ATP"]
    table = opened.workspace.table
    assert table.history.undo_label == "Open kinase.cspx"


def test_a_file_that_is_not_a_session_says_so(opened, qtbot, tmp_path) -> None:
    junk = tmp_path / "junk.cspx"
    junk.write_bytes(b"junk" * 100)
    ctx = ToolContext(opened)
    wait(qtbot, session_files.open_session(ctx, str(junk)))
    assert list(session_files.structures(ctx)) == ["1ATP"]


def test_open_sends_sessions_and_structure_files_their_ways(opened, monkeypatch) -> None:
    calls = []
    monkeypatch.setattr(session_files, "open_session", lambda ctx, p: calls.append(("s", p)))
    monkeypatch.setattr(session_files, "open_protein", lambda ctx, ps: calls.append(("f", ps)))
    session_files.open_files(ToolContext(opened), ["a.cspx", "b.cif", "C.CSPX"])
    assert calls == [("s", "a.cspx"), ("s", "C.CSPX"), ("f", ["b.cif"])]


def test_clear_session_asks_once_about_rows_and_undoes_in_one_step(with_rows, monkeypatch) -> None:
    ctx = ToolContext(with_rows)
    table = with_rows.workspace.table
    rows = len(table.oids)
    asked = []
    monkeypatch.setattr(
        session_files, "_ask_rows", lambda ctx, n, r: asked.append((n, r)) or tool.DELETE_ROWS
    )
    session = panel_of(with_rows).session
    session.labels[("1ATP", 1, "atom")] = None
    session.select_text("chain E")
    session.show_surface(True)
    looks = session.display("1ATP").state()
    assert session_files.clear_session(ctx)
    assert asked == [(1, 1)] and not session_files.structures(ctx)
    assert len(table.oids) == rows - 1 and not session.labels and not session.displays
    assert table.history.undo_label == "Clear session"
    table.undo()
    assert list(session_files.structures(ctx)) == ["1ATP"] and len(table.oids) == rows
    assert session.display("1ATP").state() == looks  # back as it looked


def test_rename_keeps_the_id_and_undoes(opened) -> None:
    ctx = ToolContext(opened)
    assert session_files.rename_structure(ctx, "1ATP", "  kinase ")
    assert opened.workspace.table.assets["1ATP"].name == "kinase"
    tree = panel_of(opened).tree
    assert tree.topLevelItem(0).text(0) == "kinase"
    assert not session_files.rename_structure(ctx, "1ATP", " ")  # a name is needed
    opened.workspace.table.undo()
    assert opened.workspace.table.assets["1ATP"].name == "1ATP"


def test_duplicate_copies_the_look_and_the_record(opened) -> None:
    ctx = ToolContext(opened)
    workspace = tool.macromolecule_workspace(ctx)
    session = workspace.session()
    session.select_text("chain E")
    session.show_surface(True)
    record = PreparedRecord("1ATP", "digest", 7.0)
    workspace.set_prepared(record)
    new = session_files.duplicate_structure(ctx, "1ATP")
    assert new == "1ATP-copy" and list(session_files.structures(ctx)) == ["1ATP", new]
    original, copy = session.display("1ATP"), session.display(new)
    assert copy.state() == original.state()
    copied = workspace.prepared(new)
    assert copied is not None and copied.protein == new and copied.ph == 7.0
    assert opened.workspace.table.history.undo_label == "Duplicate 1ATP"


def test_the_lists_menu_offers_rename_duplicate_and_delete(opened, monkeypatch) -> None:
    tree = panel_of(opened).tree
    assert tree._offered == {"Rename…": True, "Duplicate": True, "Delete…": True}
    heard = []
    monkeypatch.setattr(session_files, "rename_structure", lambda ctx, id: heard.append(id))
    monkeypatch.setattr(session_files, "duplicate_structure", lambda ctx, id: heard.append(id))
    tree.renameRequested.emit("1ATP")
    tree.duplicateRequested.emit("1ATP")
    assert heard == ["1ATP", "1ATP"]
