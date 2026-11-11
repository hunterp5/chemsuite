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
"""The Protein Viewer tool: open, fetch (a stub opener: no network), ligand rows, the table link."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from PySide6.QtWidgets import QPushButton

from chemsuite.app.tools import ToolContext
from chemsuite.core.settings import POLICY_ENV
from chemsuite.tools.protein import tool

from .conftest import gui

pytestmark = gui
DATA = Path(__file__).resolve().parents[1] / "data" / "protein"


@pytest.fixture(autouse=True)
def network(tmp_path, monkeypatch):
    """Downloads allowed, but the CCD is a stub that knows only 0G6 (and records each ask)."""
    import urllib.error

    from chemsuite.protein import components

    monkeypatch.setenv(POLICY_ENV, str(tmp_path / "no-policy.json"))
    asked: list[str] = []

    def ccd(request, timeout=None):
        import gzip
        import io

        asked.append(request.full_url)
        if request.full_url.endswith("/0G6.cif"):
            return io.BytesIO(gzip.decompress((DATA / "0g6_ccd.cif.gz").read_bytes()))
        raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(components, "fetch_bytes", _with_opener(components.fetch_bytes, ccd))
    return asked


@pytest.fixture
def opened(window, qtbot, monkeypatch):
    """1ATP opened from a file into the small ID/pIC50 table."""
    from chemsuite.app import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    job = tool.open_protein(ToolContext(window), [str(DATA / "1atp.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    return window


def test_only_the_viewer_is_in_the_main_menus(window) -> None:
    assert window.tool_action("protein").isEnabled()
    assert [t.id for t in tool.TOOLS] == ["protein"]  # files open from the viewer itself


def test_the_viewer_has_a_menu_of_its_own(window) -> None:
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    protein = window.menuBar().actions()[titles.index("Protein")].menu()
    assert window.tool_action("protein") in protein.actions()
    chemistry = window.menuBar().actions()[titles.index("Chemistry")].menu()
    assert window.tool_action("protein") not in chemistry.actions()


def test_the_viewer_opens_as_its_own_window(window, qtbot, monkeypatch) -> None:
    from chemsuite.app import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    window.resize(1400, 900)
    window.show()
    qtbot.waitExposed(window)
    panel = tool.open_viewer(ToolContext(window))
    dock = next(d for d in window._panels if d.widget() is panel)
    assert dock.isFloating()
    assert dock.width() >= window.width() * 0.8
    assert panel.views.currentWidget() is panel.empty  # no proteins yet: ways to add one
    titles = [b.text() for b in panel.empty.findChildren(QPushButton)]
    assert titles == ["Open Protein…", "Fetch from PDB…"]


def test_the_viewer_menus_open_and_fetch(window, monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(tool, "open_protein", lambda ctx: calls.append("open"))
    monkeypatch.setattr(tool, "fetch_from_pdb", lambda ctx: calls.append("fetch"))
    panel = tool.open_viewer(ToolContext(window))
    file_menu = panel.menu_bar.actions()[0].menu()
    file_menu.aboutToShow.emit()
    for action in file_menu.actions():
        if action.text() in ("Open Protein…", "Fetch from PDB…"):
            action.trigger()
    assert calls == ["open", "fetch"]
    assert [a.text() for a in panel.menu_bar.actions()] == ["File", "Edit", "View", "Pocket"]


def test_opening_adds_the_protein_and_its_ligands_in_one_step(opened) -> None:
    table = opened.workspace.table
    assert list(table.assets) == ["1ATP"]
    oid = int(table.oids[-1])
    assert table.cell("Protein", oid) == "1ATP" and table.cell("Ligand ID", oid) == "E/ATP 355"
    assert table.kind("Ligand") == "structure" and table.cell("Ligand", oid).startswith("Nc1")
    assert table.cell("Bond Orders", oid) == "File"  # the file's own chemical component tables
    panel = ToolContext(opened).find_panel("protein")
    assert panel is not None and [d.id for d in panel.session.displays] == ["1ATP"]
    table.undo()
    assert not table.assets and len(table.oids) == 3 and "Ligand" not in table


def test_a_ligand_the_file_does_not_define_comes_from_the_ccd_once(window, qtbot, network):
    ctx = ToolContext(window)
    for _ in range(2):
        job = tool.open_protein(ctx, [str(DATA / "1ppb.cif.gz")])
        qtbot.waitUntil(lambda job=job: not job.running, timeout=30_000)
    table = window.workspace.table
    oid = int(table.oids[-1])
    assert table.cell("Ligand ID", oid) == "H/0G6 1"
    assert table.cell("Bond Orders", oid) == "CCD (downloaded)"
    assert network == ["https://files.rcsb.org/ligands/download/0G6.cif"]  # then the cache
    assert "guessed" not in table.assets["1PPB"].meta


def test_a_guessed_ligand_is_marked_in_the_tree(window, qtbot, tmp_path, monkeypatch) -> None:
    from chemsuite.app import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"network_allowed": False}))
    monkeypatch.setenv(POLICY_ENV, str(policy))
    job = tool.open_protein(ToolContext(window), [str(DATA / "1ppb.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    table = window.workspace.table
    assert table.cell("Bond Orders", int(table.oids[-1])) == "Guessed from geometry"
    assert table.assets["1PPB"].meta["guessed"] == "H/0G6 1"
    tree = ToolContext(window).find_panel("protein").tree
    top = tree.topLevelItem(0)
    item = next(top.child(i) for i in range(top.childCount())
                if top.child(i).text(0).startswith("0G6"))  # fmt: skip
    assert item.text(0).endswith("⚠") and "guessed" in item.toolTip(0)


def test_a_second_copy_gets_its_own_id(opened, qtbot) -> None:
    job = tool.open_protein(ToolContext(opened), [str(DATA / "1atp.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    assert list(opened.workspace.table.assets) == ["1ATP", "1ATP-2"]


def test_an_unreadable_file_says_why(window, qtbot, tmp_path) -> None:
    bad = tmp_path / "broken.pdb"
    bad.write_text("not a protein\n")
    job = tool.open_protein(ToolContext(window), [str(bad)])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    assert "Could not read broken.pdb" in window.statusBar().currentMessage()
    assert not window.workspace.table.assets


def test_the_current_row_focuses_its_ligand(opened) -> None:
    table = opened.workspace.table
    panel = ToolContext(opened).find_panel("protein")
    ToolContext(opened).set_current_oid(int(table.oids[0]))  # a row with no protein
    assert panel.session.focused is None
    ToolContext(opened).set_current_oid(int(table.oids[-1]))
    shown = panel.session.display("1ATP")
    assert panel.session.focused_residue() == ("1ATP", 338) and shown.ghost.any()


def test_extract_ligands_fills_in_proteins_without_rows(opened) -> None:
    ctx = ToolContext(opened)
    assert tool.extract_ligands(ctx) == 0  # 1ATP has its row already
    table = opened.workspace.table
    table.remove_rows([int(table.oids[-1])])
    assert tool.extract_ligands(ctx) == 1


def test_a_ligand_can_be_added_once_from_the_viewer(opened) -> None:
    ctx = ToolContext(opened)
    table = opened.workspace.table
    panel = ctx.find_panel("protein")
    shown = panel.session.display("1ATP")
    atp = int(shown.structure.res_seq.tolist().index(355))
    assert tool.add_ligand(ctx, "1ATP", atp) == 0  # it has its row already
    assert "in the table already" in opened.statusBar().currentMessage()
    table.remove_rows([int(table.oids[-1])])
    rows = len(table.oids)
    panel.commands["ligand"]("1ATP", atp)
    assert len(table.oids) == rows + 1
    assert table.cell("Ligand ID", int(table.oids[-1])) == "E/ATP 355"
    assert tool.add_ligand(ctx, "NONE", 0) == 0


def test_fetch_downloads_by_id(window, qtbot, monkeypatch) -> None:
    from chemsuite.protein import fetch as pdb
    from chemsuite.tools.protein.dialogs import FetchDialog

    asked: list[str] = []

    def opener(request, timeout=None):
        import io

        asked.append(request.full_url)
        return io.BytesIO((DATA / "1crn.cif.gz").read_bytes())

    monkeypatch.setattr(pdb, "fetch_bytes", _with_opener(pdb.fetch_bytes, opener))
    dialog = FetchDialog(ToolContext(window))
    qtbot.addWidget(dialog)
    assert not dialog.primary.isEnabled()  # no ID yet
    dialog.code.setText("12")
    assert "not a PDB ID" in dialog.status.text()
    dialog.code.setText("1crn")
    assert dialog.primary.text() == "Fetch 1CRN"
    dialog.accept()
    qtbot.waitUntil(lambda: not dialog.job.running, timeout=30_000)
    assert asked == ["https://files.rcsb.org/download/1CRN.cif.gz"]
    assert list(window.workspace.table.assets) == ["1CRN"]


def test_fetch_says_when_downloads_are_off(window, qtbot, tmp_path, monkeypatch) -> None:
    from chemsuite.tools.protein.dialogs import FetchDialog

    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"network_allowed": False}))
    monkeypatch.setenv(POLICY_ENV, str(policy))
    dialog = FetchDialog(ToolContext(window))
    qtbot.addWidget(dialog)
    dialog.code.setText("1ATP")
    assert "turned off" in dialog.status.text() and not dialog.primary.isEnabled()


def test_copy_selection_as_pdb(opened) -> None:
    from PySide6.QtGui import QGuiApplication

    panel = ToolContext(opened).find_panel("protein")
    panel.session.select_text("ligand")
    text = panel.copy_selection("pdb")
    assert text.count("HETATM") == 31 and QGuiApplication.clipboard().text() == text
    assert panel.copy_smiles(panel.session.display("1ATP"), 338).startswith("Nc1")


def test_the_viewer_reopens_from_its_saved_state(opened, qtbot) -> None:
    panel = ToolContext(opened).find_panel("protein")
    panel.session.select_text("chain I")
    panel.session.save_selection("peptide")
    restore = next(t for t in tool.TOOLS if t.id == "protein").restore
    again = restore(ToolContext(opened), panel.state())
    assert again is not panel and again.session.named == panel.session.named


def _with_opener(fetch_bytes, opener):
    def call(url, **kwargs):
        kwargs["opener"] = opener
        return fetch_bytes(url, **kwargs)

    return call


def test_a_second_protein_adds_to_the_ligand_column(opened, qtbot) -> None:
    job = tool.open_protein(ToolContext(opened), [str(DATA / "3htb.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    table = opened.workspace.table
    oid = int(table.oids[-1])
    assert table.cell("Protein", oid) == "3HTB" and table.cell("Ligand", oid).count("\n") == 0


def test_save_selection_writes_pdb_or_mmcif(opened, tmp_path) -> None:
    from chemsuite.app.protein.menubar import save_selection
    from chemsuite.protein import io

    panel = ToolContext(opened).find_panel("protein")
    panel.session.select_text("ligand")
    for name in ("atp.pdb", "atp.cif"):
        path = save_selection(panel, str(tmp_path / name))
        back = io.read(path)
        assert back.n_atoms == 31 and str(back.res_name[0]) == "ATP"
