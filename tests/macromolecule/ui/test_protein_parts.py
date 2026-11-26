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

"""The proteins list's parts: each included in calculations or left out, renamed, duplicated
or deleted on its own."""

from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import Qt

from chemsuite.core.table import Table
from chemsuite.macromolecule.structure import parts
from chemsuite.macromolecule.structure.model import LIGAND, WATER
from chemsuite.macromolecule.ui.tree import DELETE, DUPLICATE, RENAME

from .test_protein_ui import add_protein, structure


@pytest.fixture
def table() -> Table:
    t = Table()
    add_protein(t, "1ATP", "1atp.cif.gz")
    return t


@pytest.fixture
def panel(qtbot, table, monkeypatch):
    from chemsuite.app.shell import graphics
    from chemsuite.macromolecule.ui.panel import MacromoleculePanel

    monkeypatch.setattr(graphics, "_problem", None)
    widget = MacromoleculePanel(table)
    qtbot.addWidget(widget)
    widget.resize(900, 600)
    return widget


def _row(panel, text: str):
    top = panel.tree.topLevelItem(0)
    return next(top.child(i) for i in range(top.childCount()) if top.child(i).text(0) == text)


def _menu(panel, row) -> dict[str, bool]:
    """The commands *row*'s menu offers, and whether each is enabled."""
    return {a.text(): a.isEnabled() for a in panel.tree.menu_for(row).actions()}


def test_unticking_a_part_leaves_it_out_and_an_undo_brings_it_back(panel, table, qtbot) -> None:
    session = panel.session
    shown = session.display("1ATP")
    waters = shown.structure.atom_kind == WATER
    row = _row(panel, "Waters (103)")
    assert row.checkState(0) == Qt.CheckState.Checked
    row.setCheckState(0, Qt.CheckState.Unchecked)
    qtbot.waitUntil(lambda: session.display("1ATP").excluded.any())
    shown = session.display("1ATP")
    assert np.array_equal(shown.excluded, waters)
    assert not shown.active[waters].any()  # neither drawn nor selectable
    session.select_text("water", protein="1ATP")
    assert not shown.selected.any()
    asset = table.assets["1ATP"]
    assert asset.payload is shown.structure  # the structure itself is unchanged
    assert np.array_equal(parts.excluded(asset.meta, shown.structure), waters)
    assert _row(panel, "Waters (103)").checkState(0) == Qt.CheckState.Unchecked
    table.undo()
    assert not session.display("1ATP").excluded.any()
    assert _row(panel, "Waters (103)").checkState(0) == Qt.CheckState.Checked


def test_a_parts_menu_acts_on_the_part(panel) -> None:
    assert _menu(panel, _row(panel, "ATP E355")) == {
        "Rename…": True, "Duplicate": True, "Delete": True}  # fmt: skip
    assert _menu(panel, _row(panel, "Waters (103)"))["Rename…"] is False
    top = panel.tree.topLevelItem(0)
    assert set(_menu(panel, top)) == {"Rename…", "Duplicate", "Delete…"}  # the protein's


def test_duplicating_a_ligand_adds_one_copy_in_its_chain(panel, table) -> None:
    s = structure("1atp.cif.gz")
    panel.part_command(DUPLICATE, "1ATP", "chain E and resi 355 and resn ATP", "ATP E355",
                       "residue")  # fmt: skip
    made = panel.session.display("1ATP").structure
    assert made.n_residues == s.n_residues + 1 and made.n_chains == s.n_chains
    atp = np.flatnonzero((made.res_name == "ATP") & (made.res_kind == LIGAND))
    assert len(atp) == 2
    last = int(s.res_seq[s.residues_of(0)].max())  # chain E's last number
    assert sorted(int(made.res_seq[r]) for r in atp) == [355, last + 1]
    labels = [panel.tree.topLevelItem(0).child(i).text(0)
              for i in range(panel.tree.topLevelItem(0).childCount())]  # fmt: skip
    assert sum(label.startswith("ATP E") for label in labels) == 2
    assert table.assets.__len__() == 1  # no new protein
    table.undo()
    assert panel.session.display("1ATP").structure.n_residues == s.n_residues


def test_duplicating_a_chain_adds_a_chain_and_waters_copy_the_waters(panel) -> None:
    s = structure("1atp.cif.gz")
    panel.part_command(DUPLICATE, "1ATP", "chain I and polymer", "chain I", "chain")
    made = panel.session.display("1ATP").structure
    assert made.n_chains == s.n_chains + 1
    new = made.residues_of(made.n_chains - 1)
    polymer = [r for r in range(*s.residues_of(1).indices(s.n_residues)) if s.res_kind[r] == 0]
    assert made.res_name[new].tolist() == s.res_name[polymer].tolist()  # its waters stay
    panel.part_command(DUPLICATE, "1ATP", "water", "waters", "group")
    again = panel.session.display("1ATP").structure
    assert (again.res_kind == WATER).sum() == 2 * (s.res_kind == WATER).sum()


def test_deleting_and_renaming_a_part(panel, table) -> None:
    s = structure("1atp.cif.gz")
    panel.part_command(DELETE, "1ATP", "ion", "ions", "group")
    made = panel.session.display("1ATP").structure
    assert made.n_atoms == s.n_atoms - 2 and "Deleted ions" in panel.status.text()
    dialog = panel.part_command(RENAME, "1ATP", "chain I and polymer", "chain I", "chain")
    dialog.field.setText("P")
    dialog.accept()
    assert "P" in panel.session.display("1ATP").structure.chain_name.tolist()
    dialog = panel.part_command(RENAME, "1ATP", "chain E and resi 355 and resn ATP", "ATP E355",
                                "residue")  # fmt: skip
    dialog.field.setText("ANP")
    dialog.accept()
    assert "ANP" in panel.session.display("1ATP").structure.res_name.tolist()


def test_a_reopened_document_keeps_its_parts_left_out(qtbot, monkeypatch) -> None:
    """A reopened document holds stored bytes and its exclusions in the asset's meta: they
    show at once, and including a part keeps the display (and its selection) as it was."""
    from chemsuite.app.shell import graphics
    from chemsuite.core.table import Asset, PutAsset
    from chemsuite.macromolecule.ui.panel import MacromoleculePanel

    monkeypatch.setattr(graphics, "_problem", None)
    s = structure("1atp.cif.gz")
    waters = s.atom_kind == WATER
    meta = parts.with_excluded({}, s, waters)
    table = Table()
    table.edit("Open", [PutAsset("1ATP", None, Asset("1ATP", "protein", "1ATP", s.to_bytes(),
                                                     meta))])  # fmt: skip
    widget = MacromoleculePanel(table)
    qtbot.addWidget(widget)
    session = widget.session
    shown = session.display("1ATP")
    assert np.array_equal(shown.excluded, waters)
    assert _row(widget, "Waters (103)").checkState(0) == Qt.CheckState.Unchecked
    session.select_text("resn ATP")
    _row(widget, "Waters (103)").setCheckState(0, Qt.CheckState.Checked)
    qtbot.waitUntil(lambda: not session.display("1ATP").excluded.any())
    assert session.display("1ATP") is shown and shown.selected.any()
