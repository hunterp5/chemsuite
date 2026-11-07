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

from __future__ import annotations

from rdkit import Chem

from .conftest import gui

pytestmark = gui


def _chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["A", "B", "C"], "SMILES": ["CCCCO", "c1ccccc1", "nonsense"]}, undoable=False
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_generate_dialog_end_to_end(window, qtbot) -> None:
    chem = _chem(window)
    assert chem.tool_action("conformers").isEnabled()
    dialog = chem.open_tool("conformers")
    assert dialog.structure.currentText() == "SMILES"
    dialog.count.setValue(3)
    dialog.accept()
    table = chem.workspace.table
    qtbot.waitUntil(lambda: "SMILES 3D" in table, timeout=30_000)
    assert table.kind("SMILES 3D") == "structure"
    assert [table.cell("SMILES 3D", o) for o in table.oids] == ["CCCCO", "c1ccccc1", ""]
    counts = [table.cell("Conformers", o) for o in table.oids]
    assert counts[1] == "1" and int(counts[0]) >= 1 and counts[2] == ""
    assert "Conformer Energy Range (kcal/mol)" in table
    mol = table.column("SMILES 3D").mol(table.slots_of(table.oids[:1])[0])
    assert mol.GetConformer().Is3D()
    # Running again keeps the earlier results and adds new columns.
    dialog = chem.open_tool("conformers")
    assert dialog.count.value() == 3  # remembered
    dialog.accept()
    qtbot.waitUntil(lambda: table.history.undo_label == "Generate conformers", timeout=30_000)
    qtbot.wait(50)
    assert "SMILES 3D (2)" in table and "Conformers (2)" in table
    assert table.cell("SMILES 3D", int(table.oids[0])) == "CCCCO"


def test_window_needs_optimisation(window) -> None:
    chem = _chem(window)
    dialog = chem.open_tool("conformers")
    dialog.force_field.setCurrentIndex(dialog.force_field.findData(None))
    assert not dialog.use_window.isEnabled()
    assert dialog.settings().energy_window is None and dialog.settings().optimize is None
    dialog.reject()


def test_export_conformers(window, qtbot, tmp_path) -> None:
    from chemsuite.tools.conformers.dialog import output_name

    chem = _chem(window)
    dialog = chem.open_tool("conformers")
    dialog.count.setValue(4)
    dialog.accept()
    table = chem.workspace.table
    qtbot.waitUntil(lambda: output_name("SMILES") in table, timeout=30_000)
    export = chem.open_tool("export_conformers")
    assert export.structure.currentText() == "SMILES 3D"
    path = tmp_path / "confs.sdf"
    job = export.export(path)
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    records = list(Chem.SDMolSupplier(str(path), removeHs=False))
    first = int(table.cell("Conformers", int(table.oids[0])))
    assert len(records) == first + 2  # benzene has one conformer, "nonsense" an empty record
    assert records[-1].GetNumAtoms() == 0
    assert records[0].GetProp("ID") == "A"
    export.reject()
