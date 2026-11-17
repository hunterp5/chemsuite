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

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from tests.gui import gui

pytestmark = gui

DATA = Path(__file__).resolve().parents[2] / "data" / "protein"
SMILES = ["Cc1ccc(cc1)S(=O)(=O)Nc1ncccn1", "Cc1ccc(cc1)S(=O)(=O)Nc1ccccn1", "CCO"]


def molblock(smiles: str, count: int = 1, seed: int = 7) -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, count, randomSeed=seed)
    return Chem.MolToMolBlock(mol)


def _install(window, three_d: bool = True):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    cells = {"ID": ["A", "B", "C"], "SMILES": SMILES}
    if three_d:
        cells["Mol 3D"] = [molblock(s) for s in SMILES]
    table.add_rows(cells, undoable=False)
    table.convert_column("SMILES", "structure")
    if three_d:
        table.convert_column("Mol 3D", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _ctx(window):
    from chemsuite.app.document.tools import ToolContext

    return ToolContext(window)


def test_superposes_onto_the_current_row_in_one_undo_step(window, qtbot) -> None:
    from chemsuite.superpose.model import SHAPE, reference_id_of

    chem = _install(window)
    table = chem.workspace.table
    _ctx(chem).set_current_oid(int(table.oids[0]))
    dialog = chem.open_tool("superpose")
    assert dialog.structure.currentText() == "Mol 3D"  # the 3D column, not the 2D one
    assert dialog.reference.currentText() == "Current row"
    dialog.method.setCurrentIndex(dialog.method.findData(SHAPE))
    assert dialog.problem() == ""
    dialog.accept()
    qtbot.waitUntil(lambda: "Mol 3D superposed" in table, timeout=30_000)
    assert table.history.undo_label == "Superpose"
    scores = [table.cell("Shape Tanimoto", o) for o in table.oids]
    assert float(scores[0]) > 0.99 and all(scores)
    assert table.cell("Best Conformer", table.oids[0]) == "1"
    mol = table.column("Mol 3D superposed").mol(0)
    asset = table.assets[reference_id_of(mol)]
    assert asset.kind == "superpose-reference" and asset.name == "Row 1 (Mol 3D)"
    table.undo()
    assert "Mol 3D superposed" not in table and "Shape Tanimoto" not in table
    assert not table.assets


def test_the_reference_survives_save_rename_and_reopen(window, qtbot, tmp_path) -> None:
    from chemsuite.core.document import load_document, save_document
    from chemsuite.superpose.model import reference_id_of

    chem = _install(window)
    table = chem.workspace.table
    _ctx(chem).set_current_oid(int(table.oids[1]))
    chem.open_tool("superpose").accept()
    qtbot.waitUntil(lambda: "Mol 3D superposed" in table, timeout=30_000)
    table.rename_column("Mol 3D superposed", "Aligned")
    table.duplicate_column("Aligned", "Copy")
    save_document(table, tmp_path / "doc.chemsuite")
    reopened, _extras = load_document(tmp_path / "doc.chemsuite")
    for name in ("Aligned", "Copy"):
        found = reference_id_of(reopened.column(name).mol(0))
        assert reopened.assets[found].kind == "superpose-reference"


def test_problems_block_running(window) -> None:
    from chemsuite.superpose.model import CORE

    chem = _install(window)
    table = chem.workspace.table
    dialog = chem.open_tool("superpose")
    if _ctx(chem).current_oid() is None:
        assert "current row" in dialog.problem()
    _ctx(chem).set_current_oid(int(table.oids[0]))
    dialog.reference_column.setCurrentText("SMILES")
    assert dialog.problem() == "The reference row has no 3D structure in “SMILES”."
    dialog.reference_column.setCurrentText("Mol 3D")
    dialog.method.setCurrentIndex(dialog.method.findData(CORE))
    assert "core" in dialog.problem()
    dialog.core.setText("C1CCCCC1Br")
    assert dialog.problem() == "The core is not in the reference."
    dialog.core.setText("c1ccccc1")
    assert dialog.problem() == ""


def test_offers_the_ligands_of_the_document_s_proteins(window) -> None:
    from chemsuite.core.assets import Asset, PutAsset
    from chemsuite.macromolecule.structure import io

    chem = _install(window)
    structure = io.read(DATA / "1atp.cif.gz")
    asset = Asset("1atp", "protein", "1ATP", structure)
    chem.workspace.table.edit("Add protein", [PutAsset("1atp", None, asset)])
    dialog = chem.open_tool("superpose")
    titles = [dialog.reference.itemText(i) for i in range(dialog.reference.count())]
    assert "1ATP ▸ E/ATP 355" in titles
    dialog.reference.setCurrentText("1ATP ▸ E/ATP 355")
    assert not dialog.inputs.is_row_visible(dialog.reference_column)
    assert dialog.problem() == ""
    mol = dialog.reference_mol()
    assert mol.GetNumAtoms() == 31 and mol.GetConformer().Is3D()


def test_each_row_s_conformers_onto_its_first(window, qtbot) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.mols import mol_to_blob
    from chemsuite.core.table import Table
    from chemsuite.core.table.columns import StructureColumn

    mol = Chem.MolFromMolBlock(molblock("CCCCCC(=O)NCCc1ccccc1"), removeHs=False)
    AllChem.EmbedMultipleConfs(mol, 4, randomSeed=5)
    table = Table()
    column = StructureColumn.from_parsed([mol_to_blob(mol)], ["CCCCCC(=O)NCCc1ccccc1"])
    table.add_rows({"Mol": column}, undoable=False)
    window._install(Workspace(table))
    dialog = window.open_tool("superpose")
    dialog.mode.setCurrentIndex(dialog.mode.findData("conformers"))
    assert not dialog.method_section.isVisibleTo(dialog)
    assert dialog.problem() == ""
    dialog.accept()
    qtbot.waitUntil(lambda: "Mol superposed" in table, timeout=30_000)
    assert float(table.cell("Max Conformer RMSD (Å)", table.oids[0])) >= 0
