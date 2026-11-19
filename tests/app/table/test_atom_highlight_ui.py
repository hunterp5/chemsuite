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

"""Highlight Atoms: a structure column drawn with the atoms another column names."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QSize

from chemsuite.app.table.highlights import AtomSource, resolve, set_source, sources
from tests.gui import gui

pytestmark = gui

SMILES = ["CC(=O)Nc1ccc(O)cc1", "Cc1ccccc1", ""]
#: as the Sites of Metabolism tool writes them: the three best atoms, best first
SITES = ["N3 0.81, C9 0.05, O10 0.04", "C0 0.62, C2 0.11, C4 0.09", ""]
TOP = ["N3", "C0", ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": [f"CMPD-{i}" for i in range(len(SMILES))],
            "SMILES": SMILES,
            "Sites (3A4)": SITES,
            "Top Site (3A4)": TOP,
            "Note": ["an amide", "a methyl", ""],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def submenu(menu, title: str):
    for action in menu.actions():
        if action.text() == title and action.menu() is not None:
            return action.menu()
    return None


def picture(chem, row: int = 0):
    """The drawing of one structure cell as it is on screen now."""
    column = chem.workspace.table.column("SMILES")
    delegate = chem.view.itemDelegate()
    chosen = delegate.highlights.get("SMILES")
    table = chem.workspace.table
    highlight = resolve(table, "SMILES", chosen, row) if isinstance(chosen, AtomSource) else None
    image = chem._images.image(column, row, QSize(200, 110), 1.0, False, highlight)
    assert image is not None
    return bytes(image.constBits())


def test_the_header_menu_offers_only_the_columns_that_name_atoms(chem) -> None:
    menu = chem.columns.build_menu("SMILES", [], 1)
    atoms = submenu(menu, "Highlight Atoms")
    assert atoms is not None
    assert [a.text() for a in atoms.actions() if not a.isSeparator()] == [
        "None",
        "Sites (3A4)",
        "Top Site (3A4)",
    ]  # "ID", "Note" and the structures themselves do not name atoms
    assert [a.isChecked() for a in atoms.actions() if not a.isSeparator()] == [True, False, False]
    menu.deleteLater()

    # a text column is not offered as a highlight for itself
    assert submenu(chem.columns.build_menu("Note", [], 4), "Highlight Atoms") is None


def test_choosing_a_column_draws_the_atoms_and_the_menu_remembers_it(chem) -> None:
    plain = picture(chem)
    menu = chem.columns.build_menu("SMILES", [], 1)
    atoms = submenu(menu, "Highlight Atoms")
    next(a for a in atoms.actions() if a.text() == "Sites (3A4)").trigger()
    menu.deleteLater()

    assert sources(chem.view) == {"SMILES": "Sites (3A4)"}
    marked = picture(chem)
    assert marked != plain
    # the ranked column and the top site alone are different pictures: three atoms against one
    set_source(chem.view, "SMILES", "Top Site (3A4)")
    assert picture(chem) not in (plain, marked)
    # a row the tool left blank has nothing to highlight
    delegate = chem.view.itemDelegate()
    chosen = delegate.highlights["SMILES"]
    assert resolve(chem.workspace.table, "SMILES", chosen, 2) is None

    again = chem.columns.build_menu("SMILES", [], 1)
    chosen = [a.text() for a in submenu(again, "Highlight Atoms").actions() if a.isChecked()]
    assert chosen == ["Top Site (3A4)"]
    again.deleteLater()


def test_clear_highlighting_stops_it(chem) -> None:
    plain = picture(chem)
    set_source(chem.view, "SMILES", "Sites (3A4)")
    assert picture(chem) != plain
    chem.clear_structure_highlights()
    assert not sources(chem.view)
    assert picture(chem) == plain


def test_a_column_written_for_other_structures_highlights_nothing(chem) -> None:
    table = chem.workspace.table
    plain = picture(chem)
    table.set_cells("Sites (3A4)", [int(table.oids[0])], ["Br3 0.81"])  # atom 3 is nitrogen
    set_source(chem.view, "SMILES", "Sites (3A4)")
    assert picture(chem) == plain


def test_the_highlight_is_saved_with_the_window_layout(chem) -> None:
    from chemsuite.app.document import dock_layout

    set_source(chem.view, "SMILES", "Sites (3A4)")
    saved = dock_layout.capture(chem)
    assert saved["atom_highlights"] == {"SMILES": "Sites (3A4)"}

    set_source(chem.view, "SMILES", None)
    dock_layout.restore(chem, saved)
    assert sources(chem.view) == {"SMILES": "Sites (3A4)"}

    dock_layout.restore(chem, {"atom_highlights": {"SMILES": "Gone"}})  # a column since deleted
    assert sources(chem.view) == {"SMILES": "Sites (3A4)"}
