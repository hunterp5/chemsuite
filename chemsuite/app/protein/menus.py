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
"""The protein viewer's menus: Show, Hide and Colour (on the selection), and the right-click
menu on an atom.

Every menu acts on the session's selection; the right-click menu first offers to select what
was clicked (its residue, chain or binding site), so the rest of it then applies there.
"""

from __future__ import annotations

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QColorDialog, QMenu, QWidget

from chemsuite.protein import colours, scene
from chemsuite.protein.select import SelectionError

from .session import ADD, SET, ProteinSession
from .tree import residue_expression

SHOW_STYLES = (
    ("Wire", scene.WIRE),
    ("Sticks", scene.STICK),
    ("Ball and Stick", scene.BALL_AND_STICK),
    ("Spacefill", scene.SPACEFILL),
)


def fill_show(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    for title, style in SHOW_STYLES:
        menu.addAction(title, lambda s=style: session.show_as(s))
    menu.addSeparator()
    menu.addAction("Cartoon", lambda: session.show_cartoon(True))
    menu.addAction("Surface", lambda: session.show_surface(True))
    _need_selection(menu, session)


def fill_hide(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    menu.addAction("Atoms", lambda: session.show_as(scene.HIDDEN))
    menu.addAction("Cartoon", lambda: session.show_cartoon(False))
    menu.addAction("Surface", lambda: session.show_surface(False))
    _need_selection(menu, session)


def fill_colour(menu: QMenu, session: ProteinSession, parent: QWidget) -> None:
    menu.clear()
    selected = session.selected_count() > 0
    for key, title in colours.SCHEMES.items():
        menu.addAction(f"By {title}", lambda k=key: session.colour_by(k, everything=not selected))
    menu.addSeparator()
    pick = menu.addAction("Colour…", lambda: _pick_colour(session, parent))
    pick.setEnabled(selected)


def _pick_colour(session: ProteinSession, parent: QWidget) -> None:
    colour = QColorDialog.getColor(QColor(46, 200, 74), parent, "Colour")
    if colour.isValid():
        session.colour((colour.redF(), colour.greenF(), colour.blueF()))


def _need_selection(menu: QMenu, session: ProteinSession) -> None:
    enabled = session.selected_count() > 0
    for action in menu.actions():
        action.setEnabled(enabled)


def context_menu(
    session: ProteinSession, atom: int, parent: QWidget, *, focus, copy_smiles=None
) -> QMenu:
    """The right-click menu for scene atom ``atom`` (-1: empty space)."""
    menu = QMenu(parent)
    found = session.scene.locate(atom)
    if found is not None:
        shown, index = found
        structure = shown.structure
        residue = int(structure.residue_of[index])
        title = menu.addAction(structure.atom_label(index))
        title.setEnabled(False)
        here = residue_expression(structure, residue)
        chain = f"chain {structure.chain_name[structure.chain_of[residue]]}"

        def pick(expression: str, mode: str = SET) -> None:
            try:
                session.select_text(expression, mode, protein=shown.id)
            except SelectionError:
                return

        menu.addAction("Select Residue", lambda: pick(here))
        menu.addAction("Add Residue to Selection", lambda: pick(here, ADD))
        menu.addAction("Select Chain", lambda: pick(chain))
        menu.addAction("Select Residues Within 5 Å", lambda: pick(f"byres within 5 of ({here})"))
        if copy_smiles is not None and structure.res_kind[residue] == 2:
            menu.addAction("Copy Ligand SMILES", lambda: copy_smiles(shown, residue))
        alternates = sorted(set(structure.altloc[structure.atoms_of(residue)].tolist()) - {""})
        if alternates:
            alt = menu.addMenu("Alternate Location")
            current = shown.altlocs.get(residue, structure.default_altlocs.get(residue, ""))
            for letter in alternates:
                action = alt.addAction(
                    letter, lambda x=letter: session.choose_altloc(shown.id, residue, x)
                )
                action.setCheckable(True)
                action.setChecked(letter == current)
        menu.addSeparator()
    show = menu.addMenu("Show")
    fill_show(show, session)
    hide = menu.addMenu("Hide")
    fill_hide(hide, session)
    colour = menu.addMenu("Colour")
    fill_colour(colour, session, parent)
    menu.addSeparator()
    framed = menu.addAction("Frame Selection", focus)
    framed.setEnabled(session.selected_count() > 0)
    clear = menu.addAction("Clear Selection", session.clear_selection)
    clear.setEnabled(session.selected_count() > 0)
    return menu
