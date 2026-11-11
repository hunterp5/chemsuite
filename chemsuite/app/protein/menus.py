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
"""The protein viewer's menus: Show, Hide, Colour and Label (on the selection), Hydrogens and
Transparency (for every protein), and the right-click menu on an atom.

Every menu acts on the session's selection; the right-click menu first offers to select what
was clicked (its residue, chain or binding site), so the rest of it then applies there.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtGui import QActionGroup, QColor
from PySide6.QtWidgets import QColorDialog, QMenu, QWidget

from chemsuite.molview.geometry import HYDROGENS
from chemsuite.protein import colours, scene
from chemsuite.protein.display import TRANSPARENCIES
from chemsuite.protein.io import MMCIF, PDB
from chemsuite.protein.model import LIGAND
from chemsuite.protein.select import SelectionError

from .session import ADD, ATOM_LABEL, RESIDUE_LABEL, SET, ProteinSession
from .toolbar import MEASURES
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


def fill_label(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    selected = session.selected_count() > 0
    menu.addAction("Residues", lambda: session.label(RESIDUE_LABEL)).setEnabled(selected)
    menu.addAction("Atoms", lambda: session.label(ATOM_LABEL)).setEnabled(selected)
    menu.addSeparator()
    menu.addAction("Clear Labels", session.clear_labels).setEnabled(bool(session.labels))


def fill_hydrogens(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    group = QActionGroup(menu)
    for key, title in HYDROGENS.items():
        action = menu.addAction(title, lambda k=key: session.show_hydrogens(k))
        action.setCheckable(True)
        action.setChecked(session.hydrogens == key)
        group.addAction(action)


def fill_transparency(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    group = QActionGroup(menu)
    for value in TRANSPARENCIES:
        title = f"{value:.0%}" if value else "Opaque"
        action = menu.addAction(title, lambda v=value: session.set_transparency(v))
        action.setCheckable(True)
        action.setChecked(abs(session.transparency - value) < 1e-6)
        group.addAction(action)


def _label_here(menu: QMenu, session: ProteinSession, protein: str, structure, atom: int) -> None:
    mask = np.zeros(structure.n_atoms, bool)
    mask[atom] = True
    for title, kind in (("Residue", RESIDUE_LABEL), ("Atom", ATOM_LABEL)):
        menu.addAction(title, lambda k=kind: session.label(k, {protein: mask}))


def _measure_from(session: ProteinSession, atoms: int, atom: int) -> None:
    """Start a measurement of ``atoms`` atoms with this one picked."""
    session.measure(atoms)
    session.pick(atom)


def _pick_colour(session: ProteinSession, parent: QWidget) -> None:
    colour = QColorDialog.getColor(QColor(46, 200, 74), parent, "Colour")
    if colour.isValid():
        session.colour((colour.redF(), colour.greenF(), colour.blueF()))


def _need_selection(menu: QMenu, session: ProteinSession) -> None:
    enabled = session.selected_count() > 0
    for action in menu.actions():
        action.setEnabled(enabled)


def context_menu(
    session: ProteinSession,
    atom: int,
    parent: QWidget,
    *,
    focus,
    copy_smiles=None,
    copy_as=None,
    add_ligand=None,
) -> QMenu:
    """The right-click menu for scene atom ``atom`` (-1: empty space)."""
    menu = QMenu(parent)
    found = session.scene.locate(atom)
    if found is not None and found[0] is session.pose:
        found = None  # a pose's atoms: the menu for empty space
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
        if structure.res_kind[residue] == LIGAND:
            menu.addAction("Focus This Ligand", lambda: session.focus(shown.id, residue))
        if copy_smiles is not None and structure.res_kind[residue] == LIGAND:
            menu.addAction("Copy Ligand SMILES", lambda: copy_smiles(shown, residue))
        if add_ligand is not None and structure.res_kind[residue] == LIGAND:
            menu.addAction("Add Ligand to Table", lambda: add_ligand(shown.id, residue))
        _label_here(menu.addMenu("Label"), session, shown.id, structure, index)
        measuring = menu.addMenu("Measure From Here")
        for name, atoms in MEASURES:
            measuring.addAction(name, lambda n=atoms: _measure_from(session, n, atom))
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
    if session.interactions_shown:
        menu.addAction("Hide Interactions", lambda: session.show_interactions(False))
    else:
        menu.addAction("Show Interactions", session.show_interactions)
    if session.focused is not None:
        menu.addAction("Clear Focus", session.clear_focus)
    menu.addSeparator()
    framed = menu.addAction("Frame Selection", focus)
    framed.setEnabled(session.selected_count() > 0)
    clear = menu.addAction("Clear Selection", session.clear_selection)
    clear.setEnabled(session.selected_count() > 0)
    if copy_as is not None:
        copying = menu.addMenu("Copy Selection As")
        copying.addAction("PDB", lambda: copy_as(PDB))
        copying.addAction("mmCIF", lambda: copy_as(MMCIF))
        copying.setEnabled(session.selected_count() > 0)
    return menu
