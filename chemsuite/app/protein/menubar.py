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
"""The Protein Viewer's own menu bar: it is a workspace of its own, so its files, editing,
display and pocket work live here rather than in the main window's menus.

File: Open Protein…, Fetch from PDB…, Save Selection As…, Close. Edit: copying, selecting.
View: Show, Hide, Colour, Label, Hydrogens, Surface Transparency, the camera, which areas
show.
Pocket: interactions, focus, measurements, ligands to the table. Menus are filled when they
open, so they always match the session.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import QFileDialog, QMenu, QMenuBar

from chemsuite.protein import io as protein_io

from . import menus, toolbar

if TYPE_CHECKING:
    from .panel import ProteinPanel

#: the commands a tool hands the panel (each opens its own dialog or job)
OPEN = "open"
FETCH = "fetch"
EXTRACT = "extract"
LIGAND = "ligand"  # add one ligand's row: (protein id, residue)
PROFILE = "profile"  # each ligand row's interactions, as columns
REMOVE = "remove"  # remove a protein: (protein id)
Commands = Mapping[str, Callable[..., object]]


def build(panel: ProteinPanel, commands: Commands) -> QMenuBar:
    bar = QMenuBar(panel)
    bar.setNativeMenuBar(False)  # in the viewer's own window, not the app's menu bar
    _file(bar.addMenu("File"), panel, commands)
    _edit(bar.addMenu("Edit"), panel)
    _view(bar.addMenu("View"), panel)
    _pocket(bar.addMenu("Pocket"), panel, commands)
    return bar


def _file(menu: QMenu, panel: ProteinPanel, commands: Commands) -> None:
    def fill() -> None:
        menu.clear()
        opening = menu.addAction("Open Protein…", commands.get(OPEN, lambda: None))
        opening.setShortcut(QKeySequence.StandardKey.Open)
        opening.setEnabled(OPEN in commands)
        menu.addAction("Fetch from PDB…", commands.get(FETCH, lambda: None)).setEnabled(
            FETCH in commands
        )
        removing = menu.addMenu("Remove Protein")
        for shown in panel.session.displays:
            removing.addAction(shown.name, lambda p=shown.id: commands[REMOVE](p))
        removing.setEnabled(REMOVE in commands and bool(panel.session.displays))
        menu.addSeparator()
        saving = menu.addMenu("Save Protein As")
        for shown in panel.session.displays:
            saving.addAction(f"{shown.name}…", lambda p=shown.id: save_protein(panel, p))
        saving.setEnabled(bool(panel.session.displays))
        save = menu.addAction("Save Selection As…", lambda: save_selection(panel))
        save.setEnabled(panel.session.selected_count() > 0)
        menu.addSeparator()
        menu.addAction("Close Viewer", panel.close_viewer)

    menu.aboutToShow.connect(fill)
    fill()


def _edit(menu: QMenu, panel: ProteinPanel) -> None:
    session = panel.session

    def fill() -> None:
        menu.clear()
        selected = session.selected_count() > 0
        copying = menu.addMenu("Copy Selection As")
        copying.addAction("PDB", lambda: panel.copy_selection(protein_io.PDB))
        copying.addAction("mmCIF", lambda: panel.copy_selection(protein_io.MMCIF))
        copying.setEnabled(selected)
        menu.addSeparator()
        menu.addAction("Select All", lambda: session.select_text("all"))
        clear = menu.addAction("Clear Selection", session.clear_selection)
        clear.setEnabled(selected)
        save = menu.addAction("Save Selection…", panel.selection_bar.save)
        save.setEnabled(selected)

    menu.aboutToShow.connect(fill)
    fill()


def _view(menu: QMenu, panel: ProteinPanel) -> None:
    session = panel.session

    def fill() -> None:
        menu.clear()
        menus.fill_show(menu.addMenu("Show"), session)
        menus.fill_hide(menu.addMenu("Hide"), session)
        menus.fill_colour(menu.addMenu("Colour"), session, panel)
        menus.fill_label(menu.addMenu("Label"), session)
        menus.fill_hydrogens(menu.addMenu("Hydrogens"), session)
        menus.fill_transparency(menu.addMenu("Surface Transparency"), session)
        menu.addSeparator()
        menu.addAction("Frame Selection", panel.frame_selection)
        menu.addAction("Reset View", panel.reset_view)
        menu.addSeparator()
        for title, widget in (("Proteins List", panel.tree), ("Sequence", panel.sequence)):
            action = menu.addAction(title, lambda w=widget: w.setVisible(not w.isVisible()))
            action.setCheckable(True)
            action.setChecked(widget.isVisibleTo(panel))

    menu.aboutToShow.connect(fill)
    fill()


def _pocket(menu: QMenu, panel: ProteinPanel, commands: Commands) -> None:
    session = panel.session

    def fill() -> None:
        menu.clear()
        shown = menu.addAction("Interactions", lambda: session.show_interactions(
            not session.interactions_shown))  # fmt: skip
        shown.setCheckable(True)
        shown.setChecked(session.interactions_shown)
        toolbar.fill_focus(menu.addMenu("Focus"), session)
        toolbar.fill_measure(menu.addMenu("Measure"), session)
        menu.addSeparator()
        extract = menu.addAction("Extract Ligands to Table", commands.get(EXTRACT, lambda: None))
        extract.setEnabled(EXTRACT in commands and bool(session.displays))
        profile = menu.addAction("Interactions to Table", commands.get(PROFILE, lambda: None))
        profile.setEnabled(PROFILE in commands and bool(session.displays))

    menu.aboutToShow.connect(fill)
    fill()


def save_selection(panel: ProteinPanel, path: str | None = None) -> str:
    """Write the selected atoms to a PDB or mmCIF file (by its suffix); returns the path."""
    if path is None:
        path, _ = QFileDialog.getSaveFileName(
            panel, "Save Selection", "selection.pdb", "PDB (*.pdb);;mmCIF (*.cif)"
        )
    if not path:
        return ""
    text = panel.selection_text(protein_io.format_of(path))
    _write(panel, path, text)
    return path


def save_protein(panel: ProteinPanel, protein: str, path: str | None = None) -> str:
    """Write a whole protein (every model and alternate location) to a PDB or mmCIF file."""
    shown = panel.session.display(protein)
    if shown is None:
        return ""
    if path is None:
        path, _ = QFileDialog.getSaveFileName(
            panel, "Save Protein", f"{shown.name}.cif", "mmCIF (*.cif);;PDB (*.pdb)"
        )
    if not path:
        return ""
    _write(panel, path, protein_io.to_text(shown.structure, protein_io.format_of(path)))
    return path


def _write(panel: ProteinPanel, path: str, text: str) -> None:
    try:
        protein_io.write_text(path, text)
    except OSError as error:
        panel.status.setText(f"Could not save {Path(path).name}: {error.strerror or error}")
