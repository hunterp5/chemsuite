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
"""The Protein Viewer's menu bar: the viewer is a window of its own, so everything it does is
here (and in the right-click menu), not in toolbar buttons.

- File: open, fetch, remove and save proteins; save the selection; close the window.
- Edit: copy the selection.
- Select: everything, kinds, chains, around the selection, saved selections.
- Display: Show, Hide, Colour, Label (on the selection); Hydrogens, Surface Transparency.
- Ligand: focus a ligand in its pocket, interactions, measurements, ligands to the table.
- View: the camera, which areas show, full screen.
- Window: minimise, the table's window.

Menus are filled when they open, so they always match the session.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtGui import QActionGroup, QColor, QKeySequence
from PySide6.QtWidgets import QFileDialog, QMenu, QMenuBar

from chemsuite.protein.structure import io as protein_io

from . import menus

if TYPE_CHECKING:
    from .panel import ProteinPanel

#: View ▸ Background: the theme's (light or dark), or a fixed colour for figures
THEME = "theme"
BACKGROUNDS = {
    THEME: ("Theme", None),
    "white": ("White", QColor(255, 255, 255)),
    "grey": ("Grey", QColor(128, 128, 128)),
    "black": ("Black", QColor(0, 0, 0)),
}

#: the commands a tool hands the panel (each opens its own dialog or job)
OPEN = "open"
FETCH = "fetch"
EXTRACT = "extract"
LIGAND = "ligand"  # add one ligand's row: (protein id, residue)
PROFILE = "profile"  # each ligand row's interactions, as columns
REMOVE = "remove"  # remove a protein: (protein id)
SHOW_TABLE = "show_table"  # bring the document's table window forward
Commands = Mapping[str, Callable[..., object]]
MENUS = ("File", "Edit", "Select", "Display", "Ligand", "View", "Window")


def build(panel: ProteinPanel, commands: Commands) -> QMenuBar:
    """The menu bar; its window shows it (it stays hidden while the panel has no window)."""
    bar = QMenuBar(panel)
    bar.hide()
    fills = (
        lambda m: _file(m, panel, commands),
        lambda m: _edit(m, panel),
        lambda m: _select(m, panel),
        lambda m: _display(m, panel),
        lambda m: _ligand(m, panel, commands),
        lambda m: _view(m, panel),
        lambda m: _window(m, panel, commands),
    )
    filled = []
    for title, fill in zip(MENUS, fills, strict=True):
        menu = bar.addMenu(title)
        menu.aboutToShow.connect(lambda m=menu, f=fill: (m.clear(), f(m)))
        filled.append((menu, fill))

    def refresh() -> None:
        """Fill the closed menus again, so their shortcuts are enabled as the session allows
        (a disabled action's shortcut does nothing)."""
        for menu, fill in filled:
            if not menu.isVisible():
                menu.clear()
                fill(menu)

    session = panel.session
    for signal in (session.selectionChanged, session.proteinsChanged, session.sceneChanged):
        signal.connect(refresh)
    refresh()
    return bar


def _command(menu: QMenu, title: str, commands: Commands, key: str, *args):
    action = menu.addAction(title, lambda: commands[key](*args) if key in commands else None)
    action.setEnabled(key in commands)
    return action


def _file(menu: QMenu, panel: ProteinPanel, commands: Commands) -> None:
    shown = panel.session.displays
    _command(menu, "Open Protein…", commands, OPEN).setShortcut(QKeySequence.StandardKey.Open)
    _command(menu, "Fetch from PDB…", commands, FETCH).setShortcut(QKeySequence("Ctrl+Shift+O"))
    menu.addSeparator()
    saving = menu.addMenu("Save Protein As")
    for display in shown:
        saving.addAction(f"{display.name}…", lambda p=display.id: save_protein(panel, p))
    saving.setEnabled(bool(shown))
    save = menu.addAction("Save Selection As…", lambda: save_selection(panel))
    save.setEnabled(panel.session.selected_count() > 0)
    export = menu.addAction("Export Image…", panel.export_image)
    export.setShortcut(QKeySequence("Ctrl+E"))
    export.setEnabled(panel.view is not None and bool(shown))
    menu.addSeparator()
    removing = menu.addMenu("Remove Protein")
    for display in shown:
        _command(removing, display.name, commands, REMOVE, display.id)
    removing.setEnabled(REMOVE in commands and bool(shown))
    menu.addSeparator()
    close = menu.addAction("Close Window", panel.close_viewer)
    close.setShortcut(QKeySequence.StandardKey.Close)


def _edit(menu: QMenu, panel: ProteinPanel) -> None:
    selected = panel.session.selected_count() > 0
    copy = menu.addAction("Copy", lambda: panel.copy_selection(protein_io.PDB))
    copy.setShortcut(QKeySequence.StandardKey.Copy)
    copy.setToolTip("Copy the selected atoms as PDB text")
    copy.setEnabled(selected)
    copying = menu.addMenu("Copy Selection As")
    copying.addAction("PDB", lambda: panel.copy_selection(protein_io.PDB))
    copying.addAction("mmCIF", lambda: panel.copy_selection(protein_io.MMCIF))
    copying.setEnabled(selected)
    image = menu.addAction("Copy Image", panel.copy_image)
    image.setShortcut(QKeySequence("Ctrl+Shift+C"))
    image.setEnabled(panel.view is not None and bool(panel.session.displays))


def _select(menu: QMenu, panel: ProteinPanel) -> None:
    session = panel.session
    everything = menu.addAction("All", lambda: panel.selection_bar.pick("all"))
    everything.setShortcut(QKeySequence.StandardKey.SelectAll)
    menu.addSeparator()
    panel.selection_bar.fill_menu(menu)
    menu.addSeparator()
    clear = menu.addAction("Clear Selection", session.clear_selection)
    clear.setEnabled(session.selected_count() > 0)


def _display(menu: QMenu, panel: ProteinPanel) -> None:
    session = panel.session
    menus.fill_show(menu.addMenu("Show"), session)
    menus.fill_hide(menu.addMenu("Hide"), session)
    menus.fill_colour(menu.addMenu("Colour"), session, panel)
    menus.fill_label(menu.addMenu("Label"), session)
    menu.addSeparator()
    menus.fill_hydrogens(menu.addMenu("Hydrogens"), session)
    menus.fill_transparency(menu.addMenu("Surface Transparency"), session)


def _ligand(menu: QMenu, panel: ProteinPanel, commands: Commands) -> None:
    session = panel.session
    menus.fill_focus(menu.addMenu("Focus"), session)
    skin = menu.addAction(
        "Pocket Surface", lambda: session.show_pocket_surface(not session.pocket_surface)
    )
    skin.setCheckable(True)
    skin.setChecked(session.pocket_surface)
    skin.setToolTip("Draw a see-through skin over the focused pocket, by hydrophobicity")
    clear = menu.addAction("Clear Focus", session.clear_focus)
    clear.setEnabled(session.focused is not None)
    menu.addSeparator()
    shown = menu.addAction(
        "Interactions", lambda: session.show_interactions(not session.interactions_shown)
    )
    shown.setCheckable(True)
    shown.setChecked(session.interactions_shown)
    shown.setToolTip(
        "Show the interactions of the selection, or of every ligand when nothing is selected"
    )
    menus.fill_measure(menu.addMenu("Measure"), session)
    menu.addSeparator()
    have = bool(session.displays)
    _command(menu, "Ligands to Table", commands, EXTRACT).setEnabled(EXTRACT in commands and have)
    _command(menu, "Interactions to Table", commands, PROFILE).setEnabled(
        PROFILE in commands and have
    )


def _view(menu: QMenu, panel: ProteinPanel) -> None:
    menu.addAction("Frame Selection", panel.frame_selection).setToolTip(
        "Fit the selection (everything, when nothing is selected) to the view (F)"
    )
    menu.addAction("Reset View", panel.reset_view)
    menu.addSeparator()
    background = menu.addMenu("Background")
    group = QActionGroup(background)
    for key, (title, _colour) in BACKGROUNDS.items():
        action = background.addAction(title, lambda k=key: panel.set_background(k))
        action.setCheckable(True)
        action.setChecked(panel.background == key)
        group.addAction(action)
    ortho = menu.addAction("Orthographic", lambda: panel.set_orthographic(not panel.orthographic))
    ortho.setCheckable(True)
    ortho.setChecked(panel.orthographic)
    ortho.setToolTip("Parallel projection: no perspective, sizes compare across depth")
    cue = menu.addAction("Depth Cue", lambda: panel.set_depth_cue(not panel.depth_cue))
    cue.setCheckable(True)
    cue.setChecked(panel.depth_cue)
    cue.setToolTip("Fade what lies further back into the background (fog)")
    for action in (ortho, cue):
        action.setEnabled(panel.view is not None)
    menu.addSeparator()
    for title, widget in (("Proteins List", panel.tree), ("Sequence", panel.sequence)):
        action = menu.addAction(title, lambda w=widget: w.setVisible(not w.isVisible()))
        action.setCheckable(True)
        action.setChecked(not widget.isHidden())
    menu.addSeparator()
    full = menu.addAction("Full Screen", panel.toggle_full_screen)
    full.setShortcut(QKeySequence.StandardKey.FullScreen)
    full.setCheckable(True)
    full.setChecked(panel.window().isFullScreen())


def _window(menu: QMenu, panel: ProteinPanel, commands: Commands) -> None:
    minimise = menu.addAction("Minimize", lambda: panel.window().showMinimized())
    minimise.setShortcut(QKeySequence("Ctrl+M"))
    menu.addSeparator()
    _command(menu, "Show Table", commands, SHOW_TABLE)


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
