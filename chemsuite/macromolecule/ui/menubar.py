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
"""The Macromolecule Viewer's menu bar: the viewer is a window of its own, so everything it
does is here (and in the right-click menus), not in toolbar buttons.

- File: open structures and sessions, fetch, save the shown structures, the session, the
  selection or an image; clear the session.
- Edit: undo and redo the document's steps; copy the selection; delete it, give its atoms an
  element or a charge, rename or renumber its residues or chain.
- Select: everything, kinds, chains, around the selection, saved selections.
- Display: Show, Hide, Color, Label (on the selection); Hydrogens; Surface (kind, color,
  transparency, carving, probe, quality, smoothing).
- Tools: what a tool of the app does to a structure (Prepare…), each in its own dialog;
  measurements; ligands and their interactions to the table.
- View: the camera, a ligand's pocket (focus, its surface), interactions, which areas show,
  the table's window, full screen.
- Sequence (pinned at the far right): the Sequence window, and which of its kind rows show.

Menus are filled when they open, so they always match the session.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtGui import QActionGroup, QColor, QKeySequence
from PySide6.QtWidgets import QFileDialog, QMenu, QMenuBar

from chemsuite.macromolecule.chem import interactions as contacts
from chemsuite.macromolecule.structure import io as protein_io
from chemsuite.macromolecule.structure.sequence import KIND_ROWS

from . import menus, surface_menu

if TYPE_CHECKING:
    from .panel import MacromoleculePanel

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
REMOVE = "remove"  # delete a structure: (structure id)
RENAME = "rename"  # rename a structure: (structure id)
DUPLICATE = "duplicate"  # copy a structure, look and all: (structure id)
SAVE_SESSION = "save_session"  # every structure and how it looks, as a session file
CLEAR = "clear"  # every structure removed: an empty workspace
PREPARE = "prepare"  # hydrogens, flips and protonation states (Tools ▸ Prepare…)
INSPECT = "inspect"  # every issue and choice, to tick and overrule (Tools ▸ Inspect…)
MINIMISE = "minimise"  # a prepared protein minimised, as a new one (Tools ▸ Minimise…)
DOCK = "dock"  # dock the table's ligands into a pocket (Tools ▸ Dock…)
RESCORE = "rescore"  # score poses that already exist (Tools ▸ Rescore Poses…)
SHOW_TABLE = "show_table"  # bring the document's table window forward
Commands = Mapping[str, Callable[..., object]]
MENUS = ("File", "Edit", "Select", "Display", "Tools", "View")
#: pinned at the far right of the bar (the last menu where the bar is the system's: macOS)
SEQUENCE = "Sequence"


class MenuBar(QMenuBar):
    """The viewer's menu bar, holding the Sequence menu apart from the others."""

    def __init__(self, parent) -> None:
        super().__init__(parent)
        self.sequence = QMenu(SEQUENCE, self)

    def menus(self) -> dict[str, QMenu]:
        """Every menu by its title, the Sequence menu included."""
        found = {a.text(): a.menu() for a in self.actions() if a.menu() is not None}
        return {**found, SEQUENCE: self.sequence}


def build(panel: MacromoleculePanel, commands: Commands) -> MenuBar:
    """The menu bar; its window shows it (it stays hidden while the panel has no window)."""
    bar = MenuBar(panel)
    bar.hide()
    fills = (
        lambda m: _file(m, panel, commands),
        lambda m: _edit(m, panel),
        lambda m: _select(m, panel),
        lambda m: _display(m, panel),
        lambda m: _tools(m, panel, commands),
        lambda m: _view(m, panel, commands),
    )
    filled = []
    for title, fill in zip(MENUS, fills, strict=True):
        menu = bar.addMenu(title)
        menu.aboutToShow.connect(lambda m=menu, f=fill: (m.clear(), f(m)))
        filled.append((menu, fill))
    _pin_sequence(bar)
    bar.sequence.aboutToShow.connect(lambda: (bar.sequence.clear(), _sequence(bar.sequence, panel)))
    filled.append((bar.sequence, lambda m: _sequence(m, panel)))

    def refresh() -> None:
        """Fill the closed menus again, so their shortcuts are enabled as the session allows
        (a disabled action's shortcut does nothing)."""
        for menu, fill in filled:
            if not menu.isVisible():
                menu.clear()
                fill(menu)

    # through the window's own SessionSignals, not the session's plain events: those outlive
    # the window, and a closed window's deleted menus must not be filled again
    signals = panel.signals
    for signal in (signals.selectionChanged, signals.proteinsChanged, signals.sceneChanged):
        signal.connect(refresh)
    refresh()
    return bar


def _pin_sequence(bar: MenuBar) -> None:
    """The Sequence menu at the far right, in the bar's corner, or the last menu where the
    system draws the bar (macOS, which has no corner)."""
    if bar.isNativeMenuBar():
        bar.addMenu(bar.sequence)
        return
    corner = QMenuBar(bar)  # a menu bar of its own: drawn as the titles beside it, any theme
    corner.setNativeMenuBar(False)
    corner.addMenu(bar.sequence)
    bar.setCornerWidget(corner, Qt.Corner.TopRightCorner)


def _command(menu: QMenu, title: str, commands: Commands, key: str, *args):
    action = menu.addAction(title, lambda: commands[key](*args) if key in commands else None)
    action.setEnabled(key in commands)
    return action


def _file(menu: QMenu, panel: MacromoleculePanel, commands: Commands) -> None:
    shown = panel.session.displays
    visible = any(d.visible for d in shown)
    _command(menu, "Open…", commands, OPEN).setShortcut(QKeySequence.StandardKey.Open)
    _command(menu, "Fetch from PDB…", commands, FETCH).setShortcut(QKeySequence("Ctrl+Shift+O"))
    menu.addSeparator()
    save = menu.addAction("Save…", lambda: save_structures(panel))
    save.setShortcut(QKeySequence.StandardKey.Save)
    save.setToolTip("Every shown structure in one PDB or mmCIF file")
    save.setEnabled(visible)
    session = _command(menu, "Save Session…", commands, SAVE_SESSION)
    session.setShortcut(QKeySequence("Ctrl+Shift+S"))
    session.setToolTip("Every structure, how it looks, selections, focus, labels, measurements, "
                       "reports and poses, in one file to open again")  # fmt: skip
    session.setEnabled(SAVE_SESSION in commands and bool(shown))
    selection = menu.addAction("Save Selection…", lambda: save_selection(panel))
    selection.setEnabled(panel.session.selected_count() > 0)
    export = menu.addAction("Export Image…", panel.export_image)
    export.setShortcut(QKeySequence("Ctrl+E"))
    export.setEnabled(panel.view is not None and bool(shown))
    menu.addSeparator()
    clear = _command(menu, "Clear Session…", commands, CLEAR)
    clear.setToolTip("Remove every structure: an empty workspace")
    clear.setEnabled(CLEAR in commands and bool(shown))


def _edit(menu: QMenu, panel: MacromoleculePanel) -> None:
    selected = panel.session.selected_count() > 0
    history = panel.session.store.history
    undo = menu.addAction(_step("Undo", history.undo_label), panel.undo)
    undo.setShortcut(QKeySequence.StandardKey.Undo)
    undo.setEnabled(history.can_undo)
    redo = menu.addAction(_step("Redo", history.redo_label), panel.redo)
    redo.setShortcut(QKeySequence.StandardKey.Redo)
    redo.setEnabled(history.can_redo)
    menu.addSeparator()
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
    menu.addSeparator()
    menus.fill_edit(menu, panel)


def _step(verb: str, label: str | None) -> str:
    return f"{verb} {label}" if label else verb


def _select(menu: QMenu, panel: MacromoleculePanel) -> None:
    session = panel.session
    everything = menu.addAction("All", lambda: panel.selection_bar.pick("all"))
    everything.setShortcut(QKeySequence.StandardKey.SelectAll)
    menu.addSeparator()
    panel.selection_bar.fill_menu(menu)
    menu.addSeparator()
    clear = menu.addAction("Clear Selection", session.clear_selection)
    clear.setEnabled(session.selected_count() > 0)


def _display(menu: QMenu, panel: MacromoleculePanel) -> None:
    session = panel.session
    menus.fill_show(menu.addMenu("Show"), session)
    menus.fill_hide(menu.addMenu("Hide"), session)
    menus.fill_colour(menu.addMenu("Color"), session, panel)
    menus.fill_label(menu.addMenu("Label"), session)
    menu.addSeparator()
    menus.fill_hydrogens(menu.addMenu("Hydrogens"), session)
    surface_menu.fill_surface(menu.addMenu("Surface"), session, panel)


def _view(menu: QMenu, panel: MacromoleculePanel, commands: Commands) -> None:
    session = panel.session
    menu.addAction("Frame Selection", panel.frame_selection).setToolTip(
        "Fit the selection (everything, when nothing is selected) to the view (F)"
    )
    menu.addAction("Reset View", panel.reset_view)
    menu.addSeparator()
    _ligand(menu.addMenu("Ligand"), panel)
    shown = menu.addAction(
        "Interactions", lambda: session.show_interactions(not session.interactions_shown)
    )
    shown.setCheckable(True)
    shown.setChecked(session.interactions_shown)
    shown.setToolTip(
        "Show the interactions of the selection, or of every ligand when nothing is selected"
    )
    kinds = menu.addMenu("Interaction Types")
    for kind, title in contacts.NAMES.items():
        action = kinds.addAction(
            title,
            lambda k=kind: session.show_interaction_kind(k, not session.interaction_kind_shown(k)),
        )
        action.setCheckable(True)
        action.setChecked(session.interaction_kind_shown(kind))
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
    tree = menu.addAction("Structures List", lambda: panel.tree.setVisible(panel.tree.isHidden()))
    tree.setCheckable(True)
    tree.setChecked(not panel.tree.isHidden())
    _command(menu, "Show Table", commands, SHOW_TABLE)
    menu.addSeparator()
    full = menu.addAction("Full Screen", panel.toggle_full_screen)
    full.setShortcut(QKeySequence.StandardKey.FullScreen)
    full.setCheckable(True)
    full.setChecked(panel.window().isFullScreen())


def _ligand(menu: QMenu, panel: MacromoleculePanel) -> None:
    """View ▸ Ligand: a ligand in its pocket."""
    session = panel.session
    menus.fill_focus(menu.addMenu("Focus"), session)
    skin = menu.addMenu("Pocket Surface")
    surface_menu.fill_pocket(skin, session)
    skin.menuAction().setToolTip(
        "A skin over the focused pocket, colored by hydrophobicity, of its protein's surface kind"
    )
    clear = menu.addAction("Clear Focus", session.clear_focus)
    clear.setEnabled(session.focused is not None)


def _tools(menu: QMenu, panel: MacromoleculePanel, commands: Commands) -> None:
    """What the app's tools do to a protein. Each opens a dialog of the tool's own, so the
    workbench only says the command's name."""
    prepare = _command(menu, "Prepare…", commands, PREPARE)
    prepare.setEnabled(PREPARE in commands and bool(panel.session.displays))
    prepare.setToolTip("Hydrogens, missing side chains, flips, His states and protonation")
    inspect = _command(menu, "Inspect…", commands, INSPECT)
    inspect.setEnabled(INSPECT in commands and bool(panel.session.displays))
    inspect.setToolTip("Everything a preparation would find and decide, to tick and overrule "
                       "before it is applied")  # fmt: skip
    relax = _command(menu, "Minimise…", commands, MINIMISE)
    relax.setEnabled(MINIMISE in commands and bool(panel.session.displays))
    relax.setToolTip("Relax a prepared protein in OpenMM, added as a new protein")
    menu.addSeparator()
    dock = _command(menu, "Dock…", commands, DOCK)
    dock.setEnabled(DOCK in commands and bool(panel.session.displays))
    dock.setToolTip("Dock the table's ligands into a pocket with AutoDock Vina; the poses "
                    "are listed beside the pocket")  # fmt: skip
    rescore = _command(menu, "Rescore Poses…", commands, RESCORE)
    rescore.setEnabled(RESCORE in commands and bool(panel.session.displays))
    rescore.setToolTip("Score poses that already exist where they stand")
    menu.addSeparator()
    menus.fill_measure(menu.addMenu("Measure"), panel.session)
    menu.addSeparator()
    have = bool(panel.session.displays)
    _command(menu, "Ligands to Table", commands, EXTRACT).setEnabled(EXTRACT in commands and have)
    _command(menu, "Interactions to Table", commands, PROFILE).setEnabled(
        PROFILE in commands and have
    )


def _sequence(menu: QMenu, panel: MacromoleculePanel) -> None:
    shown = panel.sequence_shown()
    viewer = menu.addAction("Sequence Viewer", lambda: panel.show_sequence(not shown))
    viewer.setShortcut(QKeySequence("Ctrl+Shift+L"))
    viewer.setCheckable(True)
    viewer.setChecked(shown)
    viewer.setToolTip("A window with each chain's sequence, and a row each for the ligands, "
                      "additives, ions and waters")  # fmt: skip
    rows = menu.addMenu("Rows")
    hidden = panel.sequence_hidden_kinds()
    for kind, title in KIND_ROWS.items():
        action = rows.addAction(title, lambda k=kind: panel.hide_sequence_kind(k, k not in hidden))
        action.setCheckable(True)
        action.setChecked(kind not in hidden)


def save_selection(panel: MacromoleculePanel, path: str | None = None) -> str:
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


def save_structures(panel: MacromoleculePanel, path: str | None = None) -> str:
    """Write every shown structure to one PDB or mmCIF file (by its suffix): one structure
    whole (every model and alternate location); several as one model, chains whose ID another
    took first renamed. Returns the path."""
    shown = [d for d in panel.session.displays if d.visible]
    if not shown:
        return ""
    if path is None:
        name = shown[0].name if len(shown) == 1 else "structures"
        path, _ = QFileDialog.getSaveFileName(
            panel, "Save", f"{name}.cif", "mmCIF (*.cif);;PDB (*.pdb)"
        )
    if not path:
        return ""
    text, renamed = protein_io.to_text_many(
        [(d.structure, None) for d in shown], protein_io.format_of(path)
    )
    _write(panel, path, text)
    if renamed:
        said = ", ".join(f"{shown[k].name} {old} → {new}" for (k, old), new in renamed.items())
        panel.status.setText(f"Saved with chains renamed: {said}")
    return path


def _write(panel: MacromoleculePanel, path: str, text: str) -> None:
    try:
        protein_io.write_text(path, text)
    except OSError as error:
        panel.status.setText(f"Could not save {Path(path).name}: {error.strerror or error}")
