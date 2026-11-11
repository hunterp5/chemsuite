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
"""The protein viewer's one source of truth: which proteins, how shown, what is selected.

The proteins are the document's ``protein`` assets; the session follows them (a protein added,
replaced, removed or undone shows at once). Each has a :class:`~chemsuite.protein.display.Display`;
the tree, the sequence strip, the selection bar and the 3D view all read and change the
session, never each other, and redraw on its signals:

- ``proteinsChanged``: proteins came or went, or one was shown or hidden;
- ``sceneChanged``: what is drawn changed (styles, cartoon, surface): geometry is rebuilt;
- ``lookChanged``: only colours or the selection changed: only the look is uploaded;
- ``hovered``: the atom under the pointer (the protein's id and atom index, or "" and -1);
- ``focusChanged``: a ligand was focused or the focus cleared (the panel moves the camera);
- ``measurementsChanged``: a measurement was added or removed, or atoms picked for one.

Selections are atom masks per protein. Saved ("named") selections, measurements and the
focused ligand are kept as atom uids, so they survive edits of the structure.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from PySide6.QtCore import QObject, Signal

from chemsuite.core.table import AssetsChanged, Table
from chemsuite.core.table.events import Changed
from chemsuite.molview.measure import describe
from chemsuite.protein import display as displays
from chemsuite.protein import scene as scenes
from chemsuite.protein import store
from chemsuite.protein.model import AMINO, LIGAND, NUCLEIC, Structure
from chemsuite.protein.select import Mask, SelectionError, select

KIND = "protein"
SET = "set"
ADD = "add"
TOGGLE = "toggle"
REMOVE = "remove"
#: Focus Ligand shows these as sticks and the rest faded
FOCUS_SITE = "byres polymer within 5 of $focus"
FOCUS_NEAR = "(water or ion) within 3.5 of $focus"
Point = tuple[str, int]  # a protein id and an atom uid


def structure_of(payload) -> Structure:
    """An asset's structure, decoding stored bytes."""
    return payload if isinstance(payload, Structure) else store.unpack(payload)


class ProteinSession(QObject):
    proteinsChanged = Signal()
    sceneChanged = Signal()
    lookChanged = Signal()
    selectionChanged = Signal()
    hovered = Signal(str, int)
    focusChanged = Signal()
    measurementsChanged = Signal()

    def __init__(self, table: Table, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self.displays: list[displays.Display] = []
        self.named: dict[str, dict[str, list[int]]] = {}
        self.focused: Point | None = None  # the focused ligand: one of its atoms
        self.measurements: list[tuple[Point, ...]] = []
        self.measuring = 0  # how many atoms the next measurement takes (0: not measuring)
        self.picked: list[Point] = []
        self._scene: displays.Scene | None = None
        self._pending: dict[str, dict] = {}  # saved display state for proteins not loaded yet
        self._unsubscribe = table.subscribe(self._on_table_event)
        self.reload()

    def close(self) -> None:
        self._unsubscribe()

    # proteins
    def reload(self) -> None:
        """Match the displays to the document's proteins, keeping how existing ones look."""
        assets = [a for a in self.table.assets.values() if a.kind == KIND]
        old = {d.id: d for d in self.displays}
        new: list[displays.Display] = []
        for asset in assets:
            structure = structure_of(asset.payload)
            current = old.get(asset.id)
            if current is not None and current.structure is structure:
                current.name = asset.name
                new.append(current)
                continue
            shown = displays.Display(asset.id, asset.name, structure)
            state = self._pending.pop(asset.id, None)
            if state is None and current is not None:
                state = current.state()  # replaced (an edit): keep the look if it fits
            if state is not None:
                shown.set_state(state)
            new.append(shown)
        self.displays = new
        self._scene = None
        self.proteinsChanged.emit()
        self.sceneChanged.emit()

    def _on_table_event(self, event) -> None:
        if isinstance(event, Changed) and any(isinstance(c, AssetsChanged) for c in event.changes):
            self.reload()

    def display(self, protein: str) -> displays.Display | None:
        return next((d for d in self.displays if d.id == protein), None)

    def set_visible(self, protein: str, visible: bool) -> None:
        shown = self.display(protein)
        if shown is not None and shown.visible != visible:
            shown.visible = visible
            self._scene = None
            self.proteinsChanged.emit()
            self.sceneChanged.emit()

    # the scene the view draws
    @property
    def scene(self) -> displays.Scene:
        if self._scene is None:
            self._scene = displays.compose(self.displays)
        return self._scene

    def _changed(self, rebuilt: bool) -> None:
        if rebuilt or self._scene is None:
            self._scene = None
            self.sceneChanged.emit()
        else:
            self._scene = displays.relook(self._scene)
            self.lookChanged.emit()

    # selection
    def selection(self) -> dict[str, Mask]:
        return {d.id: d.selected for d in self.displays}

    def selected_count(self) -> int:
        return int(sum(d.selected.sum() for d in self.displays if d.visible))

    def select(self, masks: Mapping[str, Mask], mode: str = SET) -> None:
        """Change the selection: ``masks`` per protein id (proteins left out count as
        nothing)."""
        for shown in self.displays:
            mask = masks.get(shown.id)
            mask = np.zeros(shown.structure.n_atoms, bool) if mask is None else mask
            if mode == SET:
                shown.selected = mask.copy()
            elif mode == ADD:
                shown.selected = shown.selected | mask
            elif mode == REMOVE:
                shown.selected = shown.selected & ~mask
            elif mode == TOGGLE:
                shown.selected = shown.selected ^ mask
        self._changed(rebuilt=False)
        self.selectionChanged.emit()

    def clear_selection(self) -> None:
        self.select({}, SET)

    def evaluate(self, text: str, protein: str | None = None) -> dict[str, Mask]:
        """What ``text`` selects in each visible protein, or in ``protein`` only (raises
        SelectionError)."""
        out = {}
        for shown in self.displays:
            if shown.visible and protein in (None, shown.id):
                named = self._named_masks(shown)
                out[shown.id] = select(shown.structure, text, named=named) & shown.active
        return out

    def select_text(self, text: str, mode: str = SET, protein: str | None = None) -> int:
        """Select what ``text`` selects (in ``protein`` only, if given); returns how many
        atoms. Raises SelectionError."""
        masks = self.evaluate(text, protein)
        self.select(masks, mode)
        return int(sum(m.sum() for m in masks.values()))

    def save_selection(self, name: str) -> None:
        if not name.isidentifier():
            raise SelectionError(f"'{name}' is not a name (letters, digits and _)", 0)
        self.named[name] = {
            d.id: d.structure.uid[d.selected].tolist() for d in self.displays if d.selected.any()
        }

    def _named_masks(self, shown: displays.Display) -> dict[str, Mask]:
        out = {}
        for name, by_protein in self.named.items():
            mask = np.zeros(shown.structure.n_atoms, bool)
            indices = shown.structure.indices(np.array(by_protein.get(shown.id, []), np.int64))
            mask[indices[indices >= 0]] = True
            out[name] = mask
        return out

    # changing how the selection is shown
    def _each_selected(self):
        for shown in self.displays:
            if shown.visible and shown.selected.any():
                yield shown, shown.selected

    def show_as(self, style: int) -> None:
        for shown, mask in self._each_selected():
            shown.show(mask, style)
        self._changed(rebuilt=True)

    def show_cartoon(self, on: bool) -> None:
        for shown, mask in self._each_selected():
            shown.show_cartoon(mask, on)
        self._changed(rebuilt=True)

    def show_surface(self, on: bool) -> None:
        for shown, mask in self._each_selected():
            shown.show_surface(mask, on)
        self._changed(rebuilt=True)

    def colour_by(self, scheme: str, everything: bool = False) -> None:
        for shown in self.displays:
            if everything:
                shown.colour_by(scheme)
            elif shown.visible and shown.selected.any():
                shown.colour_by(scheme, shown.selected)
        self._changed(rebuilt=False)

    def colour(self, rgb: tuple[float, float, float]) -> None:
        for shown, mask in self._each_selected():
            shown.colour(mask, rgb)
        self._changed(rebuilt=False)

    def choose_altloc(self, protein: str, residue: int, letter: str) -> None:
        shown = self.display(protein)
        if shown is not None:
            shown.choose_altloc(residue, letter)
            self._changed(rebuilt=True)

    # interactions
    @property
    def interactions_shown(self) -> bool:
        return any(d.query.any() for d in self.displays)

    def show_interactions(self, on: bool = True) -> None:
        """Show the interactions of the selection, or of every ligand when nothing is
        selected; ``on=False`` hides them."""
        selected = self.selected_count() > 0
        for shown in self.displays:
            if not on or not shown.visible:
                shown.show_interactions(None)
            elif selected:
                shown.show_interactions(shown.selected if shown.selected.any() else None)
            else:
                s = shown.structure
                shown.show_interactions(s.res_kind[s.residue_of] == LIGAND)
        self._changed(rebuilt=True)

    def interactions_of(self, protein: str, atom: int) -> list[str]:
        """The shown interactions this atom takes part in, in words."""
        from chemsuite.protein.interactions import describe as say

        shown = self.display(protein)
        if shown is None or shown not in self.scene.displays:
            return []
        found = self.scene.interactions[self.scene.displays.index(shown)]
        rows = np.flatnonzero((found.a == atom) | (found.b == atom))
        return [say(shown.structure, found, int(k)) for k in rows]

    # focus
    def ligands(self) -> list[tuple[str, int]]:
        """(protein id, residue) of every ligand in the visible proteins."""
        return [(d.id, int(r)) for d in self.displays if d.visible
                for r in np.flatnonzero(d.structure.res_kind == LIGAND)]  # fmt: skip

    def focus(self, protein: str, residue: int, restyle: bool = True) -> None:
        """Show this ligand in its pocket: the pocket's residues as sticks, the waters and ions
        touching it, its interactions; the rest of this protein faded, and its other waters,
        ions and small molecules hidden. ``restyle=False`` only
        fades (a saved focus coming back: the styles were saved already)."""
        shown = self.display(protein)
        if shown is None:
            return
        s = shown.structure
        ligand = (s.residue_of == residue) & shown.active
        named = {"focus": ligand}
        site = select(s, FOCUS_SITE, named=named) & shown.active
        waters = select(s, FOCUS_NEAR, named=named) & shown.active
        for other in self.displays:
            other.ghost = np.zeros(other.structure.n_atoms, bool)
            if restyle:
                other.show_interactions(None)
        if restyle:
            polymer = np.isin(s.res_kind[s.residue_of], (AMINO, NUCLEIC))
            shown.show(~(site | waters | ligand | polymer), scenes.HIDDEN)
            shown.show(site, scenes.STICK)
            shown.show(waters, scenes.BALL_AND_STICK)
            shown.show(ligand, scenes.BALL_AND_STICK)
            shown.show_interactions(ligand)
        shown.ghost = ~(site | waters | ligand)
        self.focused = (protein, int(s.uid[s.atoms_of(residue).start]))
        self._changed(rebuilt=True)
        self.focusChanged.emit()

    def focused_atoms(self) -> tuple[displays.Display, np.ndarray] | None:
        """The focused ligand's display and its pocket's atoms (ligand, site, waters)."""
        if self.focused is None:
            return None
        shown = self.display(self.focused[0])
        if shown is None or not shown.visible:
            return None
        return shown, ~shown.ghost

    def focused_residue(self) -> tuple[str, int] | None:
        if self.focused is None:
            return None
        shown = self.display(self.focused[0])
        index = -1 if shown is None else int(shown.structure.indices(
            np.array([self.focused[1]]))[0])  # fmt: skip
        if index < 0:
            return None
        return shown.id, int(shown.structure.residue_of[index])

    def clear_focus(self) -> None:
        if self.focused is None:
            return
        self.focused = None
        for shown in self.displays:
            shown.ghost = np.zeros(shown.structure.n_atoms, bool)
        self._changed(rebuilt=False)
        self.focusChanged.emit()

    # measurements
    def measure(self, atoms: int) -> None:
        """Measure between the next ``atoms`` atoms picked (2 distance, 3 angle, 4 dihedral;
        0 stops)."""
        self.measuring = atoms if atoms in (2, 3, 4) else 0
        self.picked = []
        self.measurementsChanged.emit()

    def pick(self, atom: int) -> bool:
        """Add scene ``atom`` to the measurement being made; False when not measuring."""
        if not self.measuring:
            return False
        found = self.scene.locate(atom)
        if found is None:
            return True
        shown, index = found
        point = (shown.id, int(shown.structure.uid[index]))
        if point in self.picked:
            self.picked.remove(point)
        else:
            self.picked.append(point)
        if len(self.picked) == self.measuring:
            self.measurements.append(tuple(self.picked))
            self.picked = []
        self.measurementsChanged.emit()
        return True

    def clear_measurements(self) -> None:
        self.measurements = []
        self.picked = []
        self.measurementsChanged.emit()

    def where(self, point: Point) -> tuple[int, np.ndarray] | None:
        """A point's scene atom and position, or None when it is not drawn."""
        shown = self.display(point[0])
        scene = self.scene
        if shown is None or shown not in scene.displays:
            return None
        index = int(shown.structure.indices(np.array([point[1]]))[0])
        if index < 0:
            return None
        return int(scene.global_atoms(shown, np.array([index]))[0]), shown.structure.xyz(0)[index]

    def measured(self) -> list[tuple[list[np.ndarray], str]]:
        """Each measurement whose atoms are all drawn: its positions and value in words."""
        out = []
        for points in self.measurements:
            places = [self.where(p) for p in points]
            if all(p is not None for p in places):
                xyz = [p[1] for p in places]
                out.append((xyz, describe(xyz)))
        return out

    # hover, from the view or the sequence
    def hover(self, atom: int) -> None:
        """``atom`` is an index into the scene (-1: nothing)."""
        found = self.scene.locate(atom)
        if found is None:
            self.hovered.emit("", -1)
        else:
            self.hovered.emit(found[0].id, found[1])

    # saving
    def state(self) -> dict:
        return {
            "proteins": {d.id: d.state() for d in self.displays},
            "named": self.named,
            "focused": list(self.focused) if self.focused else None,
            "measurements": [[list(p) for p in m] for m in self.measurements],
        }

    def set_state(self, state: Mapping) -> None:
        saved = dict(state.get("proteins", {}))
        self.named = {
            str(name): {str(p): [int(u) for u in uids] for p, uids in dict(proteins).items()}
            for name, proteins in dict(state.get("named", {})).items()
        }
        for shown in self.displays:
            if shown.id in saved:
                shown.set_state(saved.pop(shown.id))
        self._pending = saved
        self.measurements = _points(state.get("measurements"))
        focused = _points([[state.get("focused")]])
        self.focused = focused[0][0] if focused else None
        if self.focused is not None:
            place = self.focused_residue()
            if place is not None:
                self.focus(*place, restyle=False)
        self._scene = None
        self.proteinsChanged.emit()
        self.sceneChanged.emit()


def _points(data) -> list[tuple[Point, ...]]:
    """Measurements from saved state: lists of [protein id, uid] pairs; bad entries are
    dropped."""
    out = []
    for item in data if isinstance(data, list) else []:
        if not isinstance(item, list) or not 1 <= len(item) <= 4:
            continue
        points = [(str(p[0]), int(p[1])) for p in item
                  if isinstance(p, list) and len(p) == 2 and isinstance(p[1], int)]  # fmt: skip
        if len(points) == len(item):
            out.append(tuple(points))
    return out
