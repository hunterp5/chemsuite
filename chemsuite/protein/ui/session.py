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
replaced, removed or undone shows at once). Each has a
:class:`~chemsuite.protein.look.display.Display`; the tree, the sequence strip, the selection
bar and the 3D view all read and change the session, never each other, and redraw on its signals:

- ``proteinsChanged``: proteins came or went, or one was shown or hidden;
- ``sceneChanged``: what is drawn changed (styles, cartoon, surface, the model shown):
  geometry is rebuilt;
- ``lookChanged``: only colours or the selection changed: only the look is uploaded;
- ``hovered``: the atom under the pointer (the protein's id and atom index, or "" and -1);
- ``focusChanged``: a ligand was focused or the focus cleared (the panel moves the camera);
- ``measurementsChanged``: a measurement was added or removed, or atoms picked for one;
- ``labelsChanged``: labels were added or cleared;
- ``edited``: a protein's structure was replaced (an edit here, or its undo): the panel follows
  the current row again, so a pose is placed in the new structure;
- ``rowsEdited``: table rows had a cell of ``row_columns`` changed (the host's
  :class:`~.host.LinkedRow` columns: the panel follows the current row again, so a ligand
  moved in the sketcher shows as a pose).

With a way to run work in the background (``run``, as a tool's ``ctx.run`` in the quick
lane), slow parts (surfaces, interactions, large cartoons) are built off the GUI thread by a
:class:`~.builds.Builder`; the scene shows the last one until the new one arrives.

Selections are atom masks per protein. Saved ("named") selections, measurements, labels and
the focused ligand are kept as atom uids, so they survive edits of the structure.

Edits (delete, element, charge, rename, renumber: :mod:`chemsuite.protein.structure.edit`)
replace the protein's asset, each one step of the document's history; the display of the
edited structure takes the old one's look atom by atom (:meth:`Display.carry`). Renaming a
ligand that table rows name renames it there too, in the same step, through ``relabel`` (the
host's :meth:`~.host.LinkedRow.relabel`).

A pose (a row's own 3D ligand that is not one of the protein's residues: docked, or moved in
the sketcher) is drawn by a display of its own, ``pose``, holding the protein and the pose
merged (:func:`chemsuite.protein.chem.poses.merged`) with only the pose drawn, so its interactions
with the protein are found as for any ligand. It is not a document asset and is not saved: the
table's current row brings it back.
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Mapping, Sequence
from dataclasses import replace

import numpy as np
from PySide6.QtCore import QObject, Signal

from chemsuite.core.table import AssetsChanged, PutAsset, Table
from chemsuite.core.table.events import CellsChanged, Changed
from chemsuite.molview.measure import describe
from chemsuite.protein.chem import poses
from chemsuite.protein.look import display as displays
from chemsuite.protein.look import scene as scenes
from chemsuite.protein.structure import edit as edits
from chemsuite.protein.structure import store
from chemsuite.protein.structure.model import AMINO, LIGAND, NUCLEIC, Structure
from chemsuite.protein.structure.select import Mask, SelectionError, select

from .builds import Builder, Run

KIND = "protein"
SET = "set"
ADD = "add"
TOGGLE = "toggle"
REMOVE = "remove"
#: Focus Ligand shows these as sticks and the rest faded
#: whole residues: ``byres`` takes what follows it, so the distance test is in brackets
FOCUS_SITE = "byres (polymer within 5 of $focus)"
FOCUS_NEAR = "(water or ion) within 3.5 of $focus"
Point = tuple[str, int]  # a protein id and an atom uid
#: in expressions, ``$selection`` is what is selected now
SELECTION = "selection"
ATOM_LABEL = "atom"
RESIDUE_LABEL = "residue"
#: a pose's display id: this, then its protein's
POSE = "pose:"


structure_of = store.as_structure  # an asset's structure, decoding stored bytes


class ProteinSession(QObject):
    proteinsChanged = Signal()
    sceneChanged = Signal()
    lookChanged = Signal()
    selectionChanged = Signal()
    hovered = Signal(str, int)
    focusChanged = Signal()
    measurementsChanged = Signal()
    labelsChanged = Signal()
    rowsEdited = Signal(object)  # oids whose cells in row_columns changed
    edited = Signal()  # a protein's structure was replaced

    def __init__(
        self,
        table: Table,
        parent: QObject | None = None,
        run: Run | None = None,
        row_columns: Collection[str] = (),
        relabel: Callable[[str, Mapping[str, str]], Sequence[object]] | None = None,
        count_named: Callable[[str, Collection[str]], int] | None = None,
    ) -> None:
        super().__init__(parent)
        self.table = table
        self._relabel = relabel
        self._count_named = count_named
        self.row_columns = frozenset(row_columns)
        self.builder = Builder(self, run) if run is not None else None
        self.displays: list[displays.Display] = []
        self.named: dict[str, dict[str, list[int]]] = {}
        self.focused: Point | None = None  # the focused ligand: one of its atoms
        self.pose: displays.Display | None = None  # a row's own pose, in its protein
        self.pose_row: int | None = None  # that row
        self._posed_in = None  # the protein's structure the pose was placed in
        self._site: np.ndarray | None = None  # the focused pocket's residues (for its skin)
        self.measurements: list[tuple[Point, ...]] = []
        self.measuring = 0  # how many atoms the next measurement takes (0: not measuring)
        self.pocket_surface = True  # Focus draws a skin over the pocket
        #: what is labelled, in order: (protein, atom uid, ATOM_LABEL or RESIDUE_LABEL); a
        #: residue's label is kept on its first atom
        self.labels: dict[tuple[str, int, str], None] = {}
        self.picked: list[Point] = []
        self._scene: displays.Scene | None = None
        self._pending: dict[str, dict] = {}  # saved display state for proteins not loaded yet
        self._shown_assets: dict[str, object] = {}  # the asset each display was made from
        self._unsubscribe = table.subscribe(self._on_table_event)
        self.reload()

    def close(self) -> None:
        self._unsubscribe()
        if self.builder is not None:
            self.builder.cancel_all()

    # proteins
    def reload(self) -> None:
        """Match the displays to the document's proteins, keeping how existing ones look."""
        assets = [a for a in self.table.assets.values() if a.kind == KIND]
        old = {d.id: d for d in self.displays}
        new: list[displays.Display] = []
        replaced = False
        for asset in assets:
            current = old.get(asset.id)
            same = self._shown_assets.get(asset.id) is asset or (
                current is not None and current.structure is asset.payload
            )  # assets are immutable: the same one (or a renamed one) needs no new display
            if current is not None and same:
                current.name = asset.name
                new.append(current)
                continue
            structure = structure_of(asset.payload)
            shown = displays.Display(
                asset.id, asset.name, structure, deferred=self.builder is not None
            )
            state = self._pending.pop(asset.id, None)
            if state is not None:
                shown.set_state(state)
            elif current is not None:
                shown.carry(current)  # replaced (an edit, or its undo): atoms keep their look
                replaced = True
            new.append(shown)
        self.displays = new
        self._shown_assets = {asset.id: asset for asset in assets}
        if self.focused is not None and self.display(self.focused[0]) is None:
            self.focused = self._site = None  # its protein was removed
        if self.pose is not None:
            protein = self.display(self.pose.name)
            if protein is None or protein.structure is not self._posed_in:
                self.pose = self.pose_row = None
        self._scene = None
        self.proteinsChanged.emit()
        self.sceneChanged.emit()
        if replaced:
            self.edited.emit()

    def _on_table_event(self, event) -> None:
        if not isinstance(event, Changed):
            return
        if any(isinstance(c, AssetsChanged) for c in event.changes):
            self.reload()
        edited = [c for c in event.changes
                  if isinstance(c, CellsChanged) and c.column in self.row_columns]  # fmt: skip
        if edited:
            self.rowsEdited.emit(np.unique(np.concatenate([c.oids for c in edited])))

    def display(self, protein: str) -> displays.Display | None:
        return next((d for d in self.displays if d.id == protein), None)

    def shown(self, id: str) -> displays.Display | None:
        """The protein's display, or the pose's (for hovering and measuring)."""
        if self.pose is not None and id == self.pose.id:
            return self.pose
        return self.display(id)

    def all_displays(self) -> list[displays.Display]:
        """The proteins' displays and the pose's."""
        return self.displays + ([self.pose] if self.pose is not None else [])

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
            if self.pose is not None:
                protein = self.display(self.pose.name)
                self.pose.visible = protein is not None and protein.visible
            self._scene = displays.compose(self.all_displays())
            if self.builder is not None:
                self.builder.schedule(self._scene.displays)
        return self._scene

    def built(self) -> None:
        """A part built in the background arrived: draw it."""
        self._scene = None
        self.sceneChanged.emit()

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
                out[shown.id] = (
                    select(shown.structure, text, named=named, model=shown.model) & shown.active
                )
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
        if name == SELECTION:
            raise SelectionError(f"'{name}' is taken: $selection is what is selected now", 0)
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
        out[SELECTION] = shown.selected.copy()
        return out

    # editing: each edit is one step of the document's history
    def edit(
        self, label: str, structures: Mapping[str, Structure], extra_ops: Sequence = ()
    ) -> None:
        """Replace proteins' structures (by protein id) in one undoable step, with
        *extra_ops* (rows renamed with a ligand)."""
        ops = []
        for protein, structure in structures.items():
            asset = self.table.assets.get(protein)
            shown = self.display(protein)
            if asset is not None and shown is not None and structure is not shown.structure:
                ops.append(PutAsset(protein, asset, replace(asset, payload=structure)))
        if ops:
            self.table.edit(label, [*ops, *extra_ops])

    def _chosen(self) -> list[displays.Display]:
        return [d for d in self.displays if d.visible and d.selected.any()]

    def delete_selected(self) -> int:
        """Delete the selected atoms (each with its other locations); how many went. Raises
        EditError."""
        chosen = self._chosen()
        made = {d.id: edits.delete(d.structure, d.selected) for d in chosen}
        count = sum(d.structure.n_atoms - made[d.id].n_atoms for d in chosen)
        where = f" from {chosen[0].name}" if len(chosen) == 1 else ""
        self.edit(f"Delete {_atoms(count)}{where}", made)
        return count

    def set_element(self, symbol: str) -> None:
        """Make the selected atoms *symbol*. Raises EditError."""
        chosen = self._chosen()
        made = {d.id: edits.set_element(d.structure, d.selected, symbol) for d in chosen}
        count = sum(int(d.selected.sum()) for d in chosen)
        self.edit(f"Make {_atoms(count)} {symbol.strip().capitalize()}", made)

    def set_charge(self, charge: int) -> None:
        """Give the selected atoms the formal charge *charge*. Raises EditError."""
        chosen = self._chosen()
        made = {d.id: edits.set_charge(d.structure, d.selected, charge) for d in chosen}
        count = sum(int(d.selected.sum()) for d in chosen)
        self.edit(f"Set the charge of {_atoms(count)} to {charge:+d}", made)

    def selected_residues(self) -> tuple[displays.Display, np.ndarray] | None:
        """The protein holding the whole selection, and the residues selected in it (None when
        nothing is selected, or more than one protein is)."""
        chosen = self._chosen()
        if len(chosen) != 1:
            return None
        shown = chosen[0]
        return shown, np.unique(shown.structure.residue_of[shown.selected])

    def selected_chain(self) -> tuple[displays.Display, int] | None:
        """The protein and chain holding the whole selection, if it is in one chain."""
        found = self.selected_residues()
        if found is None:
            return None
        shown, residues = found
        chains = np.unique(shown.structure.chain_of[residues])
        return (shown, int(chains[0])) if len(chains) == 1 else None

    def rename_residues(self, protein: str, residues: np.ndarray, name: str) -> None:
        """Raises EditError."""
        shown = self._shown_or_raise(protein)
        made = edits.rename_residues(shown.structure, residues, name)
        if len(residues) == 1:
            what = shown.structure.residue_label(int(residues[0]))
        else:
            what = f"{len(residues)} residues"
        self._relabelled(shown, made, residues, f"Rename {what} to {name.strip().upper()}")

    def rename_chain(self, protein: str, chain: int, name: str) -> None:
        """Raises EditError."""
        shown = self._shown_or_raise(protein)
        made = edits.rename_chain(shown.structure, chain, name)
        span = shown.structure.residues_of(chain)
        residues = np.arange(span.start, span.stop)
        old = shown.structure.chain_name[chain]
        self._relabelled(shown, made, residues, f"Rename chain {old} to {name.strip()}")

    def renumber(self, protein: str, residues: np.ndarray, first: int) -> None:
        """Raises EditError."""
        shown = self._shown_or_raise(protein)
        made = edits.renumber(shown.structure, residues, first)
        label = f"Renumber {len(residues)} residue(s) from {first}"
        self._relabelled(shown, made, residues, label)

    def rows_naming(self, protein: str, atoms: np.ndarray) -> int:
        """How many table rows name a residue these atoms of *protein* belong to (0 without a
        host that links rows)."""
        shown = self.display(protein)
        if self._count_named is None or shown is None or not atoms.any():
            return 0
        s = shown.structure
        labels = {s.residue_label(int(r)) for r in np.unique(s.residue_of[atoms])}
        return self._count_named(protein, labels)

    def _shown_or_raise(self, protein: str) -> displays.Display:
        shown = self.display(protein)
        if shown is None:
            raise edits.EditError(f"There is no protein {protein}")
        return shown

    def _relabelled(self, shown, made: Structure, residues, label: str) -> None:
        """Apply *made*, with the linked rows' ligand labels following the residues."""
        before = shown.structure
        labels = {before.residue_label(int(r)): made.residue_label(int(r)) for r in residues}
        labels = {old: new for old, new in labels.items() if old != new}
        extra = self._relabel(shown.id, labels) if self._relabel and labels else ()
        self.edit(label, {shown.id: made}, extra)

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

    @property
    def hydrogens(self) -> str:
        """Which hydrogens are drawn (all, polar, none), as the first protein has it."""
        return self.displays[0].hydrogens if self.displays else displays.POLAR_H

    def show_hydrogens(self, which: str) -> None:
        for shown in self.displays:
            shown.show_hydrogens(which)
        self._changed(rebuilt=True)

    @property
    def transparency(self) -> float:
        """How see-through surfaces are (0 opaque), as the first protein has it."""
        return self.displays[0].transparency if self.displays else displays.SURFACE_TRANSPARENCY

    def set_transparency(self, transparency: float) -> None:
        for shown in self.displays:
            shown.transparency = transparency
        self._changed(rebuilt=False)

    def choose_altloc(self, protein: str, residue: int, letter: str) -> None:
        shown = self.display(protein)
        if shown is not None:
            shown.choose_altloc(residue, letter)
            self._changed(rebuilt=True)

    def set_model(self, protein: str, model: int) -> None:
        """Show another of the protein's models (NMR models; trajectory frames later)."""
        shown = self.display(protein)
        if shown is not None and shown.set_model(model):
            if self.pose is not None and self.pose.name == protein:
                self.pose.set_model(model)
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
        from chemsuite.protein.chem.interactions import describe as say

        shown = self.shown(protein)
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
        site = select(s, FOCUS_SITE, named=named, model=shown.model) & shown.active
        waters = select(s, FOCUS_NEAR, named=named, model=shown.model) & shown.active
        self.pose = self.pose_row = None
        self._pocket(shown, ligand, site, waters, restyle)
        if restyle:
            shown.show_interactions(ligand)
        shown.ghost = ~(site | waters | ligand)
        self.focused = (protein, int(s.uid[s.atoms_of(residue).start]))
        self._changed(rebuilt=True)
        self.focusChanged.emit()

    def _pocket(self, shown, ligand, site, waters, restyle: bool) -> None:
        """Fade every protein, show this one's pocket and hide its other small molecules."""
        for other in self.displays:
            other.ghost = np.zeros(other.structure.n_atoms, bool)
            other.show_pocket(None)
            if restyle:
                other.show_interactions(None)
        if restyle:
            s = shown.structure
            polymer = np.isin(s.res_kind[s.residue_of], (AMINO, NUCLEIC))
            shown.show(~(site | waters | ligand | polymer), scenes.HIDDEN)
            shown.show(site, scenes.STICK)
            shown.show(waters, scenes.BALL_AND_STICK)
            shown.show(ligand, scenes.BALL_AND_STICK)
        self._site = site
        shown.show_pocket(site if self.pocket_surface else None)

    def show_pose(self, protein: str, mol, oid: int | None = None) -> bool:
        """Show ``mol``'s own 3D pose (an RDKit molecule) in ``protein``'s pocket, as Focus
        shows a ligand: its pocket as sticks, the protein faded, its interactions. A ligand of
        the protein it overlaps (the one it stands in for) stays, faded, for comparison."""
        shown = self.display(protein)
        if shown is None:
            return False
        s = shown.structure
        pose = poses.from_mol(mol)
        heavy = pose.xyz()[pose.element != 1]
        replaced = poses.overlapping(s, heavy, shown.active, shown.model)
        both = poses.merged(s.subset(~replaced), pose)
        mine = np.zeros(both.n_atoms, bool)
        mine[s.n_atoms - int(replaced.sum()) :] = True
        display = displays.Display(
            POSE + protein,
            protein,
            both,
            model=shown.model,
            hydrogens=shown.hydrogens,
            altlocs={} if replaced.any() else dict(shown.altlocs),
            deferred=self.builder is not None,
        )
        display.show(~mine, scenes.HIDDEN)
        display.show(mine, scenes.BALL_AND_STICK)
        display.show_cartoon(np.ones(both.n_atoms, bool), False)
        display.show_interactions(mine)
        named = {"focus": mine}
        site = select(both, FOCUS_SITE, named=named, model=display.model)
        waters = select(both, FOCUS_NEAR, named=named, model=display.model)

        def in_protein(mask: np.ndarray) -> np.ndarray:
            out = np.zeros(s.n_atoms, bool)
            found = s.indices(both.uid[mask & ~mine])
            out[found[found >= 0]] = True
            return out & shown.active

        site, waters = in_protein(site), in_protein(waters)
        self.focused = None
        self._pocket(shown, replaced, site, waters, restyle=True)
        shown.show(replaced, scenes.STICK)
        shown.ghost = ~(site | waters)
        self.pose, self.pose_row, self._posed_in = display, oid, s
        self._changed(rebuilt=True)
        self.focusChanged.emit()
        return True

    def focused_atoms(self) -> tuple[displays.Display, np.ndarray] | None:
        """The focused ligand's display and its pocket's atoms (ligand, site, waters); for a
        pose, its protein's display and the pocket's."""
        if self.pose is not None:
            shown = self.display(self.pose.name)
            if shown is None or not shown.visible:
                return None
            return shown, ~shown.ghost
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
        if self.focused is None and self.pose is None:
            return
        self.focused = None
        self.pose = self.pose_row = None
        self._site = None
        for shown in self.displays:
            shown.ghost = np.zeros(shown.structure.n_atoms, bool)
            shown.show_pocket(None)
        self._changed(rebuilt=True)
        self.focusChanged.emit()

    def show_pocket_surface(self, on: bool) -> None:
        """Whether Focus draws a skin over the pocket (coloured by hydrophobicity)."""
        self.pocket_surface = on
        found = self.focused_atoms()
        if found is not None and self._site is not None:
            shown, _atoms = found
            if len(self._site) == shown.structure.n_atoms:
                shown.show_pocket(self._site if on else None)
        self._changed(rebuilt=True)

    # labels
    def label(self, kind: str, masks: Mapping[str, Mask] | None = None) -> None:
        """Label the selected residues (``RESIDUE_LABEL``) or atoms (``ATOM_LABEL``), or those
        of ``masks`` (atoms per protein)."""
        masks = self.selection() if masks is None else masks
        for shown in self.displays:
            mask = masks.get(shown.id)
            if mask is None or not mask.any():
                continue
            s = shown.structure
            if kind == RESIDUE_LABEL:
                atoms = s.res_start[:-1][np.unique(s.residue_of[mask])]
            else:
                atoms = np.flatnonzero(mask)
            for uid in s.uid[atoms].tolist():
                self.labels[(shown.id, int(uid), kind)] = None
        self.labelsChanged.emit()

    def clear_labels(self) -> None:
        self.labels = {}
        self.labelsChanged.emit()

    def labelled(self) -> tuple[np.ndarray, list[str]]:
        """Where each label goes that can be drawn, and its text: on a residue's CA (or the
        middle of its atoms) and an atom itself; none on a protein hidden, an atom faded by
        Focus or a residue not drawn at all."""
        scene = self.scene
        drawn = np.zeros(int(scene.atom_offsets[-1]), bool)
        drawn[scene.geometry.atoms] = True
        drawn[scene.pick_atoms] = True  # a cartoon's residues, through their trace atoms
        positions, texts = [], []
        for protein, uid, kind in self.labels:
            shown = self.display(protein)
            if shown is None or shown not in scene.displays:
                continue
            s = shown.structure
            index = int(s.indices(np.array([uid]))[0])
            if index < 0 or shown.ghost[index]:
                continue
            offset = int(scene.global_atoms(shown, np.zeros(1, np.int64))[0])
            if kind == RESIDUE_LABEL:
                residue = int(s.residue_of[index])
                span = s.atoms_of(residue)
                atoms = np.arange(span.start, span.stop)[shown.active[span]]
                if not len(atoms) or not drawn[atoms + offset].any():
                    continue
                ca = atoms[s.name[atoms] == "CA"]
                where = shown.xyz[ca[0]] if len(ca) else shown.xyz[atoms].mean(axis=0)
                text = f"{s.res_name[residue]} {s.res_seq[residue]}{s.icode[residue]}"
            else:
                if not drawn[index + offset]:
                    continue
                where, text = shown.xyz[index], str(s.name[index])
            positions.append(where)
            texts.append(text)
        return np.array(positions, np.float32).reshape(-1, 3), texts

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
        shown = self.shown(point[0])
        scene = self.scene
        if shown is None or shown not in scene.displays:
            return None
        index = int(shown.structure.indices(np.array([point[1]]))[0])
        if index < 0:
            return None
        return int(scene.global_atoms(shown, np.array([index]))[0]), shown.xyz[index]

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
            "pocket_surface": self.pocket_surface,
            "labels": [list(label) for label in self.labels],
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
        self.pocket_surface = state.get("pocket_surface", True) is not False
        self.labels = {
            (str(item[0]), item[1], item[2]): None
            for item in state.get("labels") or []
            if isinstance(item, list) and len(item) == 3 and isinstance(item[1], int)
            and item[2] in (ATOM_LABEL, RESIDUE_LABEL)
        }  # fmt: skip
        self.labelsChanged.emit()
        focused = _points([[state.get("focused")]])
        self.focused = focused[0][0] if focused else None
        if self.focused is not None:
            place = self.focused_residue()
            if place is not None:
                self.focus(*place, restyle=False)
        self._scene = None
        self.proteinsChanged.emit()
        self.sceneChanged.emit()


def _atoms(count: int) -> str:
    return f"{count} atom" if count == 1 else f"{count} atoms"


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
