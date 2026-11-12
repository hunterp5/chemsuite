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
"""How each protein in a view is shown, and all of them as one scene.

A :class:`Display` holds one structure's representation: a style per atom, which residues
show as cartoon and which atoms make up a surface, colours, the selection and the ghosted
atoms, the alternate locations chosen and which model is shown (an NMR model; a trajectory
can be a memory-mapped ``coords`` array, one model per frame). It builds its geometry and
meshes only when what they depend on changed (each part has its own version), so a selection
or a colour change rebuilds nothing but the look, and a model change nothing but what depends
on positions.

A display that *defers* (the viewer has a way to run work in the background) never builds its
slow parts inline: surfaces, interactions and the cartoon of a large protein. It draws the last
one it accepted while :meth:`Display.pending` hands out the work (a :class:`Build`, whose
arguments are copies taken when it was asked for) and :meth:`Display.accept` takes a result
only if that part has not changed since. Without deferring, everything builds inline, as
scripts and tests want.

:func:`compose` puts the visible displays together for the view: atoms and residues are
numbered one after the other (a display's ``atom_offset`` and ``residue_offset``), so a pick
on the combined scene maps back to its protein with :meth:`Scene.locate`.

Representation state round-trips through :meth:`Display.state` as run-length lists, small
enough for a saved document.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from chemsuite.molview.geometry import HYDROGENS, POLAR_H

from ..chem import interactions as contacts
from ..structure.model import AMINO, NUCLEIC, Structure
from . import cartoon as cartoons
from . import colours as schemes
from . import mesh as meshes
from . import scene as scenes
from . import surface as surfaces

#: how see-through a surface is drawn (0 opaque, 1 invisible), and the choices offered
SURFACE_TRANSPARENCY = 0.4
TRANSPARENCIES = (0.0, 0.2, 0.4, 0.6, 0.8)
CARTOON_PICK_RADIUS = 1.2
#: the parts a deferring display builds in the background
DEFERRED = ("surface", "pocket", "interactions", "cartoon")
#: a cartoon of more residues than this is built in the background
LARGE_CARTOON = 2000
_EMPTY: dict[str, Callable[[], Any]] = {
    "surface": meshes.Mesh.empty,
    "pocket": meshes.Mesh.empty,
    "cartoon": meshes.Mesh.empty,
    "interactions": contacts.Interactions.empty,
}


@dataclass(frozen=True, eq=False)
class Build:
    """A part to build off the GUI thread: ``fn(*args)``, safe to call on any thread (the
    structure is immutable and the arguments are copies)."""

    part: str
    version: int
    fn: Callable[..., Any]
    args: tuple

    def __call__(self) -> Any:
        return self.fn(*self.args)


@dataclass(eq=False)
class Display:
    id: str
    name: str
    structure: Structure
    visible: bool = True
    scheme: str = schemes.ELEMENT
    surface_scheme: str = schemes.HYDROPHOBICITY
    altlocs: dict[int, str] = field(default_factory=dict)
    model: int = 0  # which of the structure's coordinate sets is shown
    deferred: bool = False  # slow parts are built through pending() and accept()
    hydrogens: str = POLAR_H  # which hydrogens are drawn: all, polar, none
    transparency: float = SURFACE_TRANSPARENCY  # of the surface
    styles: np.ndarray = field(init=False)
    cartoon: np.ndarray = field(init=False)  # bool per residue
    surface: np.ndarray = field(init=False)  # bool per atom
    pocket: np.ndarray = field(init=False)  # bool per atom: the focused pocket's skin
    colours: np.ndarray = field(init=False)  # rgba per atom
    selected: np.ndarray = field(init=False)  # bool per atom
    ghost: np.ndarray = field(init=False)  # bool per atom
    query: np.ndarray = field(init=False)  # bool per atom: whose interactions are shown
    _cache: dict = field(default_factory=dict, init=False, repr=False)
    _latest: dict = field(default_factory=dict, init=False, repr=False)  # last built, per part
    _rgba: dict = field(default_factory=dict, init=False, repr=False)  # surface colours, by scheme
    _versions: dict = field(init=False, repr=False)

    def __post_init__(self) -> None:
        s = self.structure
        polymer = np.isin(s.res_kind, (AMINO, NUCLEIC))
        self.cartoon = polymer.copy()
        self.styles = scenes.default_styles(s, self.active, cartoon=bool(polymer.any()))
        self.surface = np.zeros(s.n_atoms, bool)
        self.pocket = np.zeros(s.n_atoms, bool)
        self.colours = schemes.colour(s, self.scheme)
        self.selected = np.zeros(s.n_atoms, bool)
        self.ghost = np.zeros(s.n_atoms, bool)
        self.query = np.zeros(s.n_atoms, bool)
        self._versions = {"geometry": 0, "cartoon": 0, "surface": 0, "pocket": 0, "trace": 0,
                          "interactions": 0}  # fmt: skip

    # what changed
    def touch(self, *parts: str) -> None:
        for part in parts:
            self._versions[part] += 1

    @property
    def active(self) -> np.ndarray:
        return self.structure.active(self.altlocs)

    @property
    def xyz(self) -> np.ndarray:
        """The positions shown: the current model's."""
        return self.structure.xyz(self.model)

    @property
    def n_models(self) -> int:
        return self.structure.n_models

    def set_model(self, model: int) -> bool:
        """Show another model; False when it is already shown or out of range. Only what
        depends on positions is rebuilt: colours, styles and the selection stay."""
        if not 0 <= model < self.n_models or model == self.model:
            return False
        self.model = model
        self.touch("geometry", "cartoon", "surface", "pocket", "interactions")
        return True

    # representations
    def show(self, atoms: np.ndarray, style: int) -> None:
        """Show these atoms in ``style`` (HIDDEN hides them); other locations stay hidden."""
        self.styles[atoms & self.active] = style
        self.styles[~self.active] = scenes.HIDDEN
        self.touch("geometry")

    def show_cartoon(self, atoms: np.ndarray, on: bool = True) -> None:
        residues = np.unique(self.structure.residue_of[atoms])
        polymer = np.isin(self.structure.res_kind[residues], (AMINO, NUCLEIC))
        self.cartoon[residues[polymer]] = on
        self.touch("cartoon")

    def show_surface(self, atoms: np.ndarray, on: bool = True) -> None:
        heavy = self.structure.element != 1
        self.surface[atoms & heavy & self.active] = on
        self.touch("surface")

    def show_pocket(self, atoms: np.ndarray | None) -> None:
        """Draw a skin over these atoms (a focused ligand's pocket), or none."""
        heavy = self.structure.element != 1
        pocket = np.zeros(self.structure.n_atoms, bool) if atoms is None else atoms & heavy
        if not np.array_equal(pocket, self.pocket):
            self.pocket = pocket
            self.touch("pocket")

    def show_hydrogens(self, which: str) -> None:
        """Draw all, the polar or no hydrogens (each as its heavy atom is shown)."""
        if which in HYDROGENS and which != self.hydrogens:
            self.hydrogens = which
            self.touch("geometry")

    def surface_alpha(self) -> float:
        return 1.0 - self.transparency

    def surface_colours(self, pocket: bool = False) -> np.ndarray:
        """Per atom, the surface's colours (a pocket's: hydrophobicity) at its alpha; kept, as
        they change only with the scheme or the transparency (do not modify)."""
        key = (schemes.HYDROPHOBICITY if pocket else self.surface_scheme, self.surface_alpha())
        rgba = self._rgba.get(key)
        if rgba is None:
            rgba = schemes.colour(self.structure, key[0])
            rgba[:, 3] = key[1]
            self._rgba = {key: rgba, **{k: v for k, v in self._rgba.items() if k[0] != key[0]}}
        return rgba

    def colour_by(self, scheme: str, atoms: np.ndarray | None = None) -> None:
        coloured = schemes.colour(self.structure, scheme)
        if atoms is None:
            self.scheme = scheme
            self.colours = coloured
        else:
            self.colours[atoms] = coloured[atoms]

    def colour(self, atoms: np.ndarray, rgb: tuple[float, float, float]) -> None:
        self.colours[atoms, :3] = rgb

    def choose_altloc(self, residue: int, letter: str) -> None:
        was = self.active
        self.altlocs[residue] = letter
        now = self.active
        span = self.structure.atoms_of(residue)
        # the newly shown location takes the style its partner had
        for i in range(span.start, span.stop):
            if now[i] and not was[i]:
                twins = np.flatnonzero(
                    (self.structure.name[span] == self.structure.name[i]) & was[span]
                )
                self.styles[i] = self.styles[span.start + twins[0]] if len(twins) else 0
        self.styles[~now] = scenes.HIDDEN
        self.surface &= now
        self.pocket &= now
        self.touch("geometry", "cartoon", "surface", "pocket", "trace", "interactions")

    def show_interactions(self, atoms: np.ndarray | None) -> None:
        """Show the interactions of these atoms with the rest (None: none)."""
        self.query = np.zeros(self.structure.n_atoms, bool) if atoms is None else atoms.copy()
        self.touch("interactions")

    # built parts, cached by version
    def _work(self, part: str) -> tuple[Callable[..., Any], tuple]:
        """What builds a part, with copies of what it reads."""
        s = self.structure
        if part == "geometry":
            styles = scenes.shown_styles(s, self.styles, self.active, self.hydrogens)
            return scenes.geometry, (s, styles, self.model)
        if part == "cartoon":
            return cartoons.build, (s, np.flatnonzero(self.cartoon), self.model, self.active)
        if part in ("surface", "pocket"):
            atoms = self.surface if part == "surface" else self.pocket
            return surfaces.build, (s, np.flatnonzero(atoms), self.model)
        if part == "trace":
            return cartoons.trace, (s, self.active)
        return contacts.find, (s, self.query.copy(), self.active, self.model)

    def _inline(self, part: str, args: tuple) -> bool:
        """Whether a part is built at once: always unless deferring, and then when it is quick."""
        if not self.deferred or part not in DEFERRED:
            return True
        if part == "cartoon":
            return len(args[1]) <= LARGE_CARTOON
        if part in ("surface", "pocket"):
            return not len(args[1])
        return False

    def _part(self, part: str) -> Any:
        key = (part, self._versions[part])
        if key in self._cache:
            return self._cache[key]
        fn, args = self._work(part)
        if not self._inline(part, args):
            latest = self._latest.get(part)
            return latest if latest is not None else _EMPTY[part]()
        result = fn(*args)
        self._store(key, result)
        return result

    def _store(self, key: tuple[str, int], result: Any) -> None:
        part = key[0]
        self._cache = {k: v for k, v in self._cache.items() if k[0] != part}
        self._cache[key] = result
        self._latest[part] = result

    def pending(self) -> list[Build]:
        """The deferred parts that are out of date, as work to run in the background."""
        out = []
        for part in DEFERRED:
            if part == "interactions" and not self.query.any():
                continue
            version = self._versions[part]
            if (part, version) in self._cache:
                continue
            fn, args = self._work(part)
            if not self._inline(part, args):
                out.append(Build(part, version, fn, args))
        return out

    def accept(self, part: str, version: int, result: Any) -> bool:
        """Take a part built in the background; False (dropped) when the part has changed
        since its work was handed out."""
        if version != self._versions[part]:
            return False
        self._store((part, version), result)
        return True

    def geometry(self) -> scenes.Geometry:
        return self._part("geometry")

    def cartoon_mesh(self) -> meshes.Mesh:
        return self._part("cartoon")

    def surface_mesh(self) -> meshes.Mesh:
        return self._part("surface")

    def pocket_mesh(self) -> meshes.Mesh:
        return self._part("pocket")

    def trace(self) -> tuple[np.ndarray, np.ndarray]:
        """Cartoon residues and their trace atoms (what picks the cartoon)."""
        residues, atoms = self._part("trace")
        shown = self.cartoon[residues]
        return residues[shown], atoms[shown]

    def interactions(self) -> contacts.Interactions:
        if not self.query.any():
            return contacts.Interactions.empty()
        return self._part("interactions")

    def flags(self) -> np.ndarray:
        return scenes.flags(self.structure.n_atoms, selected=self.selected, ghost=self.ghost)

    def residue_look(self) -> tuple[np.ndarray, np.ndarray]:
        """Colour and flags per residue, for the cartoon: the trace atom's colour (or the
        residue's first atom's), selected when any of its atoms is."""
        s = self.structure
        colour = self.colours[s.res_start[:-1]].copy()
        residues, atoms = self.trace()
        colour[residues] = self.colours[atoms]
        flags = np.zeros(s.n_residues, np.float32)
        flags[np.unique(s.residue_of[self.selected])] += scenes.SELECTED
        flags[np.unique(s.residue_of[self.ghost])] += scenes.GHOST
        return colour, flags

    # an edited structure
    def carry(self, old: Display) -> None:
        """Take *old*'s look for this display of an edited structure, atom by atom through
        their uids (the run lists of :meth:`state` are by position, which an edit shifts).
        Atoms new to this structure (brought back by an undo) keep the default look; an atom
        whose element changed takes its scheme's colour."""
        s, o = self.structure, old.structure
        self.visible, self.scheme, self.surface_scheme = old.visible, old.scheme, old.surface_scheme
        self.hydrogens, self.transparency = old.hydrogens, old.transparency
        self.model = old.model if old.model < self.n_models else 0
        source = o.indices(s.uid)  # each atom's index in the old structure, or -1
        kept = source >= 0
        was = source[kept]
        first = source[s.res_start[:-1]]  # each residue's first atom, in the old structure
        old_residue = np.where(first >= 0, o.residue_of[np.maximum(first, 0)], -1)
        residue_now = {int(r): i for i, r in enumerate(old_residue) if r >= 0}
        self.altlocs = {residue_now[r]: letter for r, letter in old.altlocs.items()
                        if r in residue_now}  # fmt: skip
        carried = old_residue >= 0
        self.cartoon[carried] = old.cartoon[old_residue[carried]]
        self.styles[kept] = old.styles[was]
        self.surface[kept] = old.surface[was]
        self.query[kept] = old.query[was]
        self.selected[kept] = old.selected[was]
        self.colours = schemes.colour(s, self.scheme)
        same = kept.copy()
        same[kept] = s.element[kept] == o.element[was]
        self.colours[same] = old.colours[source[same]]
        self.styles[~self.active] = scenes.HIDDEN
        self.surface &= self.active
        self.touch("geometry", "cartoon", "surface", "pocket", "trace", "interactions")

    # saving
    def state(self) -> dict:
        return {
            "visible": self.visible,
            "scheme": self.scheme,
            "surface_scheme": self.surface_scheme,
            "altlocs": {str(k): v for k, v in self.altlocs.items()},
            "model": self.model,
            "hydrogens": self.hydrogens,
            "transparency": self.transparency,
            "styles": runs(self.styles),
            "cartoon": runs(self.cartoon.astype(np.uint8)),
            "surface": runs(self.surface.astype(np.uint8)),
            "interactions": runs(self.query.astype(np.uint8)),
            "custom_colours": runs(self._custom_colours()),
        }

    def _custom_colours(self) -> np.ndarray:
        """Per atom, the colour given by hand as 0xRRGGBB, or -1 where the scheme's shows."""
        own = np.full(self.structure.n_atoms, -1, np.int64)
        scheme = schemes.colour(self.structure, self.scheme)
        changed = np.any(np.abs(self.colours[:, :3] - scheme[:, :3]) > 1e-6, axis=1)
        rgb = np.clip(np.round(self.colours[changed, :3] * 255), 0, 255).astype(np.int64)
        own[changed] = rgb[:, 0] << 16 | rgb[:, 1] << 8 | rgb[:, 2]
        return own

    def set_state(self, state: Mapping) -> None:
        n, r = self.structure.n_atoms, self.structure.n_residues
        self.visible = bool(state.get("visible", True))
        self.altlocs = {int(k): str(v) for k, v in dict(state.get("altlocs", {})).items()}
        model = state.get("model", 0)
        self.model = model if isinstance(model, int) and 0 <= model < self.n_models else 0
        hydrogens = state.get("hydrogens", POLAR_H)
        self.hydrogens = hydrogens if hydrogens in HYDROGENS else POLAR_H
        transparency = state.get("transparency", SURFACE_TRANSPARENCY)
        valid = isinstance(transparency, int | float) and 0 <= transparency < 1
        self.transparency = float(transparency) if valid else SURFACE_TRANSPARENCY
        styles = unruns(state.get("styles"), n)
        if styles is not None:
            self.styles = styles.astype(np.uint8)
        cartoon = unruns(state.get("cartoon"), r)
        if cartoon is not None:
            self.cartoon = cartoon.astype(bool)
        surface = unruns(state.get("surface"), n)
        if surface is not None:
            self.surface = surface.astype(bool)
        query = unruns(state.get("interactions"), n)
        self.query = np.zeros(n, bool) if query is None else query.astype(bool)
        scheme = str(state.get("scheme", self.scheme))
        self.scheme = scheme if scheme in schemes.SCHEMES else schemes.ELEMENT
        self.colours = schemes.colour(self.structure, self.scheme)
        own = unruns(state.get("custom_colours"), n)
        if own is not None:
            given = own >= 0
            rgb = own[given] & 0xFFFFFF
            self.colours[given, :3] = np.stack([rgb >> 16, rgb >> 8 & 255, rgb & 255], 1) / 255
        surface_scheme = str(state.get("surface_scheme", self.surface_scheme))
        if surface_scheme in schemes.SCHEMES:
            self.surface_scheme = surface_scheme
        self.touch("geometry", "cartoon", "surface", "interactions")


def runs(values: np.ndarray) -> list[int]:
    """``[value, count, value, count, ...]``: a per-atom array, small."""
    values = np.asarray(values)
    if not len(values):
        return []
    starts = np.flatnonzero(np.concatenate([[True], values[1:] != values[:-1]]))
    counts = np.diff(np.append(starts, len(values)))
    return np.stack([values[starts], counts], axis=1).ravel().astype(int).tolist()


def unruns(data, length: int) -> np.ndarray | None:
    """The array :func:`runs` gave, or None when it does not fit ``length`` atoms."""
    if not isinstance(data, list) or len(data) % 2:
        return None
    pairs = np.array(data, dtype=np.int64).reshape(-1, 2)
    if pairs[:, 1].sum() != length or (pairs[:, 1] < 0).any():
        return None
    return np.repeat(pairs[:, 0], pairs[:, 1])


@dataclass(frozen=True)
class Scene:
    """The visible displays as one: what the view draws."""

    displays: tuple[Display, ...]
    atom_offsets: np.ndarray  # int64 (d + 1)
    residue_offsets: np.ndarray  # int64 (d + 1)
    geometry: scenes.Geometry
    sphere_look: np.ndarray
    cylinder_look: np.ndarray
    cartoon: meshes.Mesh  # owners: global residues
    cartoon_look: np.ndarray
    surface: meshes.Mesh  # owners: global atoms
    surface_look: np.ndarray
    pocket: meshes.Mesh  # owners: global atoms
    pocket_look: np.ndarray
    pick_atoms: np.ndarray  # global atoms of the cartoon trace
    pick_centres: np.ndarray
    interactions: tuple[contacts.Interactions, ...] = ()  # per display
    lines: meshes.Mesh = field(default_factory=meshes.Mesh.empty)  # owners: global atoms
    lines_look: np.ndarray = field(default_factory=lambda: np.zeros((0, meshes.LOOK_FLOATS)))

    def locate(self, atom: int) -> tuple[Display, int] | None:
        """The display and its own atom index for a global atom index."""
        if atom < 0 or atom >= self.atom_offsets[-1]:
            return None
        d = int(np.searchsorted(self.atom_offsets, atom, side="right") - 1)
        return self.displays[d], atom - int(self.atom_offsets[d])

    def global_atoms(self, display: Display, atoms: np.ndarray) -> np.ndarray:
        d = self.displays.index(display)
        return np.asarray(atoms, np.int64) + int(self.atom_offsets[d])

    def global_residues(self, display: Display, residues: np.ndarray) -> np.ndarray:
        d = self.displays.index(display)
        return np.asarray(residues, np.int64) + int(self.residue_offsets[d])


def compose(displays: list[Display] | tuple[Display, ...]) -> Scene:
    shown = tuple(d for d in displays if d.visible)
    atom_offsets = np.cumsum([0] + [d.structure.n_atoms for d in shown]).astype(np.int64)
    residue_offsets = np.cumsum([0] + [d.structure.n_residues for d in shown]).astype(np.int64)
    parts = [d.geometry() for d in shown]
    geometry = scenes.Geometry(
        atoms=_cat([p.atoms + o for p, o in zip(parts, atom_offsets, strict=False)], np.int64),
        spheres=_cat([p.spheres for p in parts], np.float32, scenes.SPHERE_FLOATS),
        bonds=_cat([p.bonds for p in parts], np.int64),
        ends=_cat([p.ends + o for p, o in zip(parts, atom_offsets, strict=False)], np.int64, 2),
        cylinders=_cat([p.cylinders for p in parts], np.float32, scenes.CYLINDER_FLOATS),
    )
    colours = _cat([d.colours for d in shown], np.float32, 4)
    flags = _cat([d.flags() for d in shown], np.float32)
    cartoon = meshes.join([_shifted(d.cartoon_mesh(), o) for d, o in
                           zip(shown, residue_offsets, strict=False)])  # fmt: skip
    looks = [d.residue_look() for d in shown]
    residue_colours = _cat([c for c, _ in looks], np.float32, 4)
    residue_flags = _cat([f for _, f in looks], np.float32)
    surface = meshes.join([_shifted(d.surface_mesh(), o) for d, o in
                           zip(shown, atom_offsets, strict=False)])  # fmt: skip
    surface_colours = _surface_colours(shown, surface)
    pocket = meshes.join([_shifted(d.pocket_mesh(), o) for d, o in
                          zip(shown, atom_offsets, strict=False)])  # fmt: skip
    pocket_colours = _surface_colours(shown, pocket, pocket=True)
    traces = [d.trace()[1] for d in shown]
    pick_atoms = _cat([t + o for t, o in zip(traces, atom_offsets, strict=False)], np.int64)
    pick_centres = _cat([d.xyz[t] for d, t in zip(shown, traces, strict=True)], np.float32, 3)
    found = tuple(d.interactions() for d in shown)
    lines, lines_look = interaction_lines(found, atom_offsets)
    return Scene(
        displays=shown,
        atom_offsets=atom_offsets,
        residue_offsets=residue_offsets,
        geometry=geometry,
        sphere_look=scenes.sphere_look(geometry, colours, flags),
        cylinder_look=scenes.cylinder_look(geometry, colours, flags),
        cartoon=cartoon,
        cartoon_look=meshes.look(cartoon, residue_colours, residue_flags),
        surface=surface,
        surface_look=meshes.look(surface, surface_colours, flags),
        pocket=pocket,
        pocket_look=meshes.look(pocket, pocket_colours, flags),
        pick_atoms=pick_atoms,
        pick_centres=pick_centres,
        interactions=found,
        lines=lines,
        lines_look=lines_look,
    )


def interaction_lines(found, atom_offsets) -> tuple[meshes.Mesh, np.ndarray]:
    """Dashed lines for the interactions, coloured by kind; each line is owned by its
    query atom, so hovering that atom lights its lines."""
    found = [f for f in found]
    if not any(len(f) for f in found):
        return meshes.Mesh.empty(), np.zeros((0, meshes.LOOK_FLOATS), np.float32)
    starts = np.concatenate([f.start for f in found])
    ends = np.concatenate([f.end for f in found])
    owners = np.concatenate([f.a + o for f, o in zip(found, atom_offsets, strict=False)])
    kinds = np.concatenate([f.kind for f in found])
    lines = meshes.dashes(starts, ends)
    palette = np.array([(*contacts.COLOURS[k], 1.0) for k in sorted(contacts.COLOURS)],
                       np.float32)  # fmt: skip
    look = np.zeros((len(lines), meshes.LOOK_FLOATS), np.float32)
    look[:, :4] = palette[kinds[lines.owner]]
    return meshes.Mesh(lines.vertices, lines.triangles, owners[lines.owner]), look


def relook(scene: Scene) -> Scene:
    """The same scene with the displays' current colours and flags (nothing rebuilt)."""
    from dataclasses import replace

    colours = _cat([d.colours for d in scene.displays], np.float32, 4)
    flags = _cat([d.flags() for d in scene.displays], np.float32)
    looks = [d.residue_look() for d in scene.displays]
    surface_colours = _surface_colours(scene.displays, scene.surface)
    return replace(
        scene,
        sphere_look=scenes.sphere_look(scene.geometry, colours, flags),
        cylinder_look=scenes.cylinder_look(scene.geometry, colours, flags),
        cartoon_look=meshes.look(
            scene.cartoon,
            _cat([c for c, _ in looks], np.float32, 4),
            _cat([f for _, f in looks], np.float32),
        ),
        surface_look=meshes.look(scene.surface, surface_colours, flags),
        pocket_look=meshes.look(
            scene.pocket, _surface_colours(scene.displays, scene.pocket, pocket=True), flags
        ),
    )


def _surface_colours(
    shown: tuple[Display, ...], mesh: meshes.Mesh, pocket: bool = False
) -> np.ndarray:
    """Per atom, each display's surface scheme (a pocket's: hydrophobicity) with its surface's
    alpha; nothing when there is no surface to colour."""
    if not len(mesh):
        return np.zeros((0, 4), np.float32)
    return _cat([d.surface_colours(pocket) for d in shown], np.float32, 4)


def _shifted(mesh: meshes.Mesh, offset: int) -> meshes.Mesh:
    return meshes.Mesh(mesh.vertices, mesh.triangles, mesh.owner + int(offset))


def _cat(arrays: list[np.ndarray], dtype, width: int | None = None) -> np.ndarray:
    if arrays:
        return np.concatenate(arrays).astype(dtype, copy=False)
    return np.zeros((0, width) if width else 0, dtype)
