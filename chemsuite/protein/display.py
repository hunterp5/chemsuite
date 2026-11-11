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

:func:`compose` puts the visible displays together for the view: atoms and residues are
numbered one after the other (a display's ``atom_offset`` and ``residue_offset``), so a pick
on the combined scene maps back to its protein with :meth:`Scene.locate`.

Representation state round-trips through :meth:`Display.state` as run-length lists, small
enough for a saved document.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from . import cartoon as cartoons
from . import colours as schemes
from . import interactions as contacts
from . import mesh as meshes
from . import scene as scenes
from . import surface as surfaces
from .model import AMINO, NUCLEIC, Structure

#: how see-through a surface is drawn
SURFACE_ALPHA = 0.6
CARTOON_PICK_RADIUS = 1.2


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
    styles: np.ndarray = field(init=False)
    cartoon: np.ndarray = field(init=False)  # bool per residue
    surface: np.ndarray = field(init=False)  # bool per atom
    colours: np.ndarray = field(init=False)  # rgba per atom
    selected: np.ndarray = field(init=False)  # bool per atom
    ghost: np.ndarray = field(init=False)  # bool per atom
    query: np.ndarray = field(init=False)  # bool per atom: whose interactions are shown
    _cache: dict = field(default_factory=dict, init=False, repr=False)
    _versions: dict = field(init=False, repr=False)

    def __post_init__(self) -> None:
        s = self.structure
        polymer = np.isin(s.res_kind, (AMINO, NUCLEIC))
        self.cartoon = polymer.copy()
        self.styles = scenes.default_styles(s, self.active, cartoon=bool(polymer.any()))
        self.surface = np.zeros(s.n_atoms, bool)
        self.colours = schemes.colour(s, self.scheme)
        self.selected = np.zeros(s.n_atoms, bool)
        self.ghost = np.zeros(s.n_atoms, bool)
        self.query = np.zeros(s.n_atoms, bool)
        self._versions = {"geometry": 0, "cartoon": 0, "surface": 0, "trace": 0,
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
        self.touch("geometry", "cartoon", "surface", "interactions")
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
        self.touch("geometry", "cartoon", "surface", "trace", "interactions")

    def show_interactions(self, atoms: np.ndarray | None) -> None:
        """Show the interactions of these atoms with the rest (None: none)."""
        self.query = np.zeros(self.structure.n_atoms, bool) if atoms is None else atoms.copy()
        self.touch("interactions")

    # built parts, cached by version
    def _cached(self, part: str, build):
        key = (part, self._versions[part])
        if key not in self._cache:
            self._cache = {k: v for k, v in self._cache.items() if k[0] != part}
            self._cache[key] = build()
        return self._cache[key]

    def geometry(self) -> scenes.Geometry:
        return self._cached(
            "geometry", lambda: scenes.geometry(self.structure, self.styles, self.model)
        )

    def cartoon_mesh(self) -> meshes.Mesh:
        residues = np.flatnonzero(self.cartoon)
        return self._cached(
            "cartoon",
            lambda: cartoons.build(self.structure, residues, self.model, active=self.active),
        )

    def surface_mesh(self) -> meshes.Mesh:
        atoms = np.flatnonzero(self.surface)
        return self._cached("surface", lambda: surfaces.build(self.structure, atoms, self.model))

    def trace(self) -> tuple[np.ndarray, np.ndarray]:
        """Cartoon residues and their trace atoms (what picks the cartoon)."""
        residues, atoms = self._cached("trace", lambda: cartoons.trace(self.structure, self.active))
        shown = self.cartoon[residues]
        return residues[shown], atoms[shown]

    def interactions(self) -> contacts.Interactions:
        if not self.query.any():
            return contacts.Interactions.empty()
        return self._cached(
            "interactions",
            lambda: contacts.find(self.structure, self.query, self.active, self.model),
        )

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

    # saving
    def state(self) -> dict:
        return {
            "visible": self.visible,
            "scheme": self.scheme,
            "surface_scheme": self.surface_scheme,
            "altlocs": {str(k): v for k, v in self.altlocs.items()},
            "model": self.model,
            "styles": runs(self.styles),
            "cartoon": runs(self.cartoon.astype(np.uint8)),
            "surface": runs(self.surface.astype(np.uint8)),
            "interactions": runs(self.query.astype(np.uint8)),
        }

    def set_state(self, state: Mapping) -> None:
        n, r = self.structure.n_atoms, self.structure.n_residues
        self.visible = bool(state.get("visible", True))
        self.altlocs = {int(k): str(v) for k, v in dict(state.get("altlocs", {})).items()}
        model = state.get("model", 0)
        self.model = model if isinstance(model, int) and 0 <= model < self.n_models else 0
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
    surface_colours = _cat(
        [schemes.colour(d.structure, d.surface_scheme) for d in shown], np.float32, 4
    )
    surface_colours[:, 3] = SURFACE_ALPHA
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
    surface_colours = _cat(
        [schemes.colour(d.structure, d.surface_scheme) for d in scene.displays], np.float32, 4
    )
    surface_colours[:, 3] = SURFACE_ALPHA
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
    )


def _shifted(mesh: meshes.Mesh, offset: int) -> meshes.Mesh:
    return meshes.Mesh(mesh.vertices, mesh.triangles, mesh.owner + int(offset))


def _cat(arrays: list[np.ndarray], dtype, width: int | None = None) -> np.ndarray:
    if arrays:
        return np.concatenate(arrays).astype(dtype, copy=False)
    return np.zeros((0, width) if width else 0, dtype)
