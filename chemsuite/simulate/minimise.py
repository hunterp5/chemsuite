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

"""A prepared structure minimised by the OpenMM worker, and what came of it.

:func:`minimise` charges the small molecules (once: a spec that has its charges keeps them),
runs the staged protocol of :class:`~.spec.MinimiseSpec` and measures the result
(:func:`measure`): how far the deposited heavy atoms moved (overall, in the pocket, the
ligand), the closest contact before and after, the residues that moved most and where the
largest force is left. Heavy atoms Prepare built (``Prepared.built``) are free to move and are
left out of every RMSD: no deposited position says where they belong.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import Cancelled, CancelToken

from ..protein.structure.model import Structure
from . import build, spec
from .result import Contact, MinimiseResult
from .runner import Minimised, Progress, Simulator
from .spec import ForceFieldSpec, MinimiseSpec, SystemSpec

POCKET = 5.0  # Å: heavy atoms this near a small molecule's are its pocket
CONTACT_SEARCH = 4.0  # Å: contacts are looked for this far out
MOVED_MOST = 5  # residues listed as moving most


@dataclass(frozen=True)
class Done:
    system: SystemSpec  # with its charges, to be kept with the document and reused
    result: MinimiseResult


def minimise(
    worker: Simulator,
    s: Structure,
    system: SystemSpec,
    forcefield: ForceFieldSpec | None = None,
    protocol: MinimiseSpec | None = None,
    *,
    built: Iterable[int] = (),
    model: int = 0,
    platform: str = "",
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Done:
    """*s*'s *model* minimised as *system* (from :func:`.build.system_spec`). *built* are the
    uids Prepare built, free to move. Cancelling stops the worker between iterations and
    raises :class:`Cancelled` carrying the :class:`Done` it got to, or None before it began."""
    forcefield = forcefield or ForceFieldSpec()
    protocol = protocol or MinimiseSpec()
    versions = tuple(sorted(worker.hello(cancel=cancel).versions.items()))
    charged = worker.charges(system, forcefield, cancel=cancel)
    start = build.positions(s, charged, model)
    free = free_atoms(charged, built)
    try:
        done = worker.minimise(charged, forcefield, start, protocol, platform=platform,
                               free=free.nonzero()[0].tolist(), progress=progress,
                               cancel=cancel)  # fmt: skip
    except Cancelled as stopped:
        partial = stopped.partial
        if not isinstance(partial, Minimised):
            raise
        raise Cancelled(Done(charged, measure(charged, start, partial, free, versions))) from None
    return Done(charged, measure(charged, start, done, free, versions))


def free_atoms(system: SystemSpec, built: Iterable[int]) -> np.ndarray:
    """Which of *system*'s atoms are heavy atoms Prepare built (a mask in the spec's order)."""
    wanted = np.fromiter((int(u) for u in built), dtype=np.int64)
    return np.isin(np.array(system.uid, dtype=np.int64), wanted) & (np.array(system.element) != 1)


def measure(
    system: SystemSpec,
    start: np.ndarray,
    done: Minimised,
    free: np.ndarray,
    versions: tuple[tuple[str, str], ...] = (),
) -> MinimiseResult:
    """The result of a minimisation of *system* from *start* (Å, the spec's order), with the
    *free* atoms (a mask) left out of the RMSDs."""
    end = done.xyz
    residue_of = np.array(system.residue_of())
    kinds = np.array([r.kind for r in system.residues])
    heavy = np.array(system.element) != 1
    held = heavy & ~free
    moved = np.linalg.norm(end - start, axis=1)
    small = held & (kinds[residue_of] == spec.MOLECULE)
    pocket = None
    if small.any():
        near = cKDTree(start[small]).query(start, distance_upper_bound=POCKET)[0] <= POCKET
        pocket = _rmsd(moved, held & near & ~small)
    worst = np.zeros(len(system.residues))
    np.maximum.at(worst, residue_of[held], moved[held])
    order = np.argsort(-worst, kind="stable")[:MOVED_MOST]
    labels = [r.label for r in system.residues]
    return MinimiseResult(
        digest=system.digest(),
        uid=system.uid,
        xyz=tuple((float(x), float(y), float(z)) for x, y, z in end),
        energy_before=done.energy_before,
        energy_after=done.energy_after,
        rmsd=_rmsd(moved, held) or 0.0,
        rmsd_pocket=pocket,
        rmsd_ligand=_rmsd(moved, small),
        contact_before=closest_contact(system, start),
        contact_after=closest_contact(system, end),
        max_force=done.max_force,
        max_force_atom=atom_label(system, done.max_force_atom),
        moved_most=tuple((labels[r], float(worst[r])) for r in order if worst[r] > 0),
        charges=tuple((labels[m.residue], m.charge_model, m.charges) for m in system.molecules),
        versions=versions,
        platform=done.platform,
        seconds=tuple((st.name, st.seconds) for st in done.stages),
    )


def _rmsd(moved: np.ndarray, mask: np.ndarray) -> float | None:
    return float(np.sqrt((moved[mask] ** 2).mean())) if mask.any() else None


def atom_label(system: SystemSpec, i: int) -> str:
    """``A/HIS 57/NE2``."""
    residue = system.residues[system.residue_of()[i]]
    return f"{residue.label}/{system.name[i]}"


def closest_contact(system: SystemSpec, xyz: np.ndarray) -> Contact | None:
    """The closest two heavy atoms of different residues, neither an ion, that are not bonded
    and share no bonded neighbour (an amide O and the next N sit 2.25 Å apart); None when no
    two come within :data:`CONTACT_SEARCH`."""
    residue_of = np.array(system.residue_of())
    kinds = np.array([r.kind for r in system.residues])
    atoms = np.flatnonzero((np.array(system.element) != 1) & (kinds[residue_of] != spec.ION))
    if len(atoms) < 2:
        return None
    neighbours: list[set[int]] = [set() for _ in range(system.n_atoms)]
    for i, j, _order in system.bonds:
        neighbours[i].add(j)
        neighbours[j].add(i)
    pairs = atoms[cKDTree(xyz[atoms]).query_pairs(CONTACT_SEARCH, output_type="ndarray")]
    pairs = pairs[residue_of[pairs[:, 0]] != residue_of[pairs[:, 1]]]
    distance = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
    for k in np.argsort(distance, kind="stable").tolist():
        i, j = int(pairs[k, 0]), int(pairs[k, 1])
        if j not in neighbours[i] and not neighbours[i] & neighbours[j]:
            return Contact(atom_label(system, i), atom_label(system, j), float(distance[k]))
    return None
