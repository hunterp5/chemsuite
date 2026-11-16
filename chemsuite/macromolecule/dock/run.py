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
"""Docking many ligands: progress, cancel, and the ligands that finished kept.

Each ligand is one :class:`Docked` — its poses, or why it has none — in the order they were
given. Vina itself runs single-threaded (``--cpu 1``); the ligands run beside each other in a
pool, so a cancel stops the rest and raises
:class:`~chemsuite.core.jobs.Cancelled` carrying what was finished.
"""

from __future__ import annotations

import itertools
import threading
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken, Progress
from chemsuite.core.jobs.parallel import default_workers, ordered_map

from .box import Box
from .pdbqt import PdbqtError
from .vina import Options, Pose, VinaError, dock_one, score_only


@dataclass(frozen=True)
class Docked:
    """What came of one ligand."""

    index: int  # where it was in the input
    poses: tuple[Pose, ...] = ()
    error: str = ""

    @property
    def best(self) -> Pose | None:
        return self.poses[0] if self.poses else None


@dataclass(frozen=True)
class Job:
    """Everything one ligand needs, so a worker takes a single argument."""

    index: int
    mol: Chem.Mol
    receptor: str
    box: Box | None
    options: Options
    cache: Path | None = None
    rescore: bool = False
    cancel: CancelToken | None = field(default=None, compare=False)


def _one(job: Job) -> Docked:
    try:
        if job.rescore:
            affinity = score_only(job.receptor, job.mol, job.box, job.options, job.cancel,
                                  job.cache)  # fmt: skip
            return Docked(job.index, (Pose(job.mol, 1, affinity),))
        if job.box is None:
            raise VinaError("Docking needs a search box.")
        poses = dock_one(job.receptor, job.mol, job.box, job.options, job.cancel, job.cache)
        return Docked(job.index, tuple(poses))
    except (VinaError, PdbqtError) as err:
        return Docked(job.index, error=str(err))


def dock_many(
    receptor_pdbqt: str,
    mols: Sequence[Chem.Mol],
    box: Box,
    options: Options | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    cache: Path | None = None,
    workers: int | None = None,
) -> list[Docked]:
    """Dock every molecule of *mols* into *receptor_pdbqt* inside *box*, in order.

    Raises :class:`~chemsuite.core.jobs.Cancelled` with the ligands already finished as its
    ``partial``, so the caller can keep those rows.
    """
    return _many(receptor_pdbqt, mols, box, options, progress, cancel, cache, workers, False)


def rescore_many(
    receptor_pdbqt: str,
    mols: Sequence[Chem.Mol],
    box: Box | None = None,
    options: Options | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    cache: Path | None = None,
    workers: int | None = None,
) -> list[Docked]:
    """Score every pose of *mols* where it already sits: one :class:`Pose` each."""
    return _many(receptor_pdbqt, mols, box, options, progress, cancel, cache, workers, True)


def _many(
    receptor_pdbqt: str,
    mols: Sequence[Chem.Mol],
    box: Box | None,
    options: Options | None,
    progress: Progress | None,
    cancel: CancelToken | None,
    cache: Path | None,
    workers: int | None,
    rescore: bool,
) -> list[Docked]:
    settings = options or Options()
    jobs = [
        Job(i, mol, receptor_pdbqt, box, settings, cache, rescore, cancel)
        for i, mol in enumerate(mols)
    ]
    done: list[Docked] = []
    total = len(jobs)
    if not total:
        return done

    counter = itertools.count(1)
    lock = threading.Lock()

    def counted(_result: Docked) -> None:
        if progress is not None:
            with lock:
                progress(next(counter), total)

    # Vina is a child process, so threads are enough to run several at once
    with ThreadPoolExecutor(workers or default_workers()) as pool:
        results = ordered_map(_one, jobs, pool, cancel=cancel, finished=counted)
        try:
            for result in results:
                done.append(result)
        except Cancelled as stop:
            raise Cancelled(done) from stop
    return done
