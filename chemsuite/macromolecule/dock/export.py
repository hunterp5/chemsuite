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

"""A docking run written out as files, to run AutoDock Vina at the command line.

What Vina needs is split the way its own command line splits it: the **config file**
(:data:`CONFIG`) carries the receptor, the search box and the conditions, and the ligands are
passed as an argument. Vina 1.2's ``--batch`` takes a *directory*, so every ligand written to
:data:`LIGANDS` is docked by one command::

    vina --config vina.conf --batch ligands --dir out

Vina writes one ``<name>_out.pdbqt`` per input file into ``--dir``, so the names here are what
the results come back under. :data:`MANIFEST` records which row and conformer each file came
from, so a later run can be read back into the document.

**One file per conformer.** Vina searches a ligand's rotatable bonds itself (Meeko writes the
torsion tree, :func:`.pdbqt.ligand`), so several conformers of one ligand mostly repeat the
same search. What they do buy is the geometry Vina holds fixed: the pucker of an ordinary
saturated ring, which is never sampled. The caller decides how many to send
(:class:`LigandInput`); this module only writes what it is given.

No Qt, no table and no subprocess: this writes files. :mod:`.vina` runs them.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

from chemsuite.core.jobs import CancelToken, check

from . import pdbqt
from .box import Box
from .vina import Options

#: what is written, under the folder the chemist chose
CONFIG = "vina.conf"
RECEPTOR = "receptor.pdbqt"
LIGANDS = "ligands"
OUT = "out"
MANIFEST = "ligands.csv"
RUN_SH = "run.sh"
RUN_BAT = "run.bat"
#: the manifest's columns
MANIFEST_COLUMNS = ("file", "row", "name", "conformer")
#: how much of a row's name a file name keeps
NAME_MAX = 48


@dataclass(frozen=True)
class LigandInput:
    """One ligand file to write: a molecule with 3D coordinates, and where it came from."""

    mol: Chem.Mol
    name: str  # what the row is called
    row: int  # the row's oid, for the manifest
    conformer: int = 1  # 1-based, when a row sends several
    of: int = 1  # how many conformers this row sends in all


@dataclass(frozen=True)
class Saved:
    """What :func:`save` wrote."""

    folder: Path
    files: tuple[str, ...] = ()  # every file written, relative to the folder
    ligands: int = 0  # ligand files written
    rows: int = 0  # how many distinct rows they came from
    failed: tuple[tuple[str, str], ...] = ()  # (ligand name, why) for those that could not be
    command: str = ""


def config_text(box: Box, options: Options, receptor: str | None = RECEPTOR) -> str:
    """The Vina config file: the receptor, the search box and the conditions.

    The ligands are left out on purpose — they are the command's argument. ``cpu`` is left out
    too, so Vina uses every core it finds, which is what one ligand at a time wants (chemsuite
    docks with ``--cpu 1`` only because it runs the ligands beside each other).
    """
    lines = ["# Written by chemsuite. Run from this folder:", f"#   {command_line()}", ""]
    if receptor:
        lines.append(f"receptor = {receptor}")
        lines.append("")
    centre = ("center_x", "center_y", "center_z")
    size = ("size_x", "size_y", "size_z")
    lines += [f"{key} = {value:.3f}" for key, value in zip(centre, box.centre, strict=True)]
    lines += [f"{key} = {value:.3f}" for key, value in zip(size, box.size, strict=True)]
    lines += [
        "",
        f"scoring = {options.scoring}",
        f"exhaustiveness = {options.exhaustiveness}",
        f"num_modes = {options.poses}",
        f"energy_range = {options.energy_range:g}",
        f"seed = {options.seed}",
        "",
    ]
    return "\n".join(lines)


def command_line() -> str:
    return f"vina --config {CONFIG} --batch {LIGANDS} --dir {OUT}"


def scripts() -> dict[str, str]:
    """``run.sh`` and ``run.bat``, each running :func:`command_line` from its own folder."""
    return {
        RUN_SH: (
            "#!/bin/sh\n"
            "# Written by chemsuite. Needs AutoDock Vina on your PATH.\n"
            "set -e\n"
            f'cd "$(dirname "$0")"\n'
            f"mkdir -p {OUT}\n"
            f"{command_line()}\n"
        ),
        RUN_BAT: (
            "@echo off\r\n"
            "rem Written by chemsuite. Needs AutoDock Vina on your PATH.\r\n"
            "cd /d %~dp0\r\n"
            f"if not exist {OUT} mkdir {OUT}\r\n"
            f"{command_line()}\r\n"
        ),
    }


def file_stem(index: int, ligand: LigandInput) -> str:
    """A ligand's file name without its suffix: ``0007_imatinib``, or ``…_c02`` where a row
    sends several conformers. Numbered in the order given, so no two rows collide and the
    files sort as the table does."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", ligand.name).strip("._-")[:NAME_MAX]
    stem = f"{index:04d}_{safe}" if safe else f"{index:04d}"
    if ligand.of > 1:
        stem += f"_c{ligand.conformer:02d}"
    return stem


def save(
    folder: Path | str,
    box: Box,
    options: Options,
    receptor_pdbqt: str | None = None,
    ligands: Sequence[LigandInput] = (),
    cancel: CancelToken | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> Saved:
    """Write a docking run into *folder*, which is made if it is not there.

    Always writes the config. *receptor_pdbqt* adds :data:`RECEPTOR`; *ligands* adds a file
    each under :data:`LIGANDS`, with the manifest and the run scripts. A ligand Meeko cannot
    type is recorded in :attr:`Saved.failed` and the rest are still written.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    (folder / CONFIG).write_text(
        config_text(box, options, RECEPTOR if receptor_pdbqt is not None else None),
        encoding="utf-8",
    )
    written.append(CONFIG)
    if receptor_pdbqt is not None:
        (folder / RECEPTOR).write_text(receptor_pdbqt, encoding="utf-8")
        written.append(RECEPTOR)
    failed: list[tuple[str, str]] = []
    rows: set[int] = set()
    ligand_files = 0
    if ligands:
        where = folder / LIGANDS
        where.mkdir(exist_ok=True)
        (folder / OUT).mkdir(exist_ok=True)
        manifest: list[tuple[str, int, str, int]] = []
        if progress is not None:
            progress(0, len(ligands))
        for k, ligand in enumerate(ligands, 1):
            check(cancel)
            name = f"{file_stem(k, ligand)}.pdbqt"
            try:
                text = pdbqt.ligand(ligand.mol)
            except pdbqt.PdbqtError as err:
                failed.append((ligand.name, str(err)))
                continue
            (where / name).write_text(text, encoding="utf-8")
            written.append(f"{LIGANDS}/{name}")
            ligand_files += 1
            manifest.append((name, ligand.row, ligand.name, ligand.conformer))
            rows.add(ligand.row)
            if progress is not None:
                progress(k, len(ligands))
        _write_manifest(folder / MANIFEST, manifest)
        written.append(MANIFEST)
        for name, text in scripts().items():
            (folder / name).write_text(text, encoding="utf-8", newline="")
            written.append(name)
    return Saved(folder, tuple(written), ligand_files, len(rows), tuple(failed), command_line())


def _write_manifest(path: Path, rows: Sequence[tuple[str, int, str, int]]) -> None:
    """Which row and conformer each ligand file came from, so the results can be read back."""
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(MANIFEST_COLUMNS)
        writer.writerows(rows)


def read_manifest(path: Path | str) -> list[dict[str, str]]:
    """A manifest back as rows (``file``, ``row``, ``name``, ``conformer``)."""
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return [dict(row) for row in csv.DictReader(stream)]
