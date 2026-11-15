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
"""AutoDock Vina: the bundled binary, one ligand at a time.

The binary ships in ``bin/<platform>/`` (Apache-2.0, see the ``NOTICE`` there); chemsuite runs
it as a program and reads its output, so nothing of Vina is linked into the app.
:func:`executable` finds the one for this machine, copying it to a cache folder and marking it
executable when the install dropped the permission bit (wheels and zip files often do).

:func:`dock_one` docks a ligand and returns its :class:`Pose`s, best first; :func:`score_only`
scores a pose that already sits in the pocket. Both run the binary with ``--cpu 1``: many
ligands at once is :mod:`.run`'s business, not Vina's own threads.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

from chemsuite.core.jobs import CancelToken

from . import pdbqt
from .box import Box

log = logging.getLogger(__name__)

BIN = Path(__file__).resolve().parent / "bin"
#: ``sys.platform`` and ``platform.machine()`` → the folder in ``bin`` that holds the binary
BINARIES = {
    ("win32", "amd64"): "win_amd64/vina.exe",
    ("win32", "x86_64"): "win_amd64/vina.exe",
    ("darwin", "arm64"): "macos_arm64/vina",
    ("linux", "x86_64"): "linux_x86_64/vina",
    ("linux", "amd64"): "linux_x86_64/vina",
}
VINA = "vina"
VINARDO = "vinardo"
SCORING = (VINA, VINARDO)
#: seconds a single ligand may take before it is given up on
TIMEOUT_S = 1800
_RESULT = re.compile(r"^REMARK VINA RESULT:\s+(\S+)\s+(\S+)\s+(\S+)", re.MULTILINE)
_SCORE = re.compile(r"^Estimated Free Energy of Binding\s*:\s*(\S+)", re.MULTILINE)


class VinaError(Exception):
    """Vina could not be run, or would not dock this ligand."""


@dataclass(frozen=True)
class Options:
    scoring: str = VINA
    exhaustiveness: int = 8
    poses: int = 9
    energy_range: float = 3.0
    seed: int = 42


@dataclass(frozen=True)
class Pose:
    """One docked pose: the molecule where Vina put it, and what Vina thought of it."""

    mol: Chem.Mol
    rank: int  # 1 is the best
    affinity: float  # kcal/mol
    rmsd_lb: float = 0.0  # Å from the best pose, as Vina reports it
    rmsd_ub: float = 0.0


def executable(cache: Path | None = None) -> Path:
    """The Vina binary for this machine, ready to run.

    When the file where it ships is not executable (an install that dropped the permission
    bit), it is copied into *cache* and marked there; without a *cache* folder the system's
    temporary folder is used.
    """
    import platform

    machine = platform.machine().lower()
    key = (sys.platform, "amd64" if machine in ("amd64", "x86_64") else machine)
    name = BINARIES.get(key)
    if name is None:
        raise VinaError(f"chemsuite has no AutoDock Vina binary for {sys.platform} {machine}.")
    path = BIN / name
    if not path.exists():
        raise VinaError(f"The AutoDock Vina binary is missing from the install ({path}).")
    if sys.platform == "win32" or os.access(path, os.X_OK):
        return path
    folder = cache or Path(tempfile.gettempdir()) / "chemsuite-vina"
    folder.mkdir(parents=True, exist_ok=True)
    copy = folder / path.name
    if not copy.exists() or copy.stat().st_size != path.stat().st_size:
        shutil.copy2(path, copy)
    copy.chmod(copy.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return copy


def version(cache: Path | None = None) -> str:
    """What the bundled binary says it is (``AutoDock Vina v1.2.7``)."""
    done = _run([str(executable(cache)), "--version"], None)
    return done.stdout.splitlines()[0].strip() if done.stdout.strip() else ""


def dock_one(
    receptor_pdbqt: str,
    mol: Chem.Mol,
    box: Box,
    options: Options | None = None,
    cancel: CancelToken | None = None,
    cache: Path | None = None,
) -> list[Pose]:
    """Dock *mol* into *receptor_pdbqt* inside *box*: its poses, best first."""
    options = options or Options()
    if options.scoring not in SCORING:
        raise VinaError(f"Unknown scoring function {options.scoring!r}.")
    ligand_pdbqt = pdbqt.ligand(mol)
    arguments = [
        f"--scoring={options.scoring}",
        *box.arguments(),
        f"--exhaustiveness={options.exhaustiveness}",
        f"--num_modes={options.poses}",
        f"--energy_range={options.energy_range}",
        f"--seed={options.seed}",
    ]
    text = _vina(receptor_pdbqt, ligand_pdbqt, arguments, cancel, cache)
    return _poses(text)


def score_only(
    receptor_pdbqt: str,
    mol: Chem.Mol,
    box: Box | None = None,
    options: Options | None = None,
    cancel: CancelToken | None = None,
    cache: Path | None = None,
) -> float:
    """What Vina's scoring function makes of *mol* exactly where it is (kcal/mol).

    A grid is needed even to score, so without a *box* Vina makes one around the pose itself
    (``--autobox``), which is what rescoring a pose wants.
    """
    options = options or Options()
    if options.scoring not in SCORING:
        raise VinaError(f"Unknown scoring function {options.scoring!r}.")
    ligand_pdbqt = pdbqt.ligand(mol)
    where = box.arguments() if box is not None else ["--autobox"]
    done = _vina(
        receptor_pdbqt,
        ligand_pdbqt,
        [f"--scoring={options.scoring}", "--score_only", *where],
        cancel,
        cache,
        out=False,
    )
    found = _SCORE.search(done)
    if found is None:
        raise VinaError("Vina reported no score for this pose.")
    return float(found.group(1))


# -- running the binary -------------------------------------------------------------------------
def _vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    arguments: list[str],
    cancel: CancelToken | None,
    cache: Path | None,
    out: bool = True,
) -> str:
    """Run the binary in a folder of its own: the output PDBQT, or what it printed."""
    program = str(executable(cache))
    with tempfile.TemporaryDirectory(prefix="chemsuite-dock-") as folder:
        work = Path(folder)
        (work / "receptor.pdbqt").write_text(receptor_pdbqt, encoding="utf-8")
        (work / "ligand.pdbqt").write_text(ligand_pdbqt, encoding="utf-8")
        command = [
            program,
            "--receptor=receptor.pdbqt",
            "--ligand=ligand.pdbqt",
            *arguments,
            "--cpu=1",
        ]
        if out:
            command.append("--out=out.pdbqt")
        done = _run(command, cancel, cwd=work)
        if done.returncode != 0:
            raise VinaError(_message(done))
        if not out:
            return done.stdout
        written = work / "out.pdbqt"
        if not written.exists():
            raise VinaError(_message(done) or "Vina wrote no poses.")
        return written.read_text(encoding="utf-8")


def _run(
    command: list[str], cancel: CancelToken | None, cwd: Path | None = None
) -> subprocess.CompletedProcess[str]:
    from chemsuite.core.jobs import Cancelled

    try:
        process = subprocess.Popen(
            command,
            cwd=None if cwd is None else str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            **_no_window_kwargs(),
        )
    except OSError as err:
        raise VinaError(f"Could not start AutoDock Vina: {err}") from err
    start = time.monotonic()
    while True:
        try:
            stdout, stderr = process.communicate(timeout=0.2)
            break
        except subprocess.TimeoutExpired:
            if cancel is not None and cancel.cancelled:
                _kill(process)
                raise Cancelled from None
            if time.monotonic() - start > TIMEOUT_S:
                _kill(process)
                raise VinaError(
                    f"AutoDock Vina took longer than {TIMEOUT_S // 60} minutes on one ligand."
                ) from None
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def _no_window_kwargs() -> dict:
    """Keep Windows from flashing a console window for each ligand."""
    if sys.platform != "win32":
        return {}
    return {"creationflags": subprocess.CREATE_NO_WINDOW}


def _kill(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    try:
        process.kill()
        process.wait(timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        log.warning("could not stop AutoDock Vina (pid %s)", process.pid)


def _message(done: subprocess.CompletedProcess[str]) -> str:
    for text in (done.stderr, done.stdout):
        lines = [line.strip() for line in (text or "").splitlines() if line.strip()]
        bad = [line for line in lines if "error" in line.lower() or "usage" in line.lower()]
        if bad:
            return bad[-1]
    return f"AutoDock Vina stopped with code {done.returncode}."


def _poses(text: str) -> list[Pose]:
    """The output PDBQT's poses, with the affinities from its ``REMARK VINA RESULT`` lines."""
    mols = pdbqt.poses(text)
    numbers = [tuple(float(v) for v in found.groups()) for found in _RESULT.finditer(text)]
    if len(numbers) < len(mols):
        raise VinaError("Vina's output is missing the scores of some poses.")
    return [
        Pose(mol, rank, affinity, lower, upper)
        for rank, (mol, (affinity, lower, upper)) in enumerate(zip(mols, numbers, strict=False), 1)
    ]
