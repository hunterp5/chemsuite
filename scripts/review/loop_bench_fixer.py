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
"""PDBFixer's loop building for ``loop_bench.py pdbfixer``: run in the benchmark's own
environment (OpenMM and PDBFixer, never the app's), not imported by chemsuite.

    python loop_bench_fixer.py INPUTS OUT NAME [NAME ...]

Each ``INPUTS/NAME.cif`` gets its internal gaps built (missing chain ends are left alone, as
chemsuite leaves them) and its missing heavy atoms added, written to ``OUT/NAME.cif`` with the
input's numbering, and ``OUT/NAME.json`` with the time or the failure.
"""

import json
import sys
import time
import traceback
from pathlib import Path

from openmm.app import PDBxFile
from pdbfixer import PDBFixer


def fix(source: Path, target: Path) -> float:
    began = time.perf_counter()
    fixer = PDBFixer(filename=str(source))
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues):
        chain, at = key
        if at == 0 or at == len(list(chains[chain].residues())):
            del fixer.missingResidues[key]
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    seconds = time.perf_counter() - began
    with open(target, "w") as handle:
        PDBxFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    return seconds


def main() -> None:
    inputs, out, names = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    for k, name in enumerate(names, 1):
        numbers = {"name": name}
        try:
            numbers["seconds"] = fix(inputs / f"{name}.cif", out / f"{name}.cif")
        except Exception:  # noqa: BLE001 - a benchmark notes every failure and goes on
            numbers["failed"] = traceback.format_exc()
        (out / f"{name}.json").write_text(json.dumps(numbers))
        print(f"[{k}/{len(names)}] {name}", "failed" if "failed" in numbers else "", flush=True)


if __name__ == "__main__":
    main()
