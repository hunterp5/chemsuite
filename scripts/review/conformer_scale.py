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
"""What hundreds of conformers per molecule cost, outside the UI (ROADMAP §7.0f).

The set is 1,000 drug-like molecules from the BindingDB SMILES that ``chemspace_bench.py``
caches (``samples/chemspace/all.smi``, local data, not committed), each given CONFORGE
conformers at 10, 100 and 300 per molecule and saved as a document under
``samples/conformers/scale/`` (``scale_<n>.csd``: ID, a 2D SMILES column and "SMILES 3D"), so
the same files open in the app for the by-hand checks. Steps::

    python scripts/review/conformer_scale.py data [--counts 10,100,300] [--rows 1000]
    python scripts/review/conformer_scale.py run [--counts 10,100,300]

``run`` prints one Markdown table per count: each step on the 3D column next to the same step
on the 2D column, and the process's memory (RSS) as it goes. The app's own timings (painting,
the Browser, applying results) are the ``-m bench`` tests in
``tests/app/table/test_bench_conformers.py``.
"""

from __future__ import annotations

import argparse
import ctypes
import gc
import random
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Lipinski

from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.analysis.structure_search import SCREEN_KIND, StructureQuery, search
from chemsuite.conformers.etkdg import ConformerSettings, generate_conformers
from chemsuite.core.depict import depiction_mol
from chemsuite.core.document import load_document, save_document
from chemsuite.core.jobs.parallel import process_pool
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn
from chemsuite.descriptors.compute import compute_descriptors
from chemsuite.io.sdf import write_sdf

ROOT = Path(__file__).resolve().parents[2]
SMILES = ROOT / "samples" / "chemspace" / "all.smi"
DATA = ROOT / "samples" / "conformers" / "scale"
TWO_D, THREE_D = "SMILES", "SMILES 3D"
SEED = 2026


def rss_mb() -> float:
    """The process's resident memory in MB (no psutil: ctypes on Windows, resource elsewhere)."""
    if sys.platform == "win32":
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        kernel = ctypes.WinDLL("kernel32")
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        query = kernel.K32GetProcessMemoryInfo
        query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        counters = Counters()
        counters.cb = ctypes.sizeof(Counters)
        query(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb)
        return counters.WorkingSetSize / 1e6
    import resource  # macOS reports bytes, Linux kilobytes; this is the peak, close enough

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 1e6 if sys.platform == "darwin" else peak / 1e3


def pick(rows: int) -> list[str]:
    """*rows* distinct drug-like, flexible molecules (20–45 heavy atoms, 3–10 rotatable bonds,
    neutral), the same ones every run."""
    texts = SMILES.read_text(encoding="utf-8").split("\n")
    random.Random(SEED).shuffle(texts)
    out: list[str] = []
    with rdBase.BlockLogs():
        for text in texts:
            mol = Chem.MolFromSmiles(text)
            if mol is None or "." in text or Chem.GetFormalCharge(mol):
                continue
            heavy, rotors = mol.GetNumHeavyAtoms(), Lipinski.NumRotatableBonds(mol)
            if 20 <= heavy <= 45 and 3 <= rotors <= 10 and Descriptors.MolWt(mol) < 600:
                out.append(text)
            if len(out) == rows:
                break
    return out


def path_for(count: int) -> Path:
    return DATA / f"scale_{count}.csd"


def make_data(counts: list[int], rows: int) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    smiles = pick(rows)
    with process_pool() as pool:
        for count in counts:
            if path_for(count).exists():
                print(f"{path_for(count).name}: already made")
                continue
            settings = ConformerSettings(count=count, method="conforge", energy_window=None)
            start = time.perf_counter()
            batch = generate_conformers(smiles, settings, executor=pool)
            seconds = time.perf_counter() - start
            table = Table()
            table.add_rows({"ID": [f"M{i:04d}" for i in range(len(smiles))], TWO_D: smiles},
                           undoable=False)  # fmt: skip
            table.convert_column(TWO_D, "structure")
            table.add_column(THREE_D, StructureColumn.from_parsed(batch.blobs, batch.displays))
            save_document(table, path_for(count))
            made = batch.counts[batch.counts > 0]
            print(f"{path_for(count).name}: {len(smiles)} molecules, {made.mean():.0f} conformers "
                  f"on average (max {made.max():.0f}), made in {seconds:.0f} s")  # fmt: skip


class Report:
    def __init__(self, count: int) -> None:
        self.lines = [
            f"### {count} conformers per molecule",
            "",
            "| step | 2D column | 3D column | RSS after (MB) |",
            "|---|---|---|---|",
        ]

    def row(self, step: str, two_d: str, three_d: str) -> None:
        gc.collect()
        self.lines.append(f"| {step} | {two_d} | {three_d} | {rss_mb():,.0f} |")

    def print(self) -> None:
        print("\n".join(self.lines) + "\n")


def timed(fn) -> tuple[float, object]:
    start = time.perf_counter()
    result = fn()
    return time.perf_counter() - start, result


def seconds(value: float) -> str:
    return f"{value * 1e3:,.0f} ms" if value < 1 else f"{value:,.2f} s"


def run_count(count: int, pool) -> None:
    path = path_for(count)
    report = Report(count)
    StructureColumn._cache.clear()
    gc.collect()
    report.row("start", "", "")
    took, (table, _extras) = timed(lambda: load_document(path))
    rows = np.arange(len(table))
    columns = {name: table.column(name) for name in (TWO_D, THREE_D)}
    sizes = {n: sum(len(b) for b in c.blobs(rows) if b) / 1e6 for n, c in columns.items()}
    report.row("open the document", "", f"{seconds(took)} ({path.stat().st_size / 1e6:,.0f} MB "
               f"file, {sizes[THREE_D]:,.0f} MB of 3D blobs)")  # fmt: skip
    with tempfile.TemporaryDirectory() as folder:
        took, _ = timed(lambda: save_document(table, Path(folder) / "copy.csd"))
    report.row("save the document", "", seconds(took))

    def sweep(name: str) -> str:
        StructureColumn._cache.clear()
        took, _ = timed(lambda: [columns[name].mol(int(s)) for s in rows])
        return f"{seconds(took)} ({took / len(rows) * 1e3:.2f} ms a row)"

    report.row("load every molecule (cache cold)", sweep(TWO_D), sweep(THREE_D))
    StructureColumn._cache.clear()
    gc.collect()
    report.row("cache cleared", "", "")

    def pictures(name: str) -> str:
        mols = [columns[name].mol(int(s)) for s in rows[:200]]
        took, _ = timed(lambda: [depiction_mol(m) for m in mols if m is not None])
        return f"{took / len(mols) * 1e3:.2f} ms a cell"

    report.row("2D layout for a cell picture (200 rows)", pictures(TWO_D), pictures(THREE_D))

    def inputs(name: str) -> list:
        return columns[name].worker_inputs(rows)

    def fingerprints(name: str) -> str:
        took, _ = timed(lambda: compute_fingerprints(inputs(name), SCREEN_KIND, executor=pool))
        return seconds(took)

    report.row("substructure screen fingerprints (pool)", fingerprints(TWO_D),
               fingerprints(THREE_D))  # fmt: skip

    def substructure(name: str) -> str:
        fps = compute_fingerprints(inputs(name), SCREEN_KIND, executor=pool)
        query = StructureQuery("c1ccccc1C(=O)N")
        took, found = timed(lambda: search(query, inputs(name), fps, executor=pool))
        return f"{seconds(took)} ({int(found.hits.sum())} hits)"

    report.row("substructure search, amide (pool)", substructure(TWO_D), substructure(THREE_D))

    def descriptors(name: str) -> str:
        took, _ = timed(lambda: compute_descriptors(inputs(name), ["MolWt", "TPSA"], executor=pool))
        return seconds(took)

    report.row("MolWt + TPSA (pool)", descriptors(TWO_D), descriptors(THREE_D))
    with tempfile.TemporaryDirectory() as folder:
        out = Path(folder) / "all.sdf"
        took, _ = timed(lambda: write_sdf(table, out, structure=THREE_D, all_conformers=True))
        report.row("export every conformer (SDF)", "",
                   f"{seconds(took)} ({out.stat().st_size / 1e6:,.0f} MB)")  # fmt: skip
    report.print()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("step", choices=("data", "run"))
    parser.add_argument("--counts", default="10,100,300")
    parser.add_argument("--rows", type=int, default=1000)
    args = parser.parse_args()
    counts = [int(c) for c in args.counts.split(",")]
    if args.step == "data":
        make_data(counts, args.rows)
        return
    with process_pool() as pool:
        for count in counts:
            run_count(count, pool)


if __name__ == "__main__":
    main()
