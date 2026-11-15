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
"""How well conformer generators reproduce crystal (bioactive) ligand poses (ROADMAP §7.0b).

The set is Wang et al.'s Dataset I (J. Chem. Inf. Model. 2023, 63, 6525–6536,
doi:10.1021/acs.jcim.3c01519): 3,354 drug-like ligands from the PDB with good electron density.
Their supporting table (``ci3c01519_si_002.xlsx``) names each ligand ``<pdb>_<CCD>-<chain>-<seq>
_<n>``; the poses come from the PDB (CC0) through RCSB's ModelServer, bond orders included.
Their Table 1 gives the numbers to compare with (RDKit ETKDG with no force field and no
pruning: AMR-R 1.599 / 1.003 / 0.781 / 0.648 Å at 1 / 10 / 50 / 250 conformers).

Everything lives under ``samples/conformers/wang/`` (local data, not committed):
``ids.tsv``, ``ligands/<title>.sdf`` and ``out/<variant>.jsonl``.

Each generator starts from the crystal ligand's SMILES (stereo read from the pose) and makes at
most N conformers. A ligand's score at N is the heavy-atom RMSD, with symmetry, of its closest
conformer: AMR-R is the mean, COV-R the share within a threshold. Steps, each skipping work
already done::

    python scripts/review/conformer_bench.py ids <ci3c01519_si_002.xlsx>
    python scripts/review/conformer_bench.py fetch
    python scripts/review/conformer_bench.py run <variant> [--counts 1,10,50,250] [--jobs 8]
        [--limit N]
    python scripts/review/conformer_bench.py score [--limit N] > docs/conformers/bench.md
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign

from chemsuite.conformers.etkdg import ConformerSettings, generate

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "conformers" / "wang"
LIGANDS = DATA / "ligands"
OUT = DATA / "out"
MODEL_SERVER = "https://models.rcsb.org/v1/{pdb}/ligand?auth_asym_id={chain}&auth_seq_id={seq}"
COUNTS = (1, 10, 50, 250)
THRESHOLDS = (0.5, 1.0, 2.0)
MAX_MATCHES = 1000  # symmetry-equivalent atom maps tried per ligand

#: what is compared. "etkdg-raw" is Wang et al.'s setting for RDKit (no force field, no pruning).
VARIANTS: dict[str, ConformerSettings] = {
    "etkdg": ConformerSettings(),
    "etkdg-raw": ConformerSettings(prune_rms=0.0, optimize=None, energy_window=None),
    "etkdg-prune": ConformerSettings(optimize=None, energy_window=None),
    "etkdg-mmff": ConformerSettings(energy_window=None),
    "conforge": ConformerSettings(method="conforge"),
    "conforge-wide": ConformerSettings(method="conforge", energy_window=None),
}
#: Wang et al.'s Table 1, AMR-R (Å) at COUNTS
PUBLISHED = {
    "ConfGenX": (1.510, 0.952, 0.714, 0.597),
    "Conformator": (1.549, 0.969, 0.783, 0.669),
    "OMEGA": (1.480, 0.965, 0.739, 0.591),
    "RDKit ETKDG": (1.599, 1.003, 0.781, 0.648),
    "torsional diffusion": (1.742, 0.999, 0.731, 0.574),
}


# -- the set ------------------------------------------------------------------------------------
def ids(xlsx: Path) -> None:
    import openpyxl

    book = openpyxl.load_workbook(xlsx, read_only=True)
    rows = list(book.worksheets[0].iter_rows(values_only=True))[1:]
    DATA.mkdir(parents=True, exist_ok=True)
    lines = ["title\tpdb\tccd\tchain\tseq\trotatable"]
    for title, rotatable, *_ in rows:
        pdb, rest = title.split("_", 1)
        ccd, chain, seq = rest.rsplit("_", 1)[0].split("-")
        lines.append(f"{title}\t{pdb}\t{ccd}\t{chain}\t{seq}\t{rotatable}")
    (DATA / "ids.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(lines) - 1} ligands")


def _entries() -> list[dict[str, str]]:
    lines = (DATA / "ids.tsv").read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def _fetch_one(entry: dict[str, str]) -> str:
    """Download one pose with curl: on Windows, Python's certificate store lacks the issuer of
    models.rcsb.org's certificate (urllib cannot verify it); the system's TLS can."""
    path = LIGANDS / f"{entry['title']}.sdf"
    if path.exists():
        return ""
    url = MODEL_SERVER.format(pdb=entry["pdb"], chain=entry["chain"], seq=entry["seq"])
    part = path.with_suffix(".part")
    command = ["curl", "-sSf", "-A", "chemsuite", "--retry", "3", "--max-time", "60",
               "-o", str(part), url + "&encoding=sdf"]  # fmt: skip
    done = subprocess.run(command, capture_output=True, text=True)
    if done.returncode != 0:
        part.unlink(missing_ok=True)
        return f"{entry['title']}: {done.stderr.strip() or done.returncode}"
    if "$$$$" not in part.read_text(encoding="utf-8"):
        part.unlink()
        return f"{entry['title']}: no ligand returned"
    part.replace(path)
    return ""


def fetch(workers: int = 4) -> None:
    LIGANDS.mkdir(parents=True, exist_ok=True)
    entries = _entries()
    with ThreadPoolExecutor(workers) as pool:
        problems = [p for p in pool.map(_fetch_one, entries) if p]
    print(f"{len(entries) - len(problems)} of {len(entries)} ligands fetched")
    for problem in problems:
        print("  ", problem)


def crystal(title: str) -> Chem.Mol | None:
    """The crystal pose, heavy atoms, stereo from 3D; None if RDKit cannot read it."""
    with rdBase.BlockLogs():
        mols = [m for m in Chem.SDMolSupplier(str(LIGANDS / f"{title}.sdf")) if m is not None]
        if len(mols) != 1:  # several copies (alternate locations) or none: leave it out
            return None
        mol = Chem.RemoveHs(mols[0])
        Chem.AssignStereochemistryFrom3D(mol)
    return mol


# -- running ------------------------------------------------------------------------------------
def _maps(probe: Chem.Mol, ref: Chem.Mol) -> list[list[tuple[int, int]]]:
    matches = ref.GetSubstructMatches(probe, uniquify=False, useChirality=False,
                                      maxMatches=MAX_MATCHES)  # fmt: skip
    return [list(enumerate(m)) for m in matches]


def score_ligand(args: tuple[str, str, tuple[int, ...]]) -> list[dict]:
    """``{title, count, rmsd, made, seconds, error}`` for each count (a worker entry point)."""
    title, variant, counts = args
    ref = crystal(title)
    if ref is None:
        return [dict(title=title, count=n, error="unreadable crystal pose") for n in counts]
    smiles = Chem.MolToSmiles(ref)
    out = []
    for n in counts:
        start = time.perf_counter()
        result = generate(Chem.MolFromSmiles(smiles), replace(VARIANTS[variant], count=n))
        seconds = time.perf_counter() - start
        if result.mol is None:
            out.append(dict(title=title, count=n, seconds=seconds, error=result.error))
            continue
        probe = Chem.RemoveHs(result.mol)
        maps = _maps(probe, ref)
        if not maps:
            out.append(dict(title=title, count=n, seconds=seconds, error="no atom map"))
            continue
        best = min(
            rdMolAlign.GetBestRMS(probe, ref, c.GetId(), map=maps) for c in probe.GetConformers()
        )
        out.append(dict(title=title, count=n, rmsd=best, made=result.count, seconds=seconds))
    return out


def _done(variant: str) -> set[str]:
    path = OUT / f"{variant}.jsonl"
    if not path.exists():
        return set()
    return {json.loads(line)["title"] for line in path.read_text(encoding="utf-8").splitlines()}


def run(variant: str, counts: tuple[int, ...], jobs: int, limit: int | None) -> None:
    if variant not in VARIANTS:
        sys.exit(f"unknown variant {variant!r}; one of {', '.join(VARIANTS)}")
    OUT.mkdir(parents=True, exist_ok=True)
    titles = [e["title"] for e in _entries()][:limit]
    done = _done(variant)
    todo = [t for t in titles if t not in done and (LIGANDS / f"{t}.sdf").exists()]
    print(f"{variant}: {len(todo)} ligands to run, {len(done)} done", flush=True)
    with (
        ProcessPoolExecutor(jobs) as pool,
        (OUT / f"{variant}.jsonl").open("a", encoding="utf-8") as sink,
    ):
        futures = [pool.submit(score_ligand, (t, variant, counts)) for t in todo]
        for k, future in enumerate(as_completed(futures), 1):
            for row in future.result():
                sink.write(json.dumps(row) + "\n")
            sink.flush()
            if k % 100 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)}", flush=True)


# -- scoring ------------------------------------------------------------------------------------
def _rows(variant: str) -> Iterator[dict]:
    for line in (OUT / f"{variant}.jsonl").read_text(encoding="utf-8").splitlines():
        yield json.loads(line)


def score(limit: int | None) -> None:
    titles = [e["title"] for e in _entries()][:limit]
    keep = set(titles)
    variants = [v for v in VARIANTS if (OUT / f"{v}.jsonl").exists()]
    # only ligands every variant scored, so the columns compare like with like
    scored: dict[str, dict[tuple[str, int], dict]] = {}
    for v in variants:
        scored[v] = {(r["title"], r["count"]): r for r in _rows(v) if r["title"] in keep}
    common = {
        key
        for key in set.intersection(*(set(s) for s in scored.values()))
        if all("rmsd" in scored[v][key] for v in variants)
    }
    counts = sorted({n for _t, n in common})
    print(f"Ligands scored by every variant: {len({t for t, _n in common})} of {len(titles)}\n")
    print("AMR-R (Å, mean closest heavy-atom RMSD to the crystal pose; lower is better)\n")
    print("| method | " + " | ".join(str(n) for n in counts) + " |")
    print("|---" * (len(counts) + 1) + "|")
    for v in variants:
        cells = []
        for n in counts:
            values = [scored[v][(t, m)]["rmsd"] for t, m in common if m == n]
            cells.append(f"{np.mean(values):.3f}")
        print(f"| {v} | " + " | ".join(cells) + " |")
    for name, published in PUBLISHED.items():
        cells = [f"{published[COUNTS.index(n)]:.3f}" if n in COUNTS else "" for n in counts]
        print(f"| {name} (Wang et al.) | " + " | ".join(cells) + " |")
    for threshold in THRESHOLDS:
        print(f"\nCOV-R at {threshold} Å (% of ligands with a conformer within it)\n")
        print("| method | " + " | ".join(str(n) for n in counts) + " |")
        print("|---" * (len(counts) + 1) + "|")
        for v in variants:
            cells = []
            for n in counts:
                values = np.array([scored[v][(t, m)]["rmsd"] for t, m in common if m == n])
                cells.append(f"{100 * np.mean(values <= threshold):.1f}")
            print(f"| {v} | " + " | ".join(cells) + " |")
    print("\nSeconds per ligand (median / mean), conformers made (mean), failures\n")
    print("| method | count | median s | mean s | made | failed |")
    print("|---|---|---|---|---|---|")
    for v in variants:
        for n in counts:
            rows = [r for (t, m), r in scored[v].items() if m == n]
            seconds = [r["seconds"] for r in rows if "seconds" in r]
            made = [r["made"] for r in rows if "made" in r]
            failed = sum("error" in r for r in rows)
            print(f"| {v} | {n} | {np.median(seconds):.2f} | {np.mean(seconds):.2f} "
                  f"| {np.mean(made):.1f} | {failed} |")  # fmt: skip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    steps = parser.add_subparsers(dest="step", required=True)
    steps.add_parser("ids").add_argument("xlsx", type=Path)
    steps.add_parser("fetch").add_argument("--workers", type=int, default=4)
    runner = steps.add_parser("run")
    runner.add_argument("variant")
    runner.add_argument("--counts", default=",".join(map(str, COUNTS)))
    runner.add_argument("--jobs", type=int, default=8)
    runner.add_argument("--limit", type=int)
    steps.add_parser("score").add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.step == "ids":
        ids(args.xlsx)
    elif args.step == "fetch":
        fetch(args.workers)
    elif args.step == "run":
        counts = tuple(int(n) for n in args.counts.split(","))
        run(args.variant, counts, args.jobs, args.limit)
    else:
        score(args.limit)


if __name__ == "__main__":
    main()
