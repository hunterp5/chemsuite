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
"""Alignment A0: MSA accuracy and speed against MAFFT, FAMSA and MUSCLE (docs/ALIGN_PLAN.md).

Everything lives under ``samples/align/`` (local data, never committed: the sets state no
licence, so they are used for measuring only):

- ``download/bench.tar.gz`` and ``bench/bench1.0/<set>/{in,ref,info,qscore}``: drive5's
  collection (BAliBASE 3, its structural subset BALIS, PREFAB 4, SABRE, OXBENCH), FASTA with
  upper case marking the core columns that are scored;
- ``out/<set>/<method>/<id>``: each method's alignment, and ``out/<set>/<method>.json`` its
  wall time per set.

Steps, each skipping work already done::

    python scripts/review/align_bench.py fetch
    # the comparators, in the WSL env cs-align-bench (times measured there):
    wsl ~/miniconda3/envs/cs-align-bench/bin/python <repo>/scripts/review/align_bench.py \\
        run --set bali3 --method mafft-linsi
    python scripts/review/align_bench.py score --set bali3 > samples/align/bali3-scores.md

Scores are computed here, not by qscore: **Q** (the sum-of-pairs score on the core) is the
share of the reference's core residue pairs, both upper case and in one column, that the test
also puts in one column; **TC** is the share of the reference's core columns (two or more upper
case residues) whose core residues the test puts together in one column. The set's own
``qscore/`` files (Edgar's runs of older versions) are listed beside ours as a cross-check.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tarfile
import time
import urllib.request
import zlib
from collections import defaultdict
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2] / "samples" / "align"
BENCH = ROOT / "bench" / "bench1.0"
OUT = ROOT / "out"
URL = "https://drive5.com/bench/bench.tar.gz"
GAPS = frozenset("-.~")

#: comparator command lines; ``{inp}`` and ``{out}`` are filled in, run in cs-align-bench
COMMANDS = {
    "mafft-linsi": "mafft --localpair --maxiterate 1000 --quiet --thread 1 {inp} > {out}",
    "mafft-auto": "mafft --auto --quiet --thread 1 {inp} > {out}",
    "mafft-fftns2": "mafft --retree 2 --maxiterate 0 --quiet --thread 1 {inp} > {out}",
    "famsa": "famsa -t 1 {inp} {out} 2> /dev/null",
    "muscle5": "muscle -align {inp} -output {out} -threads 1 2> /dev/null",
}


def read_fasta(path: Path) -> dict[str, str]:
    """Names to sequences (gapped, case kept), in file order."""
    records: dict[str, list[str]] = {}
    name = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif line and name is not None:
            records[name].append(line)
    return {k: "".join(v) for k, v in records.items()}


def columns(aligned: dict[str, str]) -> list[list[tuple[str, int, bool]]]:
    """Each column as (name, residue index, upper case) for the residues in it."""
    length = max(len(s) for s in aligned.values())
    out: list[list[tuple[str, int, bool]]] = [[] for _ in range(length)]
    for name, seq in aligned.items():
        k = 0
        for i, ch in enumerate(seq):
            if ch in GAPS:
                continue
            out[i].append((name, k, ch.isupper()))
            k += 1
    return out


def placement(aligned: dict[str, str]) -> dict[tuple[str, int], int]:
    """(name, residue index) to the test column it sits in."""
    where = {}
    for col, residues in enumerate(columns(aligned)):
        for name, k, _ in residues:
            where[(name, k)] = col
    return where


def score_one(test: dict[str, str], ref: dict[str, str]) -> tuple[float, float]:
    """(Q, TC) of a test alignment against a reference with upper-case core columns."""
    where = placement(test)
    pairs = right = 0
    core_cols = right_cols = 0
    for col in columns(ref):
        core = [(n, k) for n, k, upper in col if upper]
        if len(core) < 2:
            continue
        core_cols += 1
        placed = {where.get(r) for r in core}
        if len(placed) == 1 and None not in placed:
            right_cols += 1
        for a, b in combinations(core, 2):
            pairs += 1
            if where.get(a) is not None and where.get(a) == where.get(b):
                right += 1
    return (right / pairs if pairs else 1.0), (right_cols / core_cols if core_cols else 1.0)


def ids(bench: str, half: str = "") -> list[str]:
    """The set's families; ``half`` "train" or "test" keeps half, fixed by a hash of the ID
    (the defaults are tuned on train and reported on test)."""
    every = sorted(p.name for p in (BENCH / bench / "ref").iterdir())
    if not half:
        return every
    return [f for f in every if (zlib.crc32(f.encode()) % 2 == 0) == (half == "train")]


def category(bench: str, family: str) -> str:
    if bench.startswith("bali3"):
        return family[:4]
    if bench.startswith("sabre"):
        return family.split("_")[0]
    return bench


# --- steps ---------------------------------------------------------------------------------


def fetch() -> None:
    packed = ROOT / "download" / "bench.tar.gz"
    if not packed.exists():
        packed.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(URL, packed)
    if not BENCH.exists():
        with tarfile.open(packed) as tar:
            tar.extractall(ROOT / "bench", filter="data")
    print(f"{BENCH}: {', '.join(sorted(p.name for p in BENCH.iterdir() if p.is_dir()))}")


def run(bench: str, method: str, only: list[str], half: str, label: str, options: dict) -> None:
    """One method over one set; run with the cs-align-bench interpreter (in WSL). chemsuite's
    own runs take ``options`` (``msa.align`` keywords) and a ``label`` naming them."""
    name = label or method
    folder = OUT / bench / name
    folder.mkdir(parents=True, exist_ok=True)
    times_path = OUT / bench / f"{name}.json"
    times = json.loads(times_path.read_text()) if times_path.exists() else {}
    env_bin = str(Path(sys.executable).parent)
    todo = [f for f in (only or ids(bench, half)) if f not in times]
    for n, family in enumerate(todo, 1):
        inp = BENCH / bench / "in" / family
        out = folder / family
        if method == "chemsuite":
            seconds = _run_chemsuite(inp, out, options)
        else:
            command = COMMANDS[method].format(inp=inp, out=out)
            start = time.perf_counter()
            # MAFFT leaves its scratch folders in the working directory: keep them out of the repo
            line = f'PATH="{env_bin}:$PATH" TMPDIR="{ROOT}" {command}'
            subprocess.run(["bash", "-c", line], check=True, cwd=ROOT)
            seconds = time.perf_counter() - start
        times[family] = round(seconds, 4)
        if n % 25 == 0 or n == len(todo):
            times_path.write_text(json.dumps(times, indent=0, sort_keys=True))
            print(f"{name} {bench}: {n}/{len(todo)}", flush=True)


def _run_chemsuite(inp: Path, out: Path, options: dict) -> float:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from chemsuite.macromolecule.align import msa  # built in A2

    seqs = read_fasta(inp)
    start = time.perf_counter()
    aligned = msa.align_strings(list(seqs.values()), **options)
    seconds = time.perf_counter() - start
    out.write_text("".join(f">{n}\n{s}\n" for n, s in zip(seqs, aligned, strict=True)))
    return seconds


def score(bench: str, methods: list[str], half: str) -> None:
    families = ids(bench, half)
    found = methods or sorted(p.name for p in (OUT / bench).iterdir() if p.is_dir())
    rows: dict[str, dict[str, list[tuple[float, float, float]]]] = {}
    for method in found:
        times_path = OUT / bench / f"{method}.json"
        times = json.loads(times_path.read_text()) if times_path.exists() else {}
        per: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
        for family in families:
            test = OUT / bench / method / family
            if not test.exists() or not test.stat().st_size:
                continue
            q, tc = score_one(read_fasta(test), read_fasta(BENCH / bench / "ref" / family))
            entry = (q, tc, times.get(family, float("nan")))
            per[category(bench, family)].append(entry)
            per["all"].append(entry)
        rows[method] = per
    published = _published(bench, set(families))
    cats = sorted({c for per in rows.values() for c in per} - {"all"}) + ["all"]
    print(f"## {bench}: mean Q / TC (n families), total seconds\n")
    print("| Method | " + " | ".join(cats) + " | seconds |")
    print("|---|" + "---|" * (len(cats) + 1))
    for method, per in rows.items():
        cells = []
        for c in cats:
            got = per.get(c, [])
            if got:
                q = sum(e[0] for e in got) / len(got)
                tc = sum(e[1] for e in got) / len(got)
                cells.append(f"{q:.3f} / {tc:.3f} ({len(got)})")
            else:
                cells.append("")
        total = sum(e[2] for e in per.get("all", []))
        print(f"| {method} | " + " | ".join(cells) + f" | {total:.1f} |")
    for method, per in published.items():
        cells = []
        for c in cats:
            got = per.get(c, [])
            cells.append(
                f"{sum(g[0] for g in got) / len(got):.3f} / {sum(g[1] for g in got) / len(got):.3f}"
                if got
                else ""
            )
        print(f"| {method} (qscore file) | " + " | ".join(cells) + " | |")


def _published(bench: str, keep: set[str]) -> dict[str, dict[str, list[tuple[float, float]]]]:
    """The set's own qscore results: method to category to (Q, TC) per family."""
    folder = BENCH / bench / "qscore"
    out: dict[str, dict[str, list[tuple[float, float]]]] = {}
    if not folder.exists():
        return out
    for path in sorted(folder.iterdir()):
        per: dict[str, list[tuple[float, float]]] = defaultdict(list)
        for line in path.read_text().splitlines():
            fields = dict(part.split("=", 1) for part in line.split(";") if "=" in part)
            if "Q" not in fields or "Ref" not in fields:
                continue
            family = fields["Ref"].rsplit("/", 1)[-1]
            if family not in keep:
                continue
            entry = (float(fields["Q"]), float(fields.get("TC", "nan")))
            per[category(bench, family)].append(entry)
            per["all"].append(entry)
        out[path.name] = per
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    steps = parser.add_subparsers(dest="step", required=True)
    steps.add_parser("fetch")
    p_run = steps.add_parser("run")
    p_run.add_argument("--set", required=True)
    p_run.add_argument("--method", required=True, choices=[*COMMANDS, "chemsuite"])
    p_run.add_argument("--only", default="")
    p_run.add_argument("--half", choices=["train", "test"], default="")
    p_run.add_argument("--label", default="")
    p_run.add_argument("--options", default="{}", help="msa.align keywords as JSON")
    p_score = steps.add_parser("score")
    p_score.add_argument("--set", required=True)
    p_score.add_argument("--method", action="append", default=[])
    p_score.add_argument("--half", choices=["train", "test"], default="")
    args = parser.parse_args()
    if args.step == "fetch":
        fetch()
    elif args.step == "run":
        only = [f for f in args.only.split(",") if f]
        run(args.set, args.method, only, args.half, args.label, json.loads(args.options))
    else:
        score(args.set, args.method, args.half)


if __name__ == "__main__":
    main()
