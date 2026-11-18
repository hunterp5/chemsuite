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
#
"""Summarize the som-benchmark runs and check them before they are believed.

docs/SOM_PLAN.md, M1 step 6. Reads the three result sets `run_bench.sh` leaves in
/home/hp/bench/sombench/out and prints one table, after asserting that each file holds only the
rows its own run produced.

That check is the point of this script. `03_evaluate.py` updates rows in an existing
metrics_per_run.csv keyed by (model, mode, seed) rather than rewriting it, so a file can silently
carry their committed som_all rows for anything a run did not recompute. A table built from such a
file compares their numbers against themselves. Every earlier reading of this benchmark in
SOM_PLAN.md was wrong for one reason or another, so nothing here is reported without its control:

* the control is `som_all` against their committed figures -- small drift means the environment is
  sound and the som_primary numbers from it can be trusted;
* the comparison is `som_primary`, which is the label definition ours uses (PRIMARY_SOM alone),
  against our own grouped ATTNSOM numbers, quoted from the plan.

Run in WSL with any python3. Nothing here ships.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

OUT = Path("/home/hp/bench/sombench/out")

# Ours, from docs/SOM_PLAN.md: grouped 10-fold CV and the 20-fold production run. Both are
# molecule-grouped and PRIMARY_SOM-only, which is what makes them comparable to som_primary.
OURS = {
    "ATTNSOM, grouped 10-fold": {"top1": 0.618, "top2": 0.743},
    "ATTNSOM, production 20-fold": {"top1": 0.606, "top2": 0.746},
}

LABELS = {
    "prior": "atom-type frequency",
    "smartcyp": "SMARTCyp (untrained)",
    "fp_gbm": "fingerprint + GBM",
    "basic_gbm": "atom/bond + GBM",
    "phys_gbm": "Huckel + GBM",
    "gnn": "message-passing GNN",
    "gnn_huckel": "GNN + Huckel",
    "phys_bde_gbm": "Huckel + BDE + GBM",
    "gnn_bde": "GNN + BDE",
}


def read(path: Path, mode: str = "random") -> dict[str, list[dict]]:
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["mode"] == mode]
    out: dict[str, list[dict]] = {}
    for row in rows:
        out.setdefault(row["model"], []).append(row)
    return out


def summarize(runs: list[dict], key: str) -> tuple[float, float, int]:
    values = [float(r[key]) for r in runs]
    spread = statistics.stdev(values) if len(values) > 1 else 0.0
    return statistics.fmean(values), spread, len(values)


def check_own_rows(name: str, path: Path, expect_seeds: int) -> list[str]:
    """A run's file must hold only that run's rows: one mode, the seeds asked for, no extras."""
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    modes = sorted({r["mode"] for r in rows})
    seeds = sorted({int(r["seed"]) for r in rows})
    models = sorted({r["model"] for r in rows})
    problems = []
    if modes != ["random"]:
        problems.append(f"{name}: expected only mode 'random', found {modes}")
    if seeds != list(range(expect_seeds)):
        problems.append(f"{name}: expected seeds 0..{expect_seeds - 1}, found {seeds}")
    stale = {"smartcyp", "phys_bde_gbm", "gnn_bde", "fp_gbm"} & set(models)
    if stale:
        problems.append(
            f"{name}: holds rows for {sorted(stale)}, which this run did not compute -- "
            "they are their committed som_all numbers leaking through"
        )
    return problems


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=10)
    args = parser.parse_args()

    primary_path = OUT / "som_primary" / "metrics_per_run.csv"
    all_path = OUT / "som_all" / "metrics_per_run.csv"
    base_path = OUT / "baseline-committed" / "metrics_per_run.csv"

    problems = check_own_rows("som_primary", primary_path, args.seeds)
    problems += check_own_rows("som_all", all_path, args.seeds)
    if problems:
        print("REFUSING to build a table: the result files are not clean.")
        for line in problems:
            print("  " + line)
        raise SystemExit(1)
    print(f"files are clean: one mode, seeds 0..{args.seeds - 1}, no un-recomputed models\n")

    primary = read(primary_path)
    everything = read(all_path)
    baseline = read(base_path)

    print("=== control: our som_all against their committed som_all (10-seed means) ===")
    print(f"{'model':<24} {'theirs':>8} {'ours':>8} {'drift':>8}")
    drifts = []
    for model in sorted(everything):
        if model not in baseline:
            continue
        theirs, _, _ = summarize(baseline[model], "top1")
        ours, _, _ = summarize(everything[model], "top1")
        drifts.append(abs(ours - theirs))
        print(f"{LABELS.get(model, model):<24} {theirs:>8.4f} {ours:>8.4f} {ours - theirs:>+8.4f}")
    worst = max(drifts, default=0.0)
    verdict = "sound" if worst < 0.03 else "NOT trustworthy"
    print(f"\nlargest drift {worst:.4f} -> the environment is {verdict}\n")

    print("=== the comparison: primary-only labels (PRIMARY_SOM), molecule-grouped ===")
    print(f"{'model':<28} {'top-1':>16} {'top-2':>16} {'seeds':>6}")
    ranked = []
    for model, runs in primary.items():
        top1, sd1, n = summarize(runs, "top1")
        top2, sd2, _ = summarize(runs, "top2")
        ranked.append((top1, model, sd1, top2, sd2, n))
    for top1, model, sd1, top2, sd2, n in sorted(ranked, reverse=True):
        name = LABELS.get(model, model)
        print(f"{name:<28} {top1:>8.3f} ± {sd1:<5.3f} {top2:>8.3f} ± {sd2:<5.3f} {n:>6}")
    for name, got in OURS.items():
        print(f"{name:<28} {got['top1']:>8.3f} {'':<7} {got['top2']:>8.3f}")

    print()
    print("=== what the label field alone is worth (10-seed means, same splits) ===")
    print(f"{'model':<24} {'som_all':>9} {'som_primary':>12} {'cost':>8}")
    for model in sorted(primary):
        if model not in everything:
            continue
        a, _, _ = summarize(everything[model], "top1")
        p, _, _ = summarize(primary[model], "top1")
        print(f"{LABELS.get(model, model):<24} {a:>9.4f} {p:>12.4f} {p - a:>+8.4f}")

    best = max(ranked)
    ours_best = max(v["top1"] for v in OURS.values())
    print()
    name = LABELS.get(best[1], best[1])
    print(f"best of theirs on our labels : {name} {best[0]:.3f} ± {best[2]:.3f}")
    print(f"ours                         : {ours_best:.3f}")
    gap = ours_best - best[0]
    print(f"gap (ours minus theirs)      : {gap:+.3f}, against a seed spread of ±{best[2]:.3f}")
    if abs(gap) <= best[2]:
        print("=> inside one standard deviation: a tie.")
        print("   Decide on grounds other than accuracy.")
    elif gap > 0:
        print("=> ATTNSOM is genuinely ahead. Keep it and proceed with the port.")
    else:
        print("=> their model is genuinely ahead. Reconsider the architecture (the pivot).")


if __name__ == "__main__":
    main()
