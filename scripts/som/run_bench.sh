#!/bin/bash
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
# Run the independent som-benchmark (MIT) under both label definitions, to settle whether ATTNSOM
# is the right basis for the tool (docs/SOM_PLAN.md, M1 step 6). Runs in WSL. Nothing here ships.
#
# Why: their committed numbers use --label-field som_all, the union of PRIMARY_SOM, SECONDARY_SOM
# and TERTIARY_SOM (4,163 positive atoms, rate 0.0994). Upstream ATTNSOM, and so every number of
# ours, uses PRIMARY_SOM alone (2,853 atoms, 0.0681). Their task has 46% more atoms counting as a
# hit, so their 0.684 top-1 and our 0.618 are not comparable. `label_field` reaches
# `build_atom_table`, which sets the labels used for both fitting and scoring, so re-running with
# som_primary puts their models on our definition exactly.
#
# `som_primary` runs first because it is the answer; `som_all` follows as the control, and
# reproducing their published figures within seed noise is what licenses trusting the first.
#
# Their paths.py hardcodes results/ under the repo root, so each run's output is copied out and the
# clone is restored to its committed state afterwards -- the repo stays pristine, as with train.sh.
#
# Excluded: phys_bde_gbm and gnn_bde need ALFABET bond dissociation energies, which drag in
# TensorFlow and which the plan rules out anyway; smartcyp needs an LGPL jar built separately.
set -euo pipefail

PYTHON=/home/hp/miniconda3/envs/som-bench/bin/python
REPO=/home/hp/bench/sombench/som-benchmark
OUT=/home/hp/bench/sombench/out
SEEDS="${BENCH_SEEDS:-10}"
JOBS="${BENCH_JOBS:-4}"
# fp_gbm is excluded: one seed took 96 minutes against 20 seconds for phys_gbm, and it is not a
# comparator that matters. The BDE models and smartcyp need inputs this environment does not build.
MODELS="${BENCH_MODELS:-prior basic_gbm phys_gbm gnn gnn_huckel}"

mkdir -p "$OUT"

# Keep their committed results, which the runs below would otherwise overwrite.
if [ ! -d "$OUT/baseline-committed" ]; then
    cp -r "$REPO/results" "$OUT/baseline-committed"
    echo "saved their committed results to $OUT/baseline-committed"
fi

cd "$REPO"

# Stages 1 and 2 are label-independent, so they run once.
"$PYTHON" -u scripts/01_download.py
"$PYTHON" -u scripts/02_curate.py

for FIELD in som_primary som_all; do
    echo
    echo "=== evaluating, label field = $FIELD, $SEEDS seeds, split mode random ==="

    # 03_evaluate.py *updates* rows in an existing metrics_per_run.csv keyed by
    # (model, mode, seed) rather than rewriting it, so without this the output keeps their
    # committed som_all rows for every model, mode and seed this run does not recompute -- which
    # silently mixes their numbers into ours. Clearing it means each run's file is only its own.
    rm -f "$REPO/results/metrics_per_run.csv" "$REPO/results/metrics_summary.csv"
    "$PYTHON" -u scripts/03_evaluate.py \
        --seeds "$SEEDS" \
        --modes random \
        --label-field "$FIELD" \
        --models $MODELS \
        --jobs "$JOBS" 2>&1 | tee "$OUT/evaluate-$FIELD.log"

    rm -rf "$OUT/$FIELD"
    cp -r "$REPO/results" "$OUT/$FIELD"
    echo "wrote $OUT/$FIELD"
done

# Leave the clone exactly as cloned, so a later `git pull` is clean and their committed numbers
# stay readable in place.
cd "$REPO" && git checkout -- results 2>/dev/null || true

echo BENCH_RUN_OK
