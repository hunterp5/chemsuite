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
# Train the site-of-metabolism model (docs/SOM_PLAN.md, M1). Runs in WSL, in the training env
# built by setup_env.sh, never in the app's interpreter.
#
# Usage:  train.sh <run-name> [extra main.py flags...]
#   train.sh repro-zaretzki                      # 10-fold CV on Zaretzki, the paper's Tables 2-3
#   train.sh repro-az --dataset_dir cyp_dataset/az   # AZ-ExactSOM, the paper's Table 4
#   train.sh smoke --n_splits 2 --max_epochs 1   # a fast shape-and-plumbing check
#
# Defaults come from upstream's main.py and match the paper: hidden 256, 4 message-passing
# layers, batch 32, 50 epochs, AdamW lr 1e-4 wd 1e-4, focal loss gamma 1, pos_weight 1.5,
# lambda_attn 0.5, 10 folds with 5% of each training split held out for validation.
# Checkpoints are always saved: the weight converter and the golden file need them.
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "usage: train.sh <run-name> [extra main.py flags...]" >&2
    exit 2
fi

RUN_NAME="$1"
shift

PYTHON=/home/hp/miniconda3/envs/attnsom/bin/python
# SOM_UPSTREAM points this at the working copy patched by apply_group_split.py, so a grouped
# run and an upstream-faithful run use the same driver. Defaults to the pristine clone.
UPSTREAM="${SOM_UPSTREAM:-/home/hp/bench/attnsom/upstream/ATTNSOM}"
RUNS=/home/hp/bench/attnsom/runs/"$RUN_NAME"

mkdir -p "$RUNS"

# main.py imports its siblings and resolves --dataset_dir relative to the clone, so it runs
# from there. Results and checkpoints go outside the clone, so the repo stays pristine.
cd "$UPSTREAM"

# No exec: the pipeline through tee needs pipefail (set above) to report main.py's own failure.
# -u because stdout is block-buffered through tee, which otherwise hides every print (fold
# headers, metrics) until the run ends, while stderr warnings flush at once and look like a hang.
"$PYTHON" -u main.py \
    --result_dir "$RUNS/results" \
    --save_ckpt \
    --ckpt_dir "$RUNS/checkpoint" \
    "$@" 2>&1 | tee "$RUNS/train.log"
