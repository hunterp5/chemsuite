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
# Fine-tune Uni-pKa's model from Uni-Mol's MIT backbone, one cross-validation fold, on one GPU
# (docs/UNIPKA_PLAN.md, U0/U2). Runs in WSL, in the training env (~/bench/unipka/.venv), never
# in the app's interpreter. Uni-pKa's finetune_pka.sh, less the distributed launch:
#
#   train.sh <data dir> <task name> <fold> <save dir> [epochs] [nfolds]
#
# <data dir> holds <task name>.lmdb and Uni-pKa's dict.txt and dict_charge.txt.
set -euo pipefail
data_path=$1
task_name=$2
fold=$3
save_dir=$4
epoch=${5:-20}
nfolds=${6:-5}
here=~/bench/unipka
source "$here/.venv/bin/activate"
cd "$here/src"
export OMP_NUM_THREADS=1
unicore-train "$data_path" --task-name "$task_name" --user-dir ./unimol \
    --train-subset train --valid-subset valid \
    --conf-size 11 --nfolds "$nfolds" --fold "$fold" --cv-seed 42 \
    --num-workers 6 --ddp-backend=c10d \
    --dict-name dict.txt --charge-dict-name dict_charge.txt \
    --task mol_pka --loss finetune_mse --arch unimol_pka \
    --classification-head-name chembl --num-classes 1 \
    --optimizer adam --adam-betas '(0.9, 0.99)' --adam-eps 1e-6 --clip-norm 1.0 \
    --lr-scheduler polynomial_decay --lr 3e-4 --warmup-ratio 0.06 --max-epoch "$epoch" \
    --batch-size 16 --pooler-dropout 0.1 --update-freq 2 --seed 0 \
    --fp16 --fp16-init-scale 4 --fp16-scale-window 256 \
    --log-interval 100 --log-format simple \
    --finetune-from-model "$here/weights/backbone_for_pka.pt" \
    --validate-interval 1 --keep-last-epochs 1 --all-gather-list-size 102400 \
    --save-dir "$save_dir/fold_$fold" --tmp-save-dir "$save_dir/tmp_$fold" \
    --best-checkpoint-metric valid_rmse --patience 2000 \
    --only-polar -1 --split-mode cross_valid
