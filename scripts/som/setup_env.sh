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
# Build the site-of-metabolism training environment in WSL (docs/SOM_PLAN.md, M1). Runs in WSL,
# never in the app's interpreter: the app has no torch and never will.
#
# Upstream's envs/attnsom.yaml is a full conda export, and it does not solve on this host: its
# ffmpeg and torchvision pins hit strict repo priority, and cuda-version 13 conflicts with the
# driver's __cuda=12.6. Neither torchvision, torchaudio nor ffmpeg is imported by the model, and
# wandb is optional, so this installs only what ATTNSOM actually uses, at upstream's versions.
#
# torch 2.4.0 is upstream's pin and the floor for torch.nn.RMSNorm, which the model uses.
# Layout: ~/bench/attnsom/{upstream/ATTNSOM (the MIT clone), scripts, runs}.
set -euo pipefail

/home/hp/miniconda3/bin/conda create -y -n attnsom python=3.10

/home/hp/miniconda3/envs/attnsom/bin/python -m pip install --upgrade pip

/home/hp/miniconda3/envs/attnsom/bin/python -m pip install \
    torch==2.4.0 --index-url https://download.pytorch.org/whl/cu118

/home/hp/miniconda3/envs/attnsom/bin/python -m pip install \
    'https://data.pyg.org/whl/torch-2.4.0%2Bcu118/torch_scatter-2.1.2%2Bpt24cu118-cp310-cp310-linux_x86_64.whl'

/home/hp/miniconda3/envs/attnsom/bin/python -m pip install \
    torch-geometric==2.6.1 \
    rdkit==2025.3.6 \
    numpy==2.0.1 \
    scipy==1.15.3 \
    scikit-learn==1.7.2 \
    pandas==2.3.2 \
    tqdm==4.67.1 \
    matplotlib==3.10.0 \
    levenshtein==0.27.1 \
    rapidfuzz==3.14.1

echo SETUP_OK
