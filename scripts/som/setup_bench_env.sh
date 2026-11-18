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
# Build the environment for the independent som-benchmark (MIT), which we run to settle whether
# ATTNSOM is the right basis for the tool at all (docs/SOM_PLAN.md, M1 step 6). Runs in WSL,
# never in the app's interpreter.
#
# Their environment.yml does not solve on this host: it is pinned to the macOS/arm64 versions the
# committed results were produced with, and on linux-64 conda-forge the rdkit=2025.03.3 pin wants
# python_abi 3.9 against their python=3.12.16. Every dependency their pyproject.toml declares is
# on PyPI, so this creates a bare interpreter and lets pip resolve them, which skips the conda
# solver entirely.
#
# Version drift from their committed numbers is therefore expected and is handled by measurement
# rather than by pinning: run_bench.sh evaluates `som_all` first, and reproducing their published
# figures within seed noise is what licenses trusting the `som_primary` figures from the same env.
# That comparison is internal to one environment, which is what the question actually needs.
set -euo pipefail

CONDA=/home/hp/miniconda3/bin/conda
PYTHON=/home/hp/miniconda3/envs/som-bench/bin/python
REPO=/home/hp/bench/sombench/som-benchmark

# Remove any half-built env from a failed solve, so this is safe to re-run.
"$CONDA" env remove -y -n som-bench 2>/dev/null || true
"$CONDA" create -y -n som-bench python=3.12

"$PYTHON" -m pip install --upgrade pip

# Their declared dependencies, unpinned, plus pytest for `make test`. Two calls, because
# --index-url *replaces* PyPI rather than adding to it: asking for rdkit in the same call as the
# CPU torch index makes pip look for rdkit there and fail with "from versions: none".
"$PYTHON" -m pip install rdkit scikit-learn numpy pandas scipy matplotlib pytest

# CPU-only torch, as their environment.yml pins deliberately: the GNN here is small and a CUDA
# build would be a 2.5 GB download the benchmark cannot use.
"$PYTHON" -m pip install torch --index-url https://download.pytorch.org/whl/cpu

# Their package itself, editable, from the clone.
"$PYTHON" -m pip install -e "$REPO"

"$PYTHON" -c 'import rdkit, sklearn, numpy, pandas, scipy, torch, som; print("versions:", "rdkit", rdkit.__version__, "sklearn", sklearn.__version__, "numpy", numpy.__version__, "torch", torch.__version__)'

echo BENCH_SETUP_OK
