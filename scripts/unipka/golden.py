#!/usr/bin/env python3
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

"""Free energies from the PyTorch model, for checking the numpy port (UNIPKA_PLAN.md, U3).

Runs in the training env (torch, Uni-Core, Uni-pKa's ``unimol`` package):

    python golden.py <checkpoint.pt> <dict dir> <in.json> <out.json>

*in.json* is a list of microstates ``{"atoms": [...], "charges": [...], "coords": [[[x, y,
z], ...], ...]}`` (conformers × atoms); *out.json* gets each one's free energy per conformer,
computed in fp32 with the fp16-rounded weights the archive ships, so the port must match it to
rounding.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

UNIPKA = Path.home() / "bench" / "unipka" / "src"


def main() -> None:
    checkpoint, dicts, src, dst = sys.argv[1:5]
    sys.path.insert(0, str(UNIPKA))
    from unicore import options, tasks, utils

    utils.import_user_module(argparse.Namespace(user_dir=str(UNIPKA / "unimol")))
    state = torch.load(checkpoint, map_location="cpu", weights_only=False)
    parser = options.get_training_parser()
    argv = [dicts, "--user-dir", str(UNIPKA / "unimol"), "--task", "mol_pka",
            "--arch", "unimol_pka", "--loss", "finetune_mse",
            "--classification-head-name", "chembl", "--num-classes", "1",
            "--split-mode", "predefine", "--dict-name", "dict.txt",
            "--charge-dict-name", "dict_charge.txt"]  # fmt: skip
    args = options.parse_args_and_arch(parser, argv)
    task = tasks.setup_task(args)
    model = task.build_model(args)
    weights = {k: v.half().float() for k, v in state["model"].items()}  # as the archive ships
    model.load_state_dict(weights, strict=False)
    model.eval().float()
    out = []
    for item in json.loads(Path(src).read_text()):
        atoms = [task.dictionary.index(a) for a in item["atoms"]]
        charges = [task.charge_dictionary.index(str(q)) for q in item["charges"]]
        coords = np.array(item["coords"], dtype=np.float32)
        c, n = coords.shape[0], coords.shape[1] + 2
        xyz = coords - coords.mean(axis=1, keepdims=True)
        xyz = np.concatenate([np.zeros((c, 1, 3)), xyz, np.zeros((c, 1, 3))], axis=1)
        tok = torch.tensor([task.dictionary.bos(), *atoms, task.dictionary.eos()]).repeat(c, 1)
        chg = torch.tensor([task.charge_dictionary.bos(), *charges,
                            task.charge_dictionary.eos()]).repeat(c, 1)  # fmt: skip
        coord = torch.tensor(xyz, dtype=torch.float32)
        dist = torch.cdist(coord, coord)
        edge = tok.view(c, n, 1) * len(task.dictionary) + tok.view(c, 1, n)
        with torch.no_grad():
            energy, _ = model((tok, chg, coord, dist, edge, [c]), classification_head_name="chembl",
                              features_only=True)  # fmt: skip
        out.append(energy.view(-1).tolist())
    Path(dst).write_text(json.dumps(out))
    print(len(out), "microstates")


if __name__ == "__main__":
    main()
