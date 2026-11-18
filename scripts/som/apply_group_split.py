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
"""Patch a working copy of upstream's main.py to add a molecule-grouped split.

docs/SOM_PLAN.md, M1. Upstream cross-validates over molecule-isoform *instances*, so 86.2% of
test instances share their molecule with the training split under a different isoform
(``scripts/som/leakage_check.py``). That inflates the reported metrics as estimates of
performance on a new molecule, which is what the tool will actually be asked for.

This adds ``--group_split``, which keeps every instance of one molecule in the same fold by
using ``StratifiedGroupKFold`` with the canonical SMILES as the group. Nothing else changes, so
``--group_split`` off reproduces upstream exactly.

The clone stays pristine: this copies it to <work> and patches the copy. Every replacement must
match exactly once or the script fails, so an upstream change cannot be silently mispatched.

Usage:
    apply_group_split.py --upstream ~/bench/attnsom/upstream/ATTNSOM \
                         --work ~/bench/attnsom/work/ATTNSOM
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

REPLACEMENTS: list[tuple[str, str]] = [
    (
        "from sklearn.model_selection import StratifiedKFold\n",
        "from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold\n",
    ),
    (
        "    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=args.seed)\n",
        "    # chemsuite (docs/SOM_PLAN.md M1): group every instance of a molecule into one fold,\n"
        "    # so a structure cannot be trained on for one isoform and tested on for another.\n"
        "    if args.group_split:\n"
        "        skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,\n"
        "                                   random_state=args.seed)\n"
        "        split_groups = [g.smiles for g in graphs]\n"
        "        print(f'[chemsuite] molecule-grouped split over '\n"
        "              f'{len(set(split_groups))} unique molecules')\n"
        "    else:\n"
        "        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=args.seed)\n"
        "        split_groups = None\n",
    ),
    (
        "    for fold, (trainval_idx, test_idx) in enumerate("
        "skf.split(np.arange(len(graphs)), y_strat), start=1):\n",
        "    for fold, (trainval_idx, test_idx) in enumerate(\n"
        "            skf.split(np.arange(len(graphs)), y_strat, split_groups), start=1):\n",
    ),
    (
        "    parser.add_argument('--seed', type=int, default=42)\n",
        "    parser.add_argument('--seed', type=int, default=42)\n"
        "    parser.add_argument('--group_split', default=False, action='store_true',\n"
        '                        help="chemsuite: keep all isoform instances of a molecule "\n'
        '                             "in the same CV fold (docs/SOM_PLAN.md M1)")\n',
    ),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--work", required=True, type=Path)
    args = parser.parse_args()

    if args.work.exists():
        shutil.rmtree(args.work)
    # The clone carries a .git directory; the working copy does not need it.
    shutil.copytree(args.upstream, args.work, ignore=shutil.ignore_patterns(".git"))

    target = args.work / "main.py"
    text = target.read_text()

    for old, new in REPLACEMENTS:
        found = text.count(old)
        if found != 1:
            raise SystemExit(
                f"expected exactly one occurrence of:\n{old!r}\nfound {found}. "
                "Upstream main.py has changed; update apply_group_split.py."
            )
        text = text.replace(old, new)

    target.write_text(text)
    print(f"patched {target} ({len(REPLACEMENTS)} replacements)")
    print("run it with --group_split for the honest numbers, without for upstream's")


if __name__ == "__main__":
    main()
