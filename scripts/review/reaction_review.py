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

"""Reaction review (REACTIONS_PLAN R1): reaction SMILES -> drawing -> RXN -> back.

    python scripts/review/reaction_review.py USPTO_CSV [LIMIT]

USPTO_CSV is ``data_processed.csv`` of github.com/connorcoley/retrosim (MIT): the USPTO-50k
set of Schneider et al. (J. Chem. Inf. Model. 2016, 56, 2336), 50,016 atom-mapped reactions
drawn from Lowe's CC0 extraction of US patents. Each ``rxn_smiles`` is read by RDKit (the
reference), drawn by ``chemsuite.sketch.reaction.from_text``, written back as reaction SMILES
and as an RXN file, and the RXN file read again. A reaction passes when every result equals
RDKit's own reaction SMILES of the input: same components on each side, in order, with the
same atom maps and stereo. Results print as Markdown for docs/sketcher/reaction-review.md.
"""

from __future__ import annotations

import csv
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

from rdkit import RDLogger
from rdkit.Chem import rdChemReactions

from chemsuite.sketch import convert, reaction


def side_sets(smiles: str) -> list[list[str]]:
    return [sorted(part.split(".")) for part in smiles.split(">")]


def check(text: str) -> tuple[str, str]:
    """(result, detail) for one reaction SMILES."""
    RDLogger.DisableLog("rdApp.*")
    want = rdChemReactions.ReactionToSmiles(
        rdChemReactions.ReactionFromSmarts(text, useSmiles=True)
    )
    try:
        sketch = reaction.from_text(text)
        drawn = reaction.to_smiles(sketch)
        again = reaction.to_smiles(reaction.from_text(reaction.to_rxn_block(sketch)))
    except (reaction.ReactionError, convert.SketchError, ValueError, RuntimeError) as e:
        return f"error: {type(e).__name__}", text
    if drawn == want and again == want:
        return "same", ""
    if side_sets(drawn) == side_sets(want) and side_sets(again) == side_sets(want):
        return "same components, other order", text
    return "differs", f"{text}\n    drawn: {drawn}\n    rxn:   {again}"


def main() -> None:
    path, limit = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else None
    with open(path, newline="", encoding="utf-8") as handle:
        texts = [row["rxn_smiles"] for row in csv.DictReader(handle)][:limit]
    start = time.time()
    tally: Counter[str] = Counter()
    examples: dict[str, list[str]] = {}
    with ProcessPoolExecutor() as pool:
        for kind, detail in pool.map(check, texts, chunksize=200):
            tally[kind] += 1
            if detail:
                examples.setdefault(kind, []).append(detail)
    print(f"{len(texts)} reactions in {time.time() - start:.0f} s\n")
    print("| Result | Reactions |\n|---|---|")
    for kind, count in tally.most_common():
        print(f"| {kind} | {count} |")
    for kind, items in examples.items():
        print(f"\n### {kind} ({len(items)})\n")
        for item in items:
            print(f"- `{item}`")


if __name__ == "__main__":
    main()
