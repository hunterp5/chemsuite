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

"""Sketcher review: RDKit's 2D layout against CoordGen for Clean Structure (SKETCHER_PLAN §4.6).

    python scripts/review/layout_review.py SDF [N]

Takes the first N (default 500) structures of SDF that RDKit reads, largest fragment only,
drops their coordinates and lays each out with ``rdDepictor.Compute2DCoords`` and with
``rdCoordGen.AddCoords``. For each it counts, in bond lengths (median bond = 1):

- overlaps: non-bonded atom pairs closer than 0.3;
- crossings: pairs of bonds with no shared atom whose segments intersect;
- bond spread: the largest |length - 1| of any bond;
- time per structure.

It then draws each structure, jitters its atoms and runs Clean, checking the isomeric SMILES
is unchanged.

Results print as Markdown for docs/sketcher/layout-review.md.
"""

from __future__ import annotations

import random
import sys
import time
from pathlib import Path
from statistics import mean, median

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from chemsuite.sketch import convert, layout
from chemsuite.sketch.edits import History

METHODS = {
    "RDKit Compute2DCoords": layout.rdkit_coords,
    "CoordGen": layout.coordgen,
    "Best of both (chemsuite Clean)": None,
}


def stereo_check(mols: list[Chem.Mol]) -> tuple[int, int, list[str]]:
    """Draw each structure, jitter its atoms, Clean, and compare isomeric SMILES."""
    rng = random.Random(7)
    kept, refused, changed = 0, 0, []
    for mol in mols:
        sketch = convert.from_mol(mol)
        for aid, atom in list(sketch.atoms.items()):
            sketch.atoms[aid] = atom.moved(atom.x + rng.uniform(-0.3, 0.3),
                                           atom.y + rng.uniform(-0.3, 0.3))  # fmt: skip
        try:
            before = layout.stereo_key(sketch)
        except convert.SketchError:
            continue
        history = History(sketch)
        change = history.change("Clean")
        problem = layout.clean(change)
        history.commit(change)
        if problem:
            refused += 1
        elif layout.stereo_key(history.sketch) == before:
            kept += 1
        else:
            changed.append(Chem.MolToSmiles(mol))
    return kept, refused, changed


def load(path: Path, limit: int) -> list[Chem.Mol]:
    chooser = rdMolStandardize.LargestFragmentChooser()
    out = []
    for mol in Chem.SDMolSupplier(str(path)):
        if mol is None:
            continue
        mol = chooser.choose(mol)
        if mol.GetNumAtoms() < 2:
            continue
        mol.RemoveAllConformers()
        out.append(mol)
        if len(out) >= limit:
            break
    return out


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    path = Path(sys.argv[1])
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    mols = load(path, limit)
    print(f"{len(mols)} structures from {path.name} (largest fragment), "
          f"{mean(m.GetNumAtoms() for m in mols):.1f} heavy atoms on average.\n")  # fmt: skip
    print("| Method | With overlaps | With crossings | With a bond off by >10% | Overlaps "
          "| Crossings | Bond spread, median / worst | ms / structure |")  # fmt: skip
    print("|---|---|---|---|---|---|---|---|")
    for name, method in METHODS.items():
        rows, seconds = [], 0.0
        for mol in mols:
            work = Chem.Mol(mol)
            start = time.perf_counter()
            if method is None:
                work = layout.best_layout(work)
            else:
                work.RemoveAllConformers()
                method(work)
            seconds += time.perf_counter() - start
            rows.append(layout.mol_score(work))
        print(
            f"| {name} | {sum(r['overlaps'] > 0 for r in rows)} "
            f"| {sum(r['crossings'] > 0 for r in rows)} "
            f"| {sum(r['spread'] > 0.1 for r in rows)} "
            f"| {sum(r['overlaps'] for r in rows)} | {sum(r['crossings'] for r in rows)} "
            f"| {median(r['spread'] for r in rows):.3f} / {max(r['spread'] for r in rows):.3f} "
            f"| {1000 * seconds / len(mols):.1f} |"
        )
    kept, refused, changed = stereo_check(mols)
    print(f"\nClean on each structure after jittering its atoms: {kept} kept the same isomeric "
          f"SMILES, {refused} refused, {len(changed)} changed.")  # fmt: skip
    for smiles in changed:
        print(f"- changed: `{smiles}`")


if __name__ == "__main__":
    main()
