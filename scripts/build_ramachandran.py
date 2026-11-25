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
"""Write ``chemsuite/macromolecule/prepare/data/ramachandran.json``: φ/ψ log-probability grids
for the loop modeler (docs/LOOPS_PLAN.md L1), from the loop benchmark's training deposits.

Four classes, as Lovell et al. (Proteins 2003, 50, 437) split them: general, Gly, Pro, and
pre-Pro (the residue before a Pro, Gly and Pro aside). The residues counted are those of the
training half's entries (``loop_bench.training``) in ``samples/loops/pdb/`` (CC0 PDB data,
X-ray at 1.8 Å or better, one per 30% sequence cluster): every standard amino acid of the
first model with its whole backbone in one location, peptide-bonded to both neighbours, and
every backbone atom of the three at a B-factor under :data:`B_MAX` (Lovell's filter).

Each grid holds :data:`STEP`° bins, smoothed with a Gaussian of one bin (wrapping round) plus
a floor, as natural-log probabilities. ``favored`` and ``allowed`` are the levels that enclose
98% and 99.95% of the counted residues (Lovell's contours): a φ/ψ in a bin below ``allowed``
is an outlier.

    python scripts/build_ramachandran.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import gemmi
import numpy as np
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "review"))
import loop_bench as lb  # noqa: E402

OUT = ROOT / "chemsuite" / "macromolecule" / "prepare" / "data" / "ramachandran.json"
STEP = 10  # degrees
B_MAX = 30.0  # Å²
FLOOR = 1e-6  # probability added to every bin before normalizing
FAVORED, ALLOWED = 0.98, 0.9995
CLASSES = ("general", "gly", "pro", "prepro")


def angles(path: Path) -> list[tuple[str, float, float]]:
    """(class, φ, ψ) of each residue of *path* that passes the filters."""
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    out = []
    for chain in st[0]:
        polymer = chain.get_polymer()
        for k in range(1, len(polymer) - 1):
            prev, r, nxt = polymer[k - 1], polymer[k], polymer[k + 1]
            if any(x.name not in lb.AMINO_ACIDS for x in (prev, r, nxt)):
                continue
            atoms = [_one(x, n) for x in (prev, r, nxt) for n in lb.BACKBONE]
            if any(a is None for a in atoms) or any(a.b_iso >= B_MAX for a in atoms):
                continue
            if _peptide(prev, r) > 1.5 or _peptide(r, nxt) > 1.5:
                continue
            p = {(x, n): np.array(_one(x, n).pos.tolist()) for x in (prev, r, nxt)
                 for n in ("N", "CA", "C")}  # fmt: skip
            phi = lb._dihedral(p[prev, "C"], p[r, "N"], p[r, "CA"], p[r, "C"])
            psi = lb._dihedral(p[r, "N"], p[r, "CA"], p[r, "C"], p[nxt, "N"])
            out.append((lb.rama_class(r.name, nxt.name), phi, psi))
    return out


def _one(residue: gemmi.Residue, name: str) -> gemmi.Atom | None:
    """The atom *name* of *residue* when it has one location only."""
    found = [a for a in residue if a.name == name]
    return found[0] if len(found) == 1 and found[0].altloc == "\0" else None


def _peptide(a: gemmi.Residue, b: gemmi.Residue) -> float:
    return a.find_atom("C", "*").pos.dist(b.find_atom("N", "*").pos)


def grid(points: np.ndarray) -> dict:
    bins = 360 // STEP
    edges = np.linspace(-180, 180, bins + 1)
    counts, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=[edges, edges])
    density = gaussian_filter(counts, sigma=1.0, mode="wrap")
    p = density / density.sum() + FLOOR
    p /= p.sum()
    logp = np.log(p)
    # the levels that enclose 98% and 99.95% of the residues counted (by their own bins)
    index = (((points + 180.0) % 360.0) // STEP).astype(int)
    at = np.sort(logp[index[:, 0], index[:, 1]])
    favored = float(at[int(len(at) * (1 - FAVORED))])
    allowed = float(at[int(len(at) * (1 - ALLOWED))])
    return {"step": STEP, "n": len(points), "favored": round(favored, 4),
            "allowed": round(allowed, 4),
            "logp": np.round(logp, 4).tolist()}  # fmt: skip


def main() -> None:
    found: dict[str, list[tuple[float, float]]] = {c: [] for c in CLASSES}
    entries = [e for e, _entity in lb.picked() if lb.training(e) and (lb.PDB / f"{e}.cif").exists()]
    for entry in entries:
        for kind, phi, psi in angles(lb.PDB / f"{entry}.cif"):
            found[kind].append((phi, psi))
    classes = {c: grid(np.array(found[c])) for c in CLASSES}
    OUT.write_text(json.dumps({
        "source": "loop benchmark training half (scripts/review/loop_bench.py): "
                  f"{len(entries)} PDB entries, X-ray at 1.8 Å or better, 30% clusters; "
                  f"backbone B under {B_MAX} Å²; classes after Lovell et al. 2003",
        "classes": classes,
    }, separators=(",", ":")))  # fmt: skip
    for c in CLASSES:
        print(
            c, classes[c]["n"], "favored", classes[c]["favored"], "allowed", classes[c]["allowed"]
        )


if __name__ == "__main__":
    main()
