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

"""Which surroundings carry signal about a residue's pKa shift on PKAD-R (titration round 2).

Each scored PKAD-R residue gets heavy-atom features, PROPKA's kinds of determinant:

- ``bb``: backbone N–H donors (anions) or backbone O acceptors (cations) near the group;
- ``sc``: neutral side-chain donors or acceptors (Ser, Thr, Tyr, Asn, Gln, Trp, His, Cys);
- ``ion``: charged partners within H-bond reach (Lys/Arg for an anion, Asp/Glu for a cation);
- ``buried``: PROPKA's buried fraction (heavy atoms within 15.5 Å, 280 → 0, 560 → 1);
- ``coulomb``: a unit charge of the group's sign with the charged groups within 12 Å,
  83/r² kcal/mol (ε = 4r), positive when it repels;
- ``helix``: backbone N of a helix's first turn near an anion (the helix dipole).

Then, per residue type, the shift (measured − model pKa) is fitted on these by ridge
regression, scored by leave-one-family-out cross-validation on the training half only. It is
compared with no shift (PROPKA's model pKa, the null) and with a fitted offset per type (the
fair null when offsets are fitted). The test half is not looked at.

    python scripts/review/pkad_features.py           # features.csv, then the CV table
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pkad_review as pk

OUT = pk.DATA / "features.csv"
GROUP = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "HIS": ("ND1", "NE2"), "CYS": ("SG",),
         "TYR": ("OH",), "LYS": ("NZ",)}  # fmt: skip
ANIONS = {"ASP", "GLU", "CYS", "TYR"}
SC_DONORS = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("ASN", "ND2"), ("GLN", "NE2"),
             ("TRP", "NE1"), ("HIS", "ND1"), ("HIS", "NE2"), ("CYS", "SG")}  # fmt: skip
SC_ACCEPTORS = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("ASN", "OD1"), ("GLN", "OE1"),
                ("HIS", "ND1"), ("HIS", "NE2")}  # fmt: skip
CATIONS = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}
ANION_O = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
NAMES = ("bb", "sc", "ion", "buried", "coulomb", "helix")


def ramp(d: np.ndarray, near: float = 2.9, far: float = 3.5) -> np.ndarray:
    return np.clip((far - d) / (far - near), 0, 1)


def charges(s) -> list[tuple[np.ndarray, float, int]]:
    """(centre, charge, residue) of each charged group, His left neutral."""
    from chemsuite.macromolecule.structure.model import AMINO

    x = s.xyz()
    out = []
    for c in range(s.n_chains):
        first = True
        for r in range(s.chain_start[c], s.chain_start[c + 1]):
            if s.res_kind[r] != AMINO:
                continue
            span = s.atoms_of(r)
            at = {str(s.name[i]): i for i in range(span.start, span.stop)}
            name = str(s.res_name[r])
            if first and "N" in at:
                out.append((x[at["N"]], 1.0, r))
            first = False
            if "OXT" in at:
                out.append((x[at["OXT"]], -1.0, r))
            if name in ("ASP", "GLU") and all(a in at for a in GROUP[name]):
                out.append((x[[at[a] for a in GROUP[name]]].mean(axis=0), -1.0, r))
            elif name == "LYS" and "NZ" in at:
                out.append((x[at["NZ"]], 1.0, r))
            elif name == "ARG" and "CZ" in at:
                out.append((x[at["CZ"]], 1.0, r))
    return out


def helix_starts(s) -> np.ndarray:
    """Backbone N atoms of the first four residues of each helix."""
    from chemsuite.macromolecule.structure import secondary
    from chemsuite.macromolecule.structure.model import HELIX

    ss = secondary.assign(s)
    atoms = []
    for c in range(s.n_chains):
        run = 0
        for r in range(s.chain_start[c], s.chain_start[c + 1]):
            run = run + 1 if ss[r] == HELIX else 0
            if 1 <= run <= 4:
                span = s.atoms_of(r)
                atoms += [i for i in range(span.start, span.stop) if s.name[i] == "N"]
    return np.array(atoms, dtype=np.int64)


def features(s, r: int, groups, helix: np.ndarray) -> dict[str, float] | None:
    name = str(s.res_name[r])
    span = s.atoms_of(r)
    at = {str(s.name[i]): i for i in range(span.start, span.stop)}
    if not all(a in at for a in GROUP[name]):
        return None
    x = s.xyz()
    g = x[[at[a] for a in GROUP[name]]]
    other = np.ones(s.n_atoms, dtype=bool)
    other[span.start : span.stop] = False
    heavy = (s.element > 1) & other
    d = np.linalg.norm(x[:, None] - g[None], axis=-1).min(axis=1)  # each atom to the group
    names, resn = s.name, s.res_name[s.residue_of]
    anion = name in ANIONS
    bb = sc = ion = 0.0
    for i in np.flatnonzero(heavy & (d < 3.6)).tolist():
        key = (str(resn[i]), str(names[i]))
        w = float(ramp(d[i]))
        if anion:
            if names[i] == "N" and resn[i] != "PRO":
                bb += w
            elif key in SC_DONORS:
                sc += w
            elif key in CATIONS:
                ion += w
        else:
            if names[i] == "O":
                bb += w
            elif key in SC_ACCEPTORS:
                sc += w
            elif key in ANION_O:
                ion += w
    centre = g.mean(axis=0)
    from chemsuite.macromolecule.structure.model import WATER

    near = np.linalg.norm(x - centre, axis=1) < 15.5
    count = int((near & (s.element > 1) & (s.atom_kind != WATER)).sum())
    buried = float(np.clip((count - 280) / 280, 0, 1))
    sign = -1.0 if anion else 1.0
    coulomb = 0.0
    for c, q, owner in groups:
        if owner == r:
            continue
        dist = float(np.linalg.norm(c - centre))
        if dist < 12.0:
            coulomb += sign * q * 83.0 / max(dist, 2.5) ** 2
    hx = 0.0
    if anion and len(helix):
        hx = float(ramp(d[helix][other[helix]], 4.0, 6.0).sum())
    return {"bb": bb, "sc": sc, "ion": ion, "buried": buried, "coulomb": coulomb, "helix": hx}


def build() -> None:
    from chemsuite.macromolecule.prepare import repair
    from chemsuite.macromolecule.structure import io

    by_pdb: dict[str, list] = {}
    for e in pk.entries():
        if pk.scored(e):
            by_pdb.setdefault(e.pdb, []).append(e)
    rows = []
    for pdb, wanted in sorted(by_pdb.items()):
        path = pk.PDBS / f"{pdb}.pdb"
        if not path.is_file():
            continue
        s = repair.one_location(io.read(path))
        groups, helix = charges(s), helix_starts(s)
        index = {(str(s.chain_name[s.chain_of[r]]), int(s.res_seq[r]), str(s.res_name[r])): r
                 for r in range(s.n_residues)}  # fmt: skip
        for e in wanted:
            r = index.get(e.key)
            f = None if r is None else features(s, r, groups, helix)
            if f is None:
                continue
            rows.append({"pdb": pdb, "chain": e.chain, "number": e.number, "name": e.name,
                         "family": e.family, "training": int(pk.training(e)),
                         "shift": e.pka - pk.MODEL[e.name], **f})  # fmt: skip
        print(pdb, flush=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=list(rows[0]))
        out.writeheader()
        out.writerows(rows)
    print(len(rows), "rows written to", OUT)


def ridge(x: np.ndarray, y: np.ndarray, lam: float) -> np.ndarray:
    """Coefficients (offset first, not penalised)."""
    a = np.hstack([np.ones((len(x), 1)), x])
    pen = lam * np.eye(a.shape[1])
    pen[0, 0] = 0
    return np.linalg.solve(a.T @ a + pen, a.T @ y)


def cv(rows: list[dict], cols: tuple[str, ...], lam: float) -> np.ndarray:
    """Leave-one-family-out predictions of the shift (cols empty: the offset alone)."""
    y = np.array([r["shift"] for r in rows])
    x = np.array([[r[c] for c in cols] for r in rows]).reshape(len(rows), len(cols))
    fam = np.array([r["family"] for r in rows])
    pred = np.zeros(len(rows))
    for f in set(fam):
        out = fam == f
        if out.all():
            continue
        w = ridge(x[~out], y[~out], lam)
        pred[out] = w[0] + x[out] @ w[1:]
    return pred


def report() -> None:
    with OUT.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for r in rows:
        for k in ("shift", *NAMES):
            r[k] = float(r[k])
    print("Counts per type (training / test):")
    for t in pk.MODEL:
        n = [sum(r["name"] == t and int(r["training"]) == h for r in rows) for h in (1, 0)]
        print(f"  {t}: {n[0]} / {n[1]}")
    train = [r for r in rows if int(r["training"])]
    print("\nPer type: Pearson r of each feature with the shift (training half)")
    for t in pk.MODEL:
        sub = [r for r in train if r["name"] == t]
        y = np.array([r["shift"] for r in sub])
        cells = []
        for c in NAMES:
            v = np.array([r[c] for r in sub])
            cells.append(f"{c} {np.corrcoef(v, y)[0, 1]:+.2f}" if v.std() > 0 else f"{c} –")
        print(f"  {t} (n={len(sub)}): " + ", ".join(cells))
    print("\nLeave-one-family-out RMSE on the training half (shift in pH units)")
    sets = {"null (no shift)": None, "offset per type": (), "all features": NAMES,
            "bb+sc+ion+buried": ("bb", "sc", "ion", "buried"),
            "buried+coulomb": ("buried", "coulomb")}  # fmt: skip
    for lam in (1.0, 10.0):
        print(f"  ridge λ = {lam}")
        for label, cols in sets.items():
            errors = []
            per = {}
            for t in pk.MODEL:
                sub = [r for r in train if r["name"] == t]
                y = np.array([r["shift"] for r in sub])
                if cols is None:
                    p = np.zeros(len(sub))
                else:
                    use = tuple(c for c in cols if c != "helix" or t in ANIONS)
                    p = cv(sub, use, lam)
                errors.append(p - y)
                per[t] = float(np.sqrt(np.mean((p - y) ** 2)))
            e = np.concatenate(errors)
            detail = ", ".join(f"{t} {v:.2f}" for t, v in per.items())
            print(f"    {label:20s} {np.sqrt(np.mean(e**2)):.3f}   ({detail})")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    if "--report" not in sys.argv or not OUT.exists():
        build()
    report()
