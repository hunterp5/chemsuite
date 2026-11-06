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

"""Clean Structure: a fresh 2D layout that keeps the drawing's place, turn and stereo
(SKETCHER_PLAN §4.6, docs/sketcher/layout-review.md).

- **Layout.** RDKit lays the molecule out twice, with CoordGen and with ``Compute2DCoords``,
  and the layout with fewer overlapping atoms, crossing bonds and stretched bonds is kept
  (that order). On 500 approved drugs CoordGen wins most often, the other covers several of its
  stretched macrocycles (``scripts/review/layout_review.py``).
- **Place and turn.** Each cleaned fragment is turned and moved (never mirrored) to fit where
  it was drawn as closely as possible (least-squares rigid fit, Kabsch), so Clean tidies a
  drawing without throwing it across the canvas.
- **Stereo.** The wedges as drawn are kept when they still give the same molecule on the new
  layout; otherwise RDKit chooses new wedges from the stereo read before the clean. A clean
  that would change the molecule (compared by isomeric SMILES) is refused.
"""

from __future__ import annotations

import math
from dataclasses import replace
from statistics import median

from rdkit import Chem, rdBase
from rdkit.Chem import rdCoordGen, rdDepictor

from chemsuite.sketch.convert import SketchError, atom_order, to_rdkit
from chemsuite.sketch.edits import Change
from chemsuite.sketch.model import Order, Sketch, Stereo

OVERLAP = 0.3  # bond lengths: closer non-bonded atoms overlap
STRETCHED = 0.1  # a bond more than 10% off the median length is stretched

_DIRS = {
    Chem.BondDir.BEGINWEDGE: Stereo.WEDGE,
    Chem.BondDir.BEGINDASH: Stereo.HASH,
    Chem.BondDir.UNKNOWN: Stereo.EITHER,
}

Point = tuple[float, float]


# -- scoring -----------------------------------------------------------------------------------


def _crosses(p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    d1, d2 = orient(p3, p4, p1), orient(p3, p4, p2)
    d3, d4 = orient(p1, p2, p3), orient(p1, p2, p4)
    return d1 * d2 < -1e-9 and d3 * d4 < -1e-9


def score(points: list[Point], bonds: list[tuple[int, int]]) -> dict[str, float]:
    """Layout faults, with lengths in median bond lengths: ``overlaps`` (non-bonded atom pairs
    closer than 0.3), ``crossings`` (bond pairs with no shared atom that cross), ``stretched``
    (bonds more than 10% off) and ``spread`` (the largest relative bond-length error)."""
    lengths = [math.dist(points[a], points[b]) for a, b in bonds]
    unit = median(lengths) if lengths and median(lengths) > 1e-9 else 1.0
    pts = [(x / unit, y / unit) for x, y in points]
    bonded = {frozenset(b) for b in bonds}
    overlaps = 0
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if frozenset((i, j)) not in bonded and math.dist(pts[i], pts[j]) < OVERLAP:
                overlaps += 1
    crossings = 0
    for k, (a, b) in enumerate(bonds):
        for c, d in bonds[k + 1 :]:
            if {a, b} & {c, d}:
                continue
            if _crosses(pts[a], pts[b], pts[c], pts[d]):
                crossings += 1
    errors = [abs(n / unit - 1) for n in lengths]
    return {
        "overlaps": overlaps,
        "crossings": crossings,
        "stretched": sum(e > STRETCHED for e in errors),
        "spread": max(errors, default=0.0),
    }


def mol_score(mol: Chem.Mol) -> dict[str, float]:
    conf = mol.GetConformer()
    points = [
        (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())
    ]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    return score(points, bonds)


def _rank(s: dict[str, float]) -> tuple:
    return (s["overlaps"], s["crossings"], s["stretched"], s["spread"])


# -- layout ------------------------------------------------------------------------------------


def coordgen(mol: Chem.Mol) -> None:
    rdCoordGen.AddCoords(mol, rdCoordGen.CoordGenParams())


def rdkit_coords(mol: Chem.Mol) -> None:
    """RDKit's own layout. chemsuite's depiction sets RDKit to prefer CoordGen globally, which
    would make ``Compute2DCoords`` run CoordGen too, so the preference is lifted meanwhile."""
    prefer = rdDepictor.GetPreferCoordGen()
    rdDepictor.SetPreferCoordGen(False)
    try:
        rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
    finally:
        rdDepictor.SetPreferCoordGen(prefer)


def layouts(mol: Chem.Mol) -> list[Chem.Mol]:
    """Copies of *mol* laid out by CoordGen and by RDKit, best first (fewest overlaps, then
    crossings, then stretched bonds)."""
    found = []
    for k, method in enumerate((coordgen, rdkit_coords)):
        work = Chem.Mol(mol)
        work.RemoveAllConformers()
        try:
            with rdBase.BlockLogs():
                method(work)
        except (RuntimeError, ValueError):
            continue
        if work.GetNumConformers():
            found.append((_rank(mol_score(work)), k, work))
    return [work for *_, work in sorted(found, key=lambda t: t[:2])]


def best_layout(mol: Chem.Mol) -> Chem.Mol:
    """A copy of *mol* with the better of the CoordGen and RDKit layouts as its conformer."""
    found = layouts(mol)
    if not found:
        raise SketchError("RDKit could not lay out this drawing.")
    return found[0]


def fit(new: list[Point], old: list[Point]) -> list[Point]:
    """*new* turned and moved (not mirrored) to lie as close as possible to *old* (Kabsch)."""
    n = len(new)
    if n == 0:
        return []
    ncx, ncy = sum(p[0] for p in new) / n, sum(p[1] for p in new) / n
    ocx, ocy = sum(p[0] for p in old) / n, sum(p[1] for p in old) / n
    sxx = sxy = 0.0
    for (nx, ny), (ox, oy) in zip(new, old, strict=True):
        ax, ay, bx, by = nx - ncx, ny - ncy, ox - ocx, oy - ocy
        sxx += ax * bx + ay * by  # cos term
        sxy += ax * by - ay * bx  # sin term
    angle = math.atan2(sxy, sxx) if (sxx or sxy) else 0.0
    c, s = math.cos(angle), math.sin(angle)
    return [
        (ocx + (x - ncx) * c - (y - ncy) * s, ocy + (x - ncx) * s + (y - ncy) * c) for x, y in new
    ]


# -- the edit ----------------------------------------------------------------------------------


def stereo_key(sketch: Sketch) -> str:
    """The molecule with its stereo as RDKit reads the drawing (isomeric SMILES)."""
    return Chem.MolToSmiles(to_rdkit(sketch))


def clean(change: Change, atoms=None) -> str:
    """Lay out the fragments holding *atoms* (all of them when None) afresh.

    Returns "" or why the drawing was left as it was.
    """
    sketch = change.sketch
    if sketch.is_empty():
        return ""
    try:
        mol = to_rdkit(sketch)
        before = Chem.MolToSmiles(mol)
    except SketchError:
        return "Fix the marked atoms first: RDKit cannot read this drawing."
    targets = set(sketch.atoms) if atoms is None else set(atoms) & set(sketch.atoms)
    fragments: list[list[int]] = []
    seen: set[int] = set()
    for aid in sorted(targets):
        if aid not in seen:
            frag = sketch.fragment(aid)
            seen |= frag
            if len(frag) > 1:
                fragments.append(sorted(frag))
    if not fragments:
        return ""
    for laid in layouts(mol):
        trial = sketch.copy()
        _place(trial, laid, fragments)
        if _same(trial, before):
            _apply(change, trial)
            return ""
        _rewedge(trial, mol)
        if _same(trial, before):
            _apply(change, trial)
            return ""
    return "Clean would change the stereochemistry here, so the drawing was left as it was."


def _place(trial: Sketch, laid: Chem.Mol, fragments: list[list[int]]) -> None:
    """Move *fragments* of *trial* to the layout, each fitted onto where it was drawn."""
    index = {aid: i for i, aid in enumerate(atom_order(trial))}
    conf = laid.GetConformer()
    lengths = [
        math.dist(_xy(conf, b.GetBeginAtomIdx()), _xy(conf, b.GetEndAtomIdx()))
        for b in laid.GetBonds()
    ]
    unit = median(lengths) if lengths and median(lengths) > 1e-9 else 1.0
    for ids in fragments:
        new = [(_xy(conf, index[a])[0] / unit, _xy(conf, index[a])[1] / unit) for a in ids]
        old = [trial.atoms[a].pos for a in ids]
        for aid, (x, y) in zip(ids, fit(new, old), strict=True):
            trial.atoms[aid] = trial.atoms[aid].moved(x, y)


def _rewedge(trial: Sketch, mol: Chem.Mol) -> None:
    """Replace the wedges of *trial* with RDKit's choice for *mol*'s stereo on trial's layout."""
    order = atom_order(trial)
    work = Chem.Mol(mol)
    conf = work.GetConformer()
    for i, aid in enumerate(order):
        atom = trial.atoms[aid]
        conf.SetAtomPosition(i, (atom.x, atom.y, 0.0))
    for bond in work.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            bond.SetBondDir(Chem.BondDir.NONE)
    with rdBase.BlockLogs():
        Chem.WedgeMolBonds(work, conf)
    bond_ids = sorted(trial.bonds)
    for i, bond in enumerate(work.GetBonds()):
        bid = bond_ids[i]
        old = trial.bonds[bid]
        if old.order is not Order.SINGLE or old.stereo is Stereo.EITHER:
            continue
        stereo = _DIRS.get(bond.GetBondDir(), Stereo.NONE)
        new = replace(old, stereo=stereo)
        if stereo is not Stereo.NONE and order[bond.GetBeginAtomIdx()] != old.begin:
            new = replace(new, begin=old.end, end=old.begin)
        trial.bonds[bid] = new


def _apply(change: Change, trial: Sketch) -> None:
    for aid, atom in trial.atoms.items():
        change.set_atom(aid, atom)
    for bid, bond in trial.bonds.items():
        change.set_bond(bid, bond)


def _same(sketch: Sketch, before: str) -> bool:
    try:
        return stereo_key(sketch) == before
    except SketchError:
        return False


def _xy(conf, i: int) -> Point:
    p = conf.GetAtomPosition(i)
    return (p.x, p.y)
