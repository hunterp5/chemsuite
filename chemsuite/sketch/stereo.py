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

"""CIP stereo labels (R/S, E/Z) for the drawing, and edits that set them.

Labels come from RDKit's new CIP labeller (``rdCIPLabeler``, the Hanson et al. 2018 rules) run
on the drawing exactly as RDKit reads it: wedges define tetrahedral centres and 2D geometry
defines double bonds, as in any molfile (SKETCHER_PLAN §4.4). ``r`` / ``s`` are
pseudo-asymmetric centres. Centres and double bonds that could carry stereo but are not drawn
with it are reported as *undefined*.

Setting a label edits the drawing, never a stored flag:

* **R ↔ S** at a centre with wedges swaps each wedge and hash that starts at it; at an
  undefined centre a wedge is added on its most suitable bond and turned to a hash if that
  gives the other label. The result is checked with RDKit before it is kept.
* **E ↔ Z** reflects the smaller side of the double bond across the bond's axis (only possible
  when the bond is not in a ring); wedges in that side swap so its centres keep their labels.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler

from chemsuite.sketch import edits
from chemsuite.sketch.convert import SketchError, atom_order, to_rdkit
from chemsuite.sketch.edits import Change
from chemsuite.sketch.model import Order, Sketch, Stereo

_SWAP = {Stereo.WEDGE: Stereo.HASH, Stereo.HASH: Stereo.WEDGE}
_OPPOSITE = {"R": "S", "S": "R", "r": "s", "s": "r", "E": "Z", "Z": "E"}


@dataclass(frozen=True)
class StereoLabels:
    atoms: dict[int, str] = field(default_factory=dict)  # atom id -> R, S, r, s
    bonds: dict[int, str] = field(default_factory=dict)  # bond id -> E, Z
    undefined_atoms: frozenset[int] = frozenset()
    undefined_bonds: frozenset[int] = frozenset()
    #: double bonds drawn as "either" (crossed): E/Z deliberately unknown
    unknown_bonds: frozenset[int] = frozenset()


def bond_order(sketch: Sketch) -> list[int]:
    """Sketch bond ids in molfile order (bond *i* in RDKit is ``bond_order(sketch)[i]``)."""
    return sorted(sketch.bonds)


def labels(sketch: Sketch) -> StereoLabels:
    """The CIP labels RDKit gives the drawing (empty when RDKit rejects the drawing)."""
    if sketch.is_empty():
        return StereoLabels()
    try:
        mol = to_rdkit(sketch)
    except SketchError:
        return StereoLabels()
    atoms, bonds = atom_order(sketch), bond_order(sketch)
    with rdBase.BlockLogs():
        try:
            rdCIPLabeler.AssignCIPLabels(mol)
        except (RuntimeError, ValueError):  # e.g. too complex for the labeller
            return StereoLabels()
        potential = Chem.FindPotentialStereo(mol)
    atom_labels = {
        atoms[a.GetIdx()]: a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")
    }
    bond_labels = {
        bonds[b.GetIdx()]: b.GetProp("_CIPCode") for b in mol.GetBonds() if b.HasProp("_CIPCode")
    }
    unspecified = Chem.StereoSpecified.Unspecified
    open_atoms = {
        atoms[s.centeredOn]
        for s in potential
        if s.type == Chem.StereoType.Atom_Tetrahedral and s.specified == unspecified
    }
    open_bonds = {
        bonds[s.centeredOn]
        for s in potential
        if s.type == Chem.StereoType.Bond_Double and s.specified == unspecified
    }
    unknown = {
        bonds[s.centeredOn]
        for s in potential
        if s.type == Chem.StereoType.Bond_Double and s.specified == Chem.StereoSpecified.Unknown
    }
    return StereoLabels(
        atom_labels,
        bond_labels,
        frozenset(open_atoms - set(atom_labels)),
        frozenset(open_bonds - set(bond_labels)),
        frozenset(unknown - set(bond_labels)),
    )


# -- R / S ------------------------------------------------------------------------------------
#: How far (degrees) one bond at a centre is turned to test whether its reading is robust.
NUDGE = 10.0

NO_CENTRE = "This wedge defines nothing: the atom at its narrow end is not a stereocentre."
WRONG_END = (
    "This wedge points the wrong way: its narrow end must be at the stereocentre, which is "
    "the atom at its wide end."
)
CONTRADICTORY = "The wedges at this atom contradict each other, so its configuration is unknown."
AMBIGUOUS = (
    "Ambiguous: turning one bond at this centre slightly would change its configuration. "
    "Spread its bonds more evenly."
)


def wedge_problems(sketch: Sketch) -> dict[int, str]:
    """Wedges and hashes that do not say what they seem to, as bond id -> reason.

    Judged by RDKit reading the drawing as a molfile (SKETCHER_PLAN §4.4):

    - a wedge whose narrow end RDKit does not read as a stereocentre (or atropisomer axis)
      defines nothing; when its wide end is a possible stereocentre it is drawn backwards;
    - wedges at a possible stereocentre that RDKit cannot read a configuration from contradict
      each other (for example a wedge and a hash on opposite bonds);
    - a centre whose configuration changes when any one of its bonds is turned by 10 degrees
      is drawn ambiguously (three bonds in a half plane, two bonds on top of each other).

    Wavy bonds ("either") are never a problem.
    """
    wedged = {
        bid: b
        for bid, b in sketch.bonds.items()
        if b.order is Order.SINGLE and b.stereo in (Stereo.WEDGE, Stereo.HASH)
    }
    if not wedged:
        return {}
    try:
        mol = to_rdkit(sketch)
    except SketchError:
        return {}
    order = atom_order(sketch)
    index = {aid: i for i, aid in enumerate(order)}
    with rdBase.BlockLogs():
        potential = {
            s.centeredOn
            for s in Chem.FindPotentialStereo(mol)
            if s.type == Chem.StereoType.Atom_Tetrahedral
        }
    atrop = {
        i
        for b in mol.GetBonds()
        if b.GetStereo() in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW)
        for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    }
    unspecified = Chem.ChiralType.CHI_UNSPECIFIED
    out: dict[int, str] = {}
    centres: set[int] = set()
    for bid, bond in wedged.items():
        i = index[bond.begin]
        if i in atrop:
            continue
        tag = mol.GetAtomWithIdx(i).GetChiralTag()
        if i not in potential and tag == unspecified:
            out[bid] = WRONG_END if index[bond.end] in potential else NO_CENTRE
        elif tag == unspecified:
            out[bid] = CONTRADICTORY
        else:
            centres.add(bond.begin)
    for centre in sorted(centres):
        if _ambiguous(sketch, centre, mol.GetAtomWithIdx(index[centre]).GetChiralTag(), index):
            for bid, bond in wedged.items():
                if bond.begin == centre:
                    out[bid] = AMBIGUOUS
    return out


def _ambiguous(sketch: Sketch, centre: int, tag, index: dict[int, int]) -> bool:
    """Does RDKit read *centre* differently when one of its bonds is turned by NUDGE?"""
    cx, cy = sketch.atoms[centre].pos
    for n in sketch.neighbours(centre):
        for turn in (NUDGE, -NUDGE):
            trial = sketch.copy()
            atom = trial.atoms[n]
            t = math.radians(turn)
            dx, dy = atom.x - cx, atom.y - cy
            c, s = math.cos(t), math.sin(t)
            trial.atoms[n] = atom.moved(cx + dx * c - dy * s, cy + dx * s + dy * c)
            try:
                read = to_rdkit(trial).GetAtomWithIdx(index[centre]).GetChiralTag()
            except SketchError:
                return True
            if read != tag:
                return True
    return False


def set_centre(change: Change, aid: int, want: str) -> str:
    """Make atom *aid* R or S (r/s for pseudo-asymmetric centres). Returns "" or why not."""
    sketch = change.sketch
    now = labels(sketch)
    current = now.atoms.get(aid)
    if current is not None and current.upper() == want.upper():
        return ""
    if current is None and aid not in now.undefined_atoms:
        return "This atom is not a stereocentre."
    if current is not None:
        own = [
            bid
            for bid, b in sketch.bonds.items()
            if b.begin == aid and b.stereo in _SWAP and b.order is Order.SINGLE
        ]
        for bid in own:
            bond = sketch.bonds[bid]
            change.set_bond(bid, replace(bond, stereo=_SWAP[bond.stereo]))
        if not own:
            return "Its stereo comes from wedges drawn at a neighbour; edit those bonds instead."
    else:
        bid = _wedge_candidate(sketch, aid, now)
        if bid is None:
            return "There is no single bond at this atom that can carry a wedge."
        bond = sketch.bonds[bid]
        if bond.begin != aid:
            bond = bond.flipped()
        change.set_bond(bid, replace(bond, stereo=Stereo.WEDGE))
        if labels(change.sketch).atoms.get(aid, "").upper() != want.upper():
            change.set_bond(bid, replace(bond, stereo=Stereo.HASH))
    got = labels(change.sketch).atoms.get(aid, "")
    if got.upper() != want.upper():
        return f"RDKit does not read this centre as {want} after the change."
    return ""


def _wedge_candidate(sketch: Sketch, aid: int, now: StereoLabels) -> int | None:
    """The bond at *aid* best suited to a new wedge: single, not starting a wedge elsewhere,
    preferably to a hydrogen or a chain end, out of rings, and not to another stereocentre."""
    centres = set(now.atoms) | set(now.undefined_atoms)
    scored = []
    for bid, bond in sketch.bonds_of(aid):
        if bond.order is not Order.SINGLE or bond.stereo is not Stereo.NONE:
            continue
        other = bond.other(aid)
        in_ring = _in_ring(sketch, bid)
        score = (
            sketch.atoms[other].symbol != "H",
            other in centres,
            in_ring,
            sketch.degree(other) > 1,
            bid,
        )
        scored.append((score, bid))
    return min(scored)[1] if scored else None


def _in_ring(sketch: Sketch, bid: int) -> bool:
    bond = sketch.bonds[bid]
    return bond.end in _side(sketch, bond.begin, bid)


def _side(sketch: Sketch, start: int, skip_bond: int) -> set[int]:
    """Atoms reachable from *start* without crossing bond *skip_bond*."""
    seen, stack = {start}, [start]
    while stack:
        here = stack.pop()
        for bid, bond in sketch.bonds_of(here):
            if bid == skip_bond:
                continue
            other = bond.other(here)
            if other not in seen:
                seen.add(other)
                stack.append(other)
    return seen


# -- E / Z ------------------------------------------------------------------------------------
UNKNOWN = "either"


def set_double_bond(change: Change, bid: int, want: str) -> str:
    """Make double bond *bid* E or Z by reflecting one side, or ``UNKNOWN`` (drawn crossed).

    Returns "" or why not.
    """
    sketch = change.sketch
    bond = sketch.bonds.get(bid)
    if bond is None or bond.order is not Order.DOUBLE:
        return "This is not a double bond."
    if want == UNKNOWN:
        if bond.stereo is not Stereo.EITHER:
            change.set_bond(bid, replace(bond, stereo=Stereo.EITHER))
        return ""
    if bond.stereo is Stereo.EITHER:  # back to the geometry as drawn
        bond = replace(bond, stereo=Stereo.NONE)
        change.set_bond(bid, bond)
    now = labels(change.sketch)
    current = now.bonds.get(bid)
    if current == want:
        return ""
    if current is None:
        if bid in now.undefined_bonds:
            return "Draw the double bond's substituents first; it has no E/Z yet."
        return "This double bond has no E/Z stereo (one end has two identical groups)."
    if _in_ring(sketch, bid):
        return "A double bond in a ring cannot be switched between E and Z here."
    a, b = bond.begin, bond.end
    side_a, side_b = _side(sketch, a, bid), _side(sketch, b, bid)
    pivot, moving = (b, side_b) if len(side_b) <= len(side_a) else (a, side_a)
    p, q = sketch.atoms[a].pos, sketch.atoms[b].pos
    for aid in moving:
        if aid == pivot:
            continue
        atom = change.sketch.atoms[aid]
        change.move_atom(aid, *_reflect(atom.pos, p, q))
    for other_bid, other in list(change.sketch.bonds.items()):
        if other.begin in moving and other.stereo in _SWAP and other_bid != bid:
            change.set_bond(other_bid, replace(other, stereo=_SWAP[other.stereo]))
    if labels(change.sketch).bonds.get(bid) != want:
        return f"RDKit does not read this double bond as {want} after the change."
    return ""


def _reflect(point, a, b):
    """*point* mirrored across the line through *a* and *b*."""
    (px, py), (ax, ay), (bx, by) = point, a, b
    dx, dy = bx - ax, by - ay
    length2 = dx * dx + dy * dy
    t = ((px - ax) * dx + (py - ay) * dy) / length2
    fx, fy = ax + t * dx, ay + t * dy
    return 2 * fx - px, 2 * fy - py


def opposite(label: str) -> str:
    return _OPPOSITE[label]


# -- hydrogens without losing stereo ------------------------------------------------------------
def keep_centres(change: Change, before: StereoLabels) -> list[int]:
    """Put back R/S labels that an edit changed (drawing or removing a hydrogen next to a
    wedge changes how the wedge reads). Returns the centres that could not be restored."""
    after = labels(change.sketch)
    lost = []
    for aid, code in before.atoms.items():
        changed = aid in change.sketch.atoms and after.atoms.get(aid) != code
        if changed and set_centre(change, aid, code):
            lost.append(aid)
    return lost


def show_hydrogens(change: Change, atoms, counts: dict[int, int]) -> int:
    """``edits.show_hydrogens`` that keeps every stereocentre's R/S."""
    before = labels(change.sketch)
    added = edits.show_hydrogens(change, atoms, counts)
    keep_centres(change, before)
    return added


def hide_hydrogens(change: Change, atoms) -> tuple[int, int]:
    """``edits.hide_hydrogens`` that keeps every stereocentre's R/S."""
    before = labels(change.sketch)
    result = edits.hide_hydrogens(change, atoms)
    keep_centres(change, before)
    return result
