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

"""Ligand sites: each ligand's protonation states and tautomers, to be scored in its pocket.

A ligand starts as its component draws it (:mod:`.hydrogens`). Its forms are:

- **Protomers.** MolGpKa (:mod:`chemsuite.ionisation.molgpka`) gives the pKa values of the
  neutral parent's sites. Each ionisable atom takes the charge levels (−1, 0, +1) whose free
  energy at the pH is within :data:`WINDOW` pH units of its best. The :data:`MAX_SITES` atoms
  closest to a tie take all of them; the rest stay in their best level. A form's free energy in
  water is −RT ln10 · Σ pKa over the protons it carries, as for the protein's groups
  (:attr:`.sites.State.intrinsic`), plus the Coulomb energy between the charges it adds
  (:attr:`.sites.State.within`): MolGpKa gives each site's pKa with the rest neutral, and a
  second charge on one molecule costs. A base on an amine N conjugated with another base
  (adenine's N6, beside N1) is left out.
- **Tautomers** of the neutral parent (RDKit's enumerator) that move hydrogens only between N,
  O and S, the :data:`MAX_TAUTOMERS` best by RDKit's score, each carrying how much worse it
  scores than the drawn form (:attr:`.sites.State.tautomer_cost`).

The drawn form is always one of them, and every form carries how many polar atoms it protonates
unlike the drawing (:attr:`.sites.State.unlike`): the drawing is right most of the time, so
another form must win by a margin. The atoms that move are the polar atoms whose hydrogens or
charge differ between the forms; their hydrogens come from RDKit's placement on the pose, and a
terminal OH, SH or NH turns about its bond (in :data:`STEPS`, as fine as :data:`MAX_STATES`
allows). A charge is shared by the polar atoms of its conjugated group (a carboxylate's two O,
an amidinium's two N) for the Coulomb term.

A ligand bonded to another residue keeps the atoms of that bond as they are; a ligand whose
bond orders were guessed, or that RDKit cannot read, stays as drawn.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, replace

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from chemsuite.core.jobs import CancelToken, check
from chemsuite.ionisation import molgpka

from ..chem import ligands
from ..structure.components import CcdCache
from ..structure.model import LIGAND, Structure
from . import sites
from .hydrogens import HYDROGEN, Note, Protonated, ligand_names
from .sites import KT, POLAR, RT_LN10, Site, State

WINDOW = 4.0  # pH units: charge levels this much worse than an atom's best are left out
MAX_SITES = 4  # ionisable atoms that take more than one level
MAX_FORMS = 16  # protomers kept, the best at the pH first
MAX_TAUTOMERS = 4
TAUTOMER_GAP = 10  # RDKit tautomer score units below the drawn form's: further ones dropped
MAX_STATES = 400  # per ligand, rotor variants included
STEPS = (30.0, 60.0, 120.0, 180.0)  # degrees, the finest that fits first
#: charges on one ligand: 332/(εr) with ε = 4r as in :func:`.energy.coulomb`, scaled to the
#: shifts seen in small molecules (a phosphate's second pKa, glycine's amine), the distance
#: held to at least :data:`WITHIN_FLOOR` Å
WITHIN_SCALE, WITHIN_FLOOR = 0.5, 3.0
#: elements MolGpKa knows; a ligand with any other (a haem's iron) stays as drawn
ORGANIC = frozenset([1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53])
AS_DRAWN = "as drawn"


@dataclass(frozen=True)
class Form:
    """One protomer or tautomer: hydrogens and formal charge on each atom of the ligand mol."""

    mol: Chem.Mol  # heavy atoms, explicit hydrogen counts, sanitised
    hydrogens: tuple[int, ...]
    charges: tuple[int, ...]
    protons: int
    intrinsic: float
    tautomer_cost: float = 0.0
    within: float = 0.0  # kcal/mol between the like charges the form adds to the neutral parent
    attract: float = 0.0  # and between its opposite ones (negative)


def find(
    p: Protonated,
    ph: float,
    ccd: CcdCache | None = None,
    skip: frozenset[int] = frozenset(),
    cancel: CancelToken | None = None,
) -> tuple[list[Site], list[Note]]:
    """A site for each ligand (other than those in *skip*) with more than one form at *ph*."""
    s = p.structure
    graph = sites.bond_graph(s)
    found, notes = [], []
    for r in np.flatnonzero(s.res_kind == LIGAND).tolist():
        check(cancel)
        if r in skip:
            continue
        try:
            ligand = ligands.to_mol(s, r, ccd=ccd, cancel=cancel)
        except (ValueError, RuntimeError):
            continue  # hydrogens.py has said why
        if ligand.guessed:
            notes.append(Note(r, "Kept as drawn: its bond orders were guessed"))
            continue
        if any(a.GetAtomicNum() not in ORGANIC for a in ligand.mol.GetAtoms()):
            continue  # a metal complex: no pKa model for it
        if ligand.mol.GetNumAtoms() < 2:
            continue  # an ion (OH, D3O, ...): the component names its form
        try:
            site = _site(s, r, ligand.mol, graph, ph)
        except (ValueError, RuntimeError) as error:
            notes.append(Note(r, f"Kept as drawn: {error}"))
            continue
        if site is not None:
            found.append(site)
    return found, notes


# -- forms -----------------------------------------------------------------------------------------
def _drawn(s: Structure, mol: Chem.Mol, index: list[int], graph: list[list[int]]) -> Chem.Mol:
    """*mol* with the hydrogens and charges the structure gives each atom."""
    out = Chem.RWMol(mol)
    for atom, i in zip(out.GetAtoms(), index, strict=True):
        atom.SetNumExplicitHs(sum(1 for j in graph[i] if s.element[j] == HYDROGEN))
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(int(s.charge[i]))
    drawn = out.GetMol()
    Chem.SanitizeMol(drawn)
    return drawn


def _counts(mol: Chem.Mol) -> tuple[tuple[int, ...], tuple[int, ...]]:
    atoms = list(mol.GetAtoms())
    return tuple(a.GetTotalNumHs() for a in atoms), tuple(a.GetFormalCharge() for a in atoms)


def _with(mol: Chem.Mol, hydrogens: list[int], charges: list[int]) -> Chem.Mol | None:
    out = Chem.RWMol(mol)
    for atom, h, q in zip(out.GetAtoms(), hydrogens, charges, strict=True):
        atom.SetNumExplicitHs(h)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(q)
    made = out.GetMol()
    try:
        Chem.SanitizeMol(made)
    except (ValueError, RuntimeError):
        return None
    return made


def _levels(found: list[molgpka.Site]) -> dict[int, list[tuple[int, float]]]:
    """Each ionisable atom's charge levels (−1, 0, +1), lowest first, with their free energies
    in water less the pH term (−RT ln10 · pKa per proton above the lowest)."""
    by_atom: dict[int, dict[str, float]] = {}
    for site in found:
        by_atom.setdefault(site.atom, {})[site.kind] = site.pka
    out = {}
    for atom, pkas in by_atom.items():
        acid, base = pkas.get(molgpka.ACID), pkas.get(molgpka.BASE)
        levels = [(-1, 0.0)] if acid is not None else []
        g0 = -RT_LN10 * acid if acid is not None else 0.0
        levels.append((0, g0))
        if base is not None:
            levels.append((1, g0 - RT_LN10 * base))
        out[atom] = levels
    return out


def _plausible(neutral: Chem.Mol, found: list[molgpka.Site]) -> list[molgpka.Site]:
    """*found* less the sites MolGpKa gets wrong by chemistry:

    - bases whose lone pair is not free: an aromatic N with three neighbours (a pyrrole N,
      adenine's N9), and an amine N whose lone pair is shared with a conjugated system that
      has another base (adenine's N6, beside N1);
    - acids on an alcohol, an O whose one neighbour is an sp³ carbon (pKa ≈ 15; MolGpKa puts
      NADP's ribose hydroxyls at 2–7)."""
    supplier = Chem.ResonanceMolSupplier(neutral)
    group = [supplier.GetAtomConjGrpIdx(a.GetIdx()) for a in neutral.GetAtoms()]
    bases = [site.atom for site in found if site.kind == molgpka.BASE]

    def amine(k: int) -> bool:
        atom = neutral.GetAtomWithIdx(k)
        single = all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds())
        return (atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() > 0 and not atom.IsInRing()
                and group[k] >= 0 and single)  # fmt: skip

    def shared(k: int) -> bool:
        return any(b != k and group[b] == group[k] for b in bases)

    def pyrrole(k: int) -> bool:
        atom = neutral.GetAtomWithIdx(k)
        return atom.GetIsAromatic() and atom.GetDegree() + atom.GetTotalNumHs() >= 3

    def alcohol(k: int) -> bool:
        atom = neutral.GetAtomWithIdx(k)
        heavy = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        return (atom.GetAtomicNum() == 8 and len(heavy) == 1 and heavy[0].GetAtomicNum() == 6
                and heavy[0].GetHybridization() == Chem.HybridizationType.SP3)  # fmt: skip

    def wrong(site: molgpka.Site) -> bool:
        if site.kind == molgpka.ACID:
            return alcohol(site.atom)
        return pyrrole(site.atom) or (amine(site.atom) and shared(site.atom))

    return [site for site in found if not wrong(site)]


def _couplings(neutral: Chem.Mol, atoms: list[int], xyz: np.ndarray) -> np.ndarray:
    """The Coulomb energy (ε = 4r, as :func:`.energy.coulomb`) between like unit charges on
    each pair of *atoms*, each charge shared by the polar atoms of its conjugated group. Only
    like charges count (``level_a · level_b > 0``): MolGpKa's pKa values already carry most of
    the pull between opposite ones, which zwitterions in its training set show."""
    supplier = Chem.ResonanceMolSupplier(neutral)
    spread = []
    for k in atoms:
        group = supplier.GetAtomConjGrpIdx(k)
        members = [a.GetIdx() for a in neutral.GetAtoms()
                   if group >= 0 and supplier.GetAtomConjGrpIdx(a.GetIdx()) == group
                   and a.GetAtomicNum() in POLAR] or [k]  # fmt: skip
        spread.append(members)
    out = np.zeros((len(atoms), len(atoms)))
    for a, b in itertools.combinations(range(len(atoms)), 2):
        d = np.linalg.norm(xyz[spread[a]][:, None] - xyz[spread[b]][None], axis=-1)
        e = float((83.0 / np.maximum(d, WITHIN_FLOOR) ** 2).sum()) * WITHIN_SCALE
        e /= len(spread[a]) * len(spread[b])
        out[a, b] = out[b, a] = e
    return out


def _repulsion(row: np.ndarray, level: int, current: list[int], k: int, like: bool = True) -> float:
    """What atom *k* at *level* pays for the like (or, with *like* False, gains from the
    opposite) charges of the others at *current*."""
    return float(sum(row[j] * level * q for j, q in enumerate(current)
                     if j != k and (level * q > 0 if like else level * q < 0)))  # fmt: skip


def _protomers(
    neutral: Chem.Mol, frozen: set[int], ph: float, xyz: np.ndarray, drawn: Chem.Mol | None = None
) -> tuple[list[Form], dict[int, int], Form | None]:
    """The protomers of *neutral* at *ph* (*xyz*: its pose), the level each atom outside the
    choice keeps, and *drawn* as a form scored the same way (every form's free energy sums
    every ionisable atom's level, so the drawing competes on equal terms).

    Each atom's level is first settled in the field of the others' (iterated to a fixed
    point), so a phosphate's second oxygen stays protonated beside its first; the
    :data:`MAX_SITES` atoms closest to a tie in that field, within :data:`WINDOW`, stay open.
    """
    levels = _levels(_plausible(neutral, molgpka.predict_on(neutral)))
    atoms = [a for a in levels if a not in frozen]
    coupling = _couplings(neutral, atoms, xyz)

    def energies(k: int, current: list[int]) -> np.ndarray:
        lv = levels[atoms[k]]
        return np.array([g + RT_LN10 * (level - lv[0][0]) * ph
                         + _repulsion(coupling[k], level, current, k)
                         for level, g in lv])  # fmt: skip

    current = [0] * len(atoms)
    for _ in range(20):
        moved = False
        for k in range(len(atoms)):
            best = levels[atoms[k]][int(np.argmin(energies(k, current)))][0]
            if best != current[k]:
                current[k], moved = best, True
        if not moved:
            break
    allowed, gaps = {}, {}
    for k in range(len(atoms)):
        e = energies(k, current)
        keep = [lv for lv, x in zip(levels[atoms[k]], e, strict=True)
                if x - e.min() <= RT_LN10 * WINDOW]  # fmt: skip
        if len(keep) > 1:
            allowed[k] = keep
            gaps[k] = float(np.sort(e)[1] - e.min())
    open_k = sorted(allowed, key=lambda k: gaps[k])[:MAX_SITES]
    fixed = {atoms[k]: current[k] for k in range(len(atoms)) if k not in open_k}
    h0, q0 = _counts(neutral)

    def form(mol: Chem.Mol, level: list[int]) -> Form:
        protons = sum(lv - levels[a][0][0] for a, lv in zip(atoms, level, strict=True))
        intrinsic = sum(_free_energy(levels[a], lv) for a, lv in zip(atoms, level, strict=True))
        within = sum(_repulsion(coupling[k], lv, level, k) for k, lv in enumerate(level)) / 2
        attract = sum(_repulsion(coupling[k], lv, level, k, like=False)
                      for k, lv in enumerate(level)) / 2  # fmt: skip
        h, q = _counts(mol)
        return Form(mol, h, q, protons, intrinsic, within=within, attract=attract)

    forms = []
    for chosen in itertools.product(*(allowed[k] for k in open_k)):
        level = list(current)
        for k, (lv, _g) in zip(open_k, chosen, strict=True):
            level[k] = lv
        h, q = list(h0), list(q0)
        for k, lv in enumerate(level):
            h[atoms[k]] += lv
            q[atoms[k]] += lv
        mol = _with(neutral, h, q) if min(h) >= 0 else None
        if mol is not None:
            forms.append(form(mol, level))
    forms.sort(key=lambda f: f.intrinsic + f.within + RT_LN10 * f.protons * ph)
    as_drawn = None
    if drawn is not None:
        level = [int(np.clip(drawn.GetAtomWithIdx(a).GetFormalCharge()
                             - neutral.GetAtomWithIdx(a).GetFormalCharge(), -1, 1))
                 for a in atoms]  # fmt: skip
        as_drawn = form(drawn, level)
    return forms[:MAX_FORMS], fixed, as_drawn


def _free_energy(levels: list[tuple[int, float]], level: int) -> float:
    """The free energy of *level* among an atom's *levels* (the nearest one it has)."""
    return min(levels, key=lambda lv: abs(lv[0] - level))[1]


def _tautomers(
    neutral: Chem.Mol, frozen: set[int], fixed: dict[int, int], base: Form | None
) -> list[Form]:
    """The neutral parent's best tautomers that move hydrogens only between N, O and S, with
    the levels *fixed* applied and *base*'s free energy (the protomer with every open atom
    neutral)."""
    if base is None:
        return []
    enumerator = rdMolStandardize.TautomerEnumerator()
    reference = enumerator.ScoreTautomer(neutral)
    h0, q0 = _counts(neutral)
    polar = [a.GetAtomicNum() in POLAR for a in neutral.GetAtoms()]
    scored = []
    for t in enumerator.Enumerate(neutral):
        if t.GetNumAtoms() != neutral.GetNumAtoms():
            continue
        h, q = _counts(t)
        moved = [k for k in range(len(h)) if h[k] != h0[k]]
        if not moved or any(not polar[k] or k in frozen for k in moved) or q != q0:
            continue
        score = enumerator.ScoreTautomer(t)
        if score >= reference - TAUTOMER_GAP:
            scored.append((score, t, h))
    scored.sort(key=lambda x: -x[0])
    forms = []
    for score, t, h in scored[:MAX_TAUTOMERS]:
        h, q = list(h), list(q0)
        for atom, level in fixed.items():
            h[atom] += level
            q[atom] += level
        mol = _with(t, h, q) if min(h) >= 0 else None
        if mol is not None:
            forms.append(Form(mol, tuple(h), tuple(q), base.protons, base.intrinsic,
                              float(reference - score), base.within, base.attract))  # fmt: skip
    return forms


def _shared(mol: Chem.Mol) -> list[float]:
    """Each atom's charge, a formal charge shared by the polar atoms of its conjugated group."""
    out = [0.0] * mol.GetNumAtoms()
    supplier = Chem.ResonanceMolSupplier(mol)
    for atom in mol.GetAtoms():
        q = atom.GetFormalCharge()
        if not q:
            continue
        group = supplier.GetAtomConjGrpIdx(atom.GetIdx())
        members = [a.GetIdx() for a in mol.GetAtoms()
                   if group >= 0 and supplier.GetAtomConjGrpIdx(a.GetIdx()) == group
                   and a.GetAtomicNum() in POLAR]  # fmt: skip
        if atom.GetIdx() not in members:
            members = [atom.GetIdx()]
        for k in members:
            out[k] += q / len(members)
    return out


# -- the site --------------------------------------------------------------------------------------
def _site(s: Structure, r: int, mol: Chem.Mol, graph: list[list[int]], ph: float) -> Site | None:
    index = [s.index_of_uid[a.GetIntProp("uid")] for a in mol.GetAtoms()]
    outside = [any(s.residue_of[j] != r and s.element[j] != HYDROGEN for j in graph[i])
               for i in index]  # fmt: skip
    frozen = {k for k, bonded in enumerate(outside) if bonded}
    drawn = _drawn(s, mol, index, graph)
    neutral = rdMolStandardize.Uncharger().uncharge(drawn)
    if neutral.GetNumAtoms() != drawn.GetNumAtoms():
        return None
    Chem.SanitizeMol(neutral)
    pose = s.coords[0, index].astype(np.float64)
    protomers, fixed, as_drawn = _protomers(neutral, frozen, ph, pose, drawn)
    base = next((f for f in protomers if all(f.charges[a] == neutral.GetAtomWithIdx(a)
                                             .GetFormalCharge() + fixed.get(a, 0)
                                             for a in range(len(index)))), None)  # fmt: skip
    drawn_form = [as_drawn] if as_drawn is not None else []
    forms = _unique(drawn_form + protomers + _tautomers(neutral, frozen, fixed, base))
    if not forms:
        return None
    h_drawn, q_drawn = _counts(drawn)
    shared = [_shared(f.mol) for f in forms]
    polar = [a.GetAtomicNum() in POLAR for a in mol.GetAtoms()]
    moving = [k for k in range(len(index)) if polar[k] and k not in frozen and (
        len({(f.hydrogens[k], f.charges[k]) for f in forms} | {(h_drawn[k], q_drawn[k])}) > 1
        or len({round(q[k], 3) for q in shared}) > 1)]  # fmt: skip
    if not moving or (len(forms) == 1 and all(
            forms[0].hydrogens[k] == h_drawn[k] and forms[0].charges[k] == q_drawn[k]
            for k in moving)):  # fmt: skip
        return None
    unlike = [sum(f.hydrogens[k] != h_drawn[k] for k in moving) for f in forms]
    names = _hydrogen_names(s, r, index, moving, forms, graph)
    budget = max(1, MAX_STATES // len(forms))
    states: list[State] = []
    for f, q, extra in zip(forms, shared, unlike, strict=True):
        label = _label(s, index, f, h_drawn)
        states += _states(s, index, f, q, moving, names, label, extra, budget)
    if len(states) < 2:
        return None
    lowest = min(st.protons for st in states)
    states = [replace(st, protons=st.protons - lowest) for st in states]
    heavy = [index[k] for k in moving]
    atoms = heavy + [j for i in heavy for j in graph[i] if s.element[j] == HYDROGEN]
    inside = set(heavy)
    anchors = sorted({j for i in heavy for j in graph[i]
                      if j not in inside and s.element[j] != HYDROGEN}) or heavy  # fmt: skip
    site = sites.make_site(s, r, sites.LIGAND, np.array(atoms), anchors, states, graph)
    return replace(site, titrating=True)


def _unique(forms: list[Form]) -> list[Form]:
    seen, out = set(), []
    for f in forms:
        key = (f.hydrogens, f.charges)
        if key not in seen:
            seen.add(key)
            out.append(f)
    return out


def _label(s: Structure, index: list[int], f: Form, h_drawn: tuple[int, ...]) -> str:
    """What the form changes from the drawing: "+H N1, −H OXT", or :data:`AS_DRAWN`."""
    parts = []
    for k, i in enumerate(index):
        d = f.hydrogens[k] - h_drawn[k]
        if d:
            sign = "+" if d > 0 else "−"
            parts.append(f"{sign}{abs(d) if abs(d) > 1 else ''}H {s.name[i]}")
    return ", ".join(parts) or AS_DRAWN


def _hydrogen_names(
    s: Structure,
    r: int,
    index: list[int],
    moving: list[int],
    forms: list[Form],
    graph: list[list[int]],
) -> dict[int, list[str]]:
    """Names for the hydrogens on each moving atom: those it has, then new ones."""
    span = s.atoms_of(r)
    taken = {str(s.name[i]) for i in range(span.start, span.stop)}
    out = {}
    for k in moving:
        own = [str(s.name[j]) for j in graph[index[k]] if s.element[j] == HYDROGEN]
        most = max(f.hydrogens[k] for f in forms)
        out[k] = own + ligand_names(None, str(s.name[index[k]]), max(0, most - len(own)), taken)
    return out


def _states(
    s: Structure,
    index: list[int],
    f: Form,
    shared: list[float],
    moving: list[int],
    names: dict[int, list[str]],
    label: str,
    unlike: int,
    budget: int,
) -> list[State]:
    """The form's states: its hydrogens placed by RDKit in every model, then each terminal
    group turned."""
    x = s.coords.astype(np.float64)
    placed = []
    for model in range(s.n_models):
        mol = Chem.Mol(f.mol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for k, i in enumerate(index):
            conf.SetAtomPosition(k, x[model, i].tolist())
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        full = Chem.AddHs(mol, addCoords=True)
        placed.append(full)
    first = placed[0]
    positions = np.stack([m.GetConformer().GetPositions() for m in placed])  # (k, n, 3)
    hs = {k: [n.GetIdx() for n in first.GetAtomWithIdx(k).GetNeighbors()
              if n.GetAtomicNum() == HYDROGEN] for k in moving}  # fmt: skip
    turns = _turns(f.mol, moving, hs, budget)
    out = []
    for angles in itertools.product(*(a for _k, a in turns)):
        xyz = positions.copy()
        for (k, _a), angle in zip(turns, angles, strict=True):
            if angle:
                axis_atom = next(n.GetIdx() for n in f.mol.GetAtomWithIdx(k).GetNeighbors())
                xyz[:, hs[k]] = _turned(xyz[:, hs[k]], xyz[:, axis_atom], xyz[:, k], angle)
        atom_names, parents, elements, coords, formal, partial = [], [], [], [], [], []
        for k in moving:
            atom_names.append(str(s.name[index[k]]))
            parents.append("")
            elements.append(int(f.mol.GetAtomWithIdx(k).GetAtomicNum()))
            coords.append(x[:, index[k]])
            formal.append(f.charges[k])
            partial.append(shared[k])
        for k in moving:
            for h, name in zip(hs[k], names[k], strict=False):
                atom_names.append(name)
                parents.append(str(s.name[index[k]]))
                elements.append(HYDROGEN)
                coords.append(xyz[:, h])
                formal.append(0)
                partial.append(0.0)
        stacked = np.stack(coords, axis=1)
        donor = np.full((len(atom_names), 3), np.nan)
        for m, parent in enumerate(parents):
            if parent:
                donor[m] = stacked[0, atom_names.index(parent)]
        suffix = " ".join(f"{a:.0f}°" for a in angles if a)
        out.append(State(f"{label} {suffix}".strip(), tuple(atom_names), tuple(parents),
                         np.array(elements), stacked, donor, form=label, protons=f.protons,
                         formal=tuple(formal), partial=tuple(partial),
                         unlike=float(unlike), tautomer_cost=f.tautomer_cost,
                         within=f.within, attract=f.attract))  # fmt: skip
    variants = KT * np.log(len(out))
    return [replace(st, intrinsic=f.intrinsic + variants) for st in out]


def _turns(
    mol: Chem.Mol, moving: list[int], hs: dict[int, list[int]], budget: int
) -> list[tuple[int, list[float]]]:
    """The angles (degrees) each terminal group of *mol* turns through, as fine as the budget
    of states allows: a planar OH flips, an sp³ group turns in steps."""
    groups = []
    for k in moving:
        atom = mol.GetAtomWithIdx(k)
        heavy = list(atom.GetNeighbors())
        if not hs[k] or len(heavy) != 1:
            continue
        bond = mol.GetBondBetweenAtoms(k, heavy[0].GetIdx())
        planar = atom.GetHybridization() == Chem.HybridizationType.SP2 or bond.GetIsConjugated()
        if planar:
            if len(hs[k]) == 1:
                groups.append((k, None))
            continue
        groups.append((k, 120.0 if len(hs[k]) == 3 else 360.0))
    for step in STEPS:
        angles = [(k, [0.0, 180.0] if period is None else list(np.arange(0.0, period - 1e-6, step)))
                  for k, period in groups]  # fmt: skip
        if int(np.prod([len(a) for _k, a in angles])) <= budget:
            return angles
    return []


def _turned(h: np.ndarray, axis: np.ndarray, centre: np.ndarray, degrees: float) -> np.ndarray:
    """Hydrogens *h* (k, n, 3) turned by *degrees* about the bond *axis* → *centre* (k, 3)."""
    u = centre - axis
    u = u / np.linalg.norm(u, axis=-1, keepdims=True)
    theta = np.radians(degrees)
    v = h - centre[:, None]
    u = u[:, None]
    cross = np.cross(u, v)
    dot = (u * v).sum(-1, keepdims=True)
    return (
        centre[:, None] + v * np.cos(theta) + cross * np.sin(theta) + u * dot * (1 - np.cos(theta))
    )
