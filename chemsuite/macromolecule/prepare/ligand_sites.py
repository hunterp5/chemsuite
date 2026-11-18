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

A ligand starts as its component draws it (:mod:`.hydrogens`). With chemsuite pKa (a
:class:`.ligand_pka.LigandPka`), its forms are the model's protonation microstates and the best
one's tautomers (RDKit's enumerator, moving hydrogens only between N, O and S), each with the
model's free energy, which holds what the charges do to each other (:func:`_model_forms`).
Without a model (``pka`` None) every ligand stays as drawn, as does one with more than
:data:`~chemsuite.ionisation.unipka.ranges.MAX_MICROSTATES` microstates.

The drawn form is always one of them, and every form carries how many polar atoms it protonates
unlike the drawing (:attr:`.sites.State.unlike`): the drawing is right most of the time, so
another form must win by a margin. The atoms that move are the polar atoms whose hydrogens or
charge differ between the forms; their hydrogens come from RDKit's placement on the pose, and a
terminal OH, SH or NH turns about its bond (in :data:`STEPS`, as fine as :data:`MAX_STATES`
allows). A charge is shared by the polar atoms of its conjugated group (a carboxylate's two O,
an amidinium's two N).

A ligand bonded to another residue keeps the atoms of that bond as they are; a ligand whose
bond orders were guessed, or that RDKit cannot read, stays as drawn.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from chemsuite.core.jobs import CancelToken, check
from chemsuite.ionisation.unipka import ranges

from ..chem import ligands
from ..structure.components import CcdCache
from ..structure.model import LIGAND, Structure
from . import sites
from .hydrogens import HYDROGEN, Note, Protonated, ligand_names
from .ligand_pka import LigandPka
from .sites import KT, POLAR, RT_LN10, Site, State

WINDOW = 4.0  # pH units: forms this much worse than the best at the pH are left out
MAX_FORMS = 16  # forms kept, the best at the pH first
MAX_TAUTOMERS = 4
TAUTOMER_GAP = 10  # RDKit tautomer score units below the best form's: further ones dropped
MAX_STATES = 400  # per ligand, rotor variants included
STEPS = (30.0, 60.0, 120.0, 180.0)  # degrees, the finest that fits first
#: elements the model is used for; a ligand with any other (a haem's iron) stays as drawn
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


def find(
    p: Protonated,
    ph: float,
    ccd: CcdCache | None = None,
    skip: frozenset[int] = frozenset(),
    cancel: CancelToken | None = None,
    pka: LigandPka | None = None,
) -> tuple[list[Site], list[Note]]:
    """A site for each ligand (other than those in *skip*) with more than one form at *ph*,
    its forms from *pka*'s microstates (chemsuite pKa). Without *pka* none: every ligand stays
    as drawn."""
    if pka is None:
        return [], []
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
            site = _site(s, r, ligand.mol, graph, ph, pka, cancel)
        except (ValueError, RuntimeError, ranges.TooManyMicrostates) as error:
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


def _tautomer_mols(
    mol: Chem.Mol, frozen: set[int]
) -> list[tuple[float, Chem.Mol, tuple[int, ...]]]:
    """The :data:`MAX_TAUTOMERS` best tautomers of *mol* by RDKit's score (within
    :data:`TAUTOMER_GAP` of *mol*'s) that move hydrogens only between N, O and S, leave the
    *frozen* atoms and every charge alone: (how much worse it scores, tautomer, hydrogens)."""
    enumerator = rdMolStandardize.TautomerEnumerator()
    reference = enumerator.ScoreTautomer(mol)
    h0, q0 = _counts(mol)
    polar = [a.GetAtomicNum() in POLAR for a in mol.GetAtoms()]
    scored = []
    for t in enumerator.Enumerate(mol):
        if t.GetNumAtoms() != mol.GetNumAtoms():
            continue
        h, q = _counts(t)
        moved = [k for k in range(len(h)) if h[k] != h0[k]]
        if not moved or any(not polar[k] or k in frozen for k in moved) or q != q0:
            continue
        score = enumerator.ScoreTautomer(t)
        if score >= reference - TAUTOMER_GAP:
            scored.append((score, t, h))
    scored.sort(key=lambda x: -x[0])
    return [(float(reference - score), t, h) for score, t, h in scored[:MAX_TAUTOMERS]]


def _model_forms(
    neutral: Chem.Mol,
    frozen: set[int],
    ph: float,
    drawn: Chem.Mol,
    pka: LigandPka,
    cancel: CancelToken | None = None,
) -> list[Form]:
    """*drawn*, then the microstates of *neutral* (and the best one's tautomers) within
    :data:`WINDOW` pH units of the best at *ph*, best first, each with the model's free
    energy: ``intrinsic = kT · G − RT ln10 · n · μ``, so with the solver's pH term a form
    costs ``G + ln10 · n · (pH − μ)`` kT. The model's G already holds what the charges on
    the molecule do to each other, so ``within``, ``attract`` and ``tautomer_cost`` are 0.

    The charge range widens while the best microstate holds the lowest or highest charge
    found (:mod:`chemsuite.ionisation.unipka.ranges`). Microstates that change a *frozen*
    atom are left out. Raises :class:`~chemsuite.ionisation.unipka.ranges.TooManyMicrostates`.
    """
    h0, q0 = _counts(neutral)
    scored: dict[tuple[tuple[int, ...], tuple[int, ...]], Form] = {}

    def form(mol: Chem.Mol) -> Form:
        h, q = _counts(mol)
        if (h, q) not in scored:
            check(cancel)
            protons = sum(h) - sum(h0)
            g = KT * pka.free_energy(mol) - RT_LN10 * protons * pka.mean
            scored[h, q] = Form(mol, h, q, protons, g)
        return scored[h, q]

    def at_ph(f: Form) -> float:
        return f.intrinsic + RT_LN10 * f.protons * ph

    found = _microstates(neutral, frozen, pka, form, at_ph)
    if found:
        best = min(found, key=at_ph)
        for _gap, t, h in _tautomer_mols(best.mol, frozen):
            mol = _with(t, list(h), list(best.charges))
            if mol is not None:
                found.append(form(mol))
    as_drawn = form(drawn)
    floor = min(at_ph(f) for f in [as_drawn, *found])
    kept = sorted((f for f in found if at_ph(f) - floor <= RT_LN10 * WINDOW), key=at_ph)
    return [as_drawn, *kept[:MAX_FORMS]]


def _microstates(
    neutral: Chem.Mol,
    frozen: set[int],
    pka: LigandPka,
    form: Callable[[Chem.Mol], Form],
    at_ph: Callable[[Form], float],
) -> list[Form]:
    """The microstates of *neutral* that leave the *frozen* atoms alone, as forms; the charge
    range widens (:func:`chemsuite.ionisation.unipka.ranges.widen`) while the best (by
    *at_ph*) holds the lowest or highest charge found."""
    h0, q0 = _counts(neutral)

    def accept(pools: ranges.Pools) -> list[Form]:
        found = []
        for m in (m for pool in pools.values() for m in pool):
            if any(m.hydrogens[k] != h0[k] or m.charges[k] != q0[k] for k in frozen):
                continue
            mol = _with(neutral, list(m.hydrogens), list(m.charges))
            if mol is not None:
                found.append(form(mol))
        return found

    def ends(found: list[Form], low: int, high: int) -> tuple[bool, bool]:
        charge = sum(min(found, key=at_ph).charges)
        return charge == low, charge == high

    return ranges.widen(lambda lo, hi: pka.microstates(neutral, lo, hi), accept, ends)


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
def _site(
    s: Structure,
    r: int,
    mol: Chem.Mol,
    graph: list[list[int]],
    ph: float,
    pka: LigandPka,
    cancel: CancelToken | None = None,
) -> Site | None:
    index = [s.index_of_uid[a.GetIntProp("uid")] for a in mol.GetAtoms()]
    outside = [any(s.residue_of[j] != r and s.element[j] != HYDROGEN for j in graph[i])
               for i in index]  # fmt: skip
    frozen = {k for k, bonded in enumerate(outside) if bonded}
    drawn = _drawn(s, mol, index, graph)
    neutral = rdMolStandardize.Uncharger().uncharge(drawn)
    if neutral.GetNumAtoms() != drawn.GetNumAtoms():
        return None
    Chem.SanitizeMol(neutral)
    forms = _unique(_model_forms(neutral, frozen, ph, drawn, pka, cancel))
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
                         unlike=float(unlike)))  # fmt: skip
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
