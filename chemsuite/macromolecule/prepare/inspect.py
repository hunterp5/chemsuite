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

"""What a structure needs before its hydrogens can be placed, as report rows (:class:`Finding`).

- alternate locations (the one shown is kept);
- chain breaks: neighbouring residues of a chain with no peptide bond, and residues missing
  from the deposited sequence, inside a chain or at its ends (no loops are built: those ends
  take neutral ACE and NME caps, :func:`ends`);
- missing heavy atoms of polymer residues, against their templates;
- residues with no template (not built, hydrogens from RDKit);
- covalent links between a ligand or glycan and the polymer, disulfides, metal sites.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from ..structure.model import AMINO, ION, KIND_NAMES, NUCLEIC, Structure
from .hydrogens import HYDROGEN, METAL_CONTACT, SOFT_METALS
from .templates import Library

ALTLOCS = "altlocs"
BREAK = "break"
GAP = "gap"
TRUNCATED = "truncated"
CAP = "cap"
FILLED = "filled"  # a residue missing alone, built (:func:`.repair.fill_single_gaps`)
REMOVED = "removed"  # too incomplete to build on or cap (:func:`.repair.remove_stubs`)
MISSING = "missing"
UNKNOWN = "unknown"
LINK = "link"
DISULFIDE = "disulfide"
METAL = "metal"
POLYMER = (AMINO, NUCLEIC)
#: residues that end a chain as caps, not as termini
CAPS = frozenset(["ACE", "NME", "NH2"])
#: atoms a residue needs before its missing atoms are built
BACKBONE = {AMINO: ("N", "CA", "C"), NUCLEIC: ("C1'", "C4'")}


@dataclass(frozen=True)
class Finding:
    kind: str
    label: str  # the residue (or residues) it concerns, as the viewer names them
    text: str
    warn: bool = False  # needs the user's eye (⚠)
    residue: int = -1  # in the structure it was found on (-1 for none)
    atoms: tuple[str, ...] = ()  # the atoms concerned, by name (MISSING)


def inspect(s: Structure, library: Library, cancel: CancelToken | None = None) -> list[Finding]:
    found: list[Finding] = []
    found += _altlocs(s)
    found += _breaks(s)
    found += _missing(s, library, cancel)
    found += _links(s)
    found += _metals(s)
    return found


def _altlocs(s: Structure) -> list[Finding]:
    chosen = s.default_altlocs
    if not chosen:
        return []
    labels = [s.residue_label(r) for r in sorted(chosen) if s.altloc_group[r] == r]
    shown = ", ".join(labels[:6]) + (" …" if len(labels) > 6 else "")
    text = f"{len(labels)} residues with alternate locations: the location shown kept"
    return [Finding(ALTLOCS, shown, text)]


def _breaks(s: Structure) -> list[Finding]:
    found = []
    linked = _peptide_links(s)
    for c in range(s.n_chains):
        span = s.residues_of(c)
        for r in range(span.start, span.stop - 1):
            nxt = r + 1
            if s.res_kind[r] not in POLYMER or s.res_kind[nxt] != s.res_kind[r]:
                continue
            label = f"{s.residue_label(r)} – {s.residue_label(nxt)}"
            missing = _missing_between(s, r, nxt)
            if missing > 0:
                found.append(Finding(GAP, label, f"{missing} residues not modelled", warn=True,
                                     residue=r))  # fmt: skip
            elif (r, nxt) not in linked:
                found.append(Finding(BREAK, label, "No bond between these residues", warn=True,
                                     residue=r))  # fmt: skip
    found += _truncated(s)
    return found


def _truncated(s: Structure) -> list[Finding]:
    """The chain ends that are not the ends of the deposited sequence."""
    found = []
    for c in range(s.n_chains):
        span = s.residues_of(c)
        amino = [r for r in range(span.start, span.stop) if s.res_kind[r] == AMINO]
        if not amino:
            continue
        head, tail = amino[0], amino[-1]
        if s.seqres[head] > 0 and str(s.res_name[head]) not in CAPS:
            found.append(Finding(TRUNCATED, s.residue_label(head),
                                 f"{int(s.seqres[head])} N-terminal residues not modelled",
                                 residue=head))  # fmt: skip
        if _more_after(s, tail) and str(s.res_name[tail]) not in CAPS:
            more = len(s.sequence_of(int(s.entity[tail]))) - 1 - int(s.seqres[tail])
            found.append(
                Finding(
                    TRUNCATED,
                    s.residue_label(tail),
                    f"{more} C-terminal residues not modelled",
                    residue=tail,
                )
            )
    return found


def _peptide_links(s: Structure) -> set[tuple[int, int]]:
    heads = {"C", "O3'"}
    tails = {"N", "P"}
    found = set()
    for a, b in s.bonds.tolist():
        ra, rb = int(s.residue_of[a]), int(s.residue_of[b])
        if rb == ra + 1 and s.name[a] in heads and s.name[b] in tails:
            found.add((ra, rb))
        elif ra == rb + 1 and s.name[b] in heads and s.name[a] in tails:
            found.add((rb, ra))
    return found


@dataclass(frozen=True)
class Ends:
    """Where a protein's chains end: the true termini, and the ends of what was modelled."""

    n_termini: frozenset[int]  # residues that start their sequence
    c_termini: frozenset[int]  # residues that end it (they take an OXT)
    open_n: frozenset[int]  # residues after a gap or a missing N-terminal stretch (ACE)
    open_c: frozenset[int]  # residues before a gap or a missing C-terminal stretch (NME)


def ends(s: Structure) -> Ends:
    """The ends of every chain's stretches of amino acids. A stretch ends at a gap in the
    deposited sequence (``seqres``) or at the chain's end; a chain's end is a true terminus
    unless the sequence goes on beyond it, or it is a cap (ACE, NME, NH2) already. Without a
    deposited sequence every chain end is a true terminus."""
    n_term, c_term, open_n, open_c = set(), set(), set(), set()
    for c in range(s.n_chains):
        span = s.residues_of(c)
        amino = [r for r in range(span.start, span.stop) if s.res_kind[r] == AMINO]
        if not amino:
            continue
        stretches = [[amino[0]]]
        for r, nxt in zip(amino, amino[1:], strict=False):
            if _missing_between(s, r, nxt) > 0:
                stretches.append([])
            stretches[-1].append(nxt)
        for k, stretch in enumerate(stretches):
            head, tail = stretch[0], stretch[-1]
            if str(s.res_name[head]) not in CAPS:
                truncated = k > 0 or s.seqres[head] > 0
                (open_n if truncated else n_term).add(head)
            if str(s.res_name[tail]) not in CAPS:
                truncated = k < len(stretches) - 1 or _more_after(s, tail)
                (open_c if truncated else c_term).add(tail)
    return Ends(frozenset(n_term), frozenset(c_term), frozenset(open_n), frozenset(open_c))


def _missing_between(s: Structure, r: int, nxt: int) -> int:
    """How many residues of the deposited sequence lie between *r* and *nxt*."""
    if s.seqres[r] < 0 or s.seqres[nxt] < 0 or s.entity[r] != s.entity[nxt]:
        return 0
    return int(s.seqres[nxt] - s.seqres[r] - 1)


def _more_after(s: Structure, r: int) -> bool:
    """Whether the deposited sequence goes on after residue *r*."""
    e = int(s.entity[r])
    if s.seqres[r] < 0 or not 0 <= e < len(s.entity_id):
        return False
    return int(s.seqres[r]) < len(s.sequence_of(e)) - 1


def last_in_chain(s: Structure) -> set[int]:
    """The residues that end their chain's sequence, which take an OXT (not those before a gap
    or a missing C-terminal stretch, which take an NME cap)."""
    return set(ends(s).c_termini)


def missing_atoms(
    s: Structure, r: int, library: Library, ends: set[int], cancel: CancelToken | None = None
) -> list[str] | None:
    """The template's heavy atoms residue *r* lacks (None when it has no template); *ends*
    are the residues that end their chains (:func:`last_in_chain`), which need an OXT."""
    template = library.get(str(s.res_name[r]), cancel)
    if template is None:
        return None
    active = s.active()
    span = s.atoms_of(r)
    here = {str(s.name[i]) for i in range(span.start, span.stop) if active[i]}
    leaving = set(template.leaving)
    if s.res_kind[r] == AMINO and r in ends:
        leaving.discard("OXT")
    return [a for a in template.heavy if a not in here and a not in leaving]


def _missing(s: Structure, library: Library, cancel: CancelToken | None) -> list[Finding]:
    found = []
    ends = last_in_chain(s)
    for r in range(s.n_residues):
        if s.res_kind[r] not in POLYMER:
            continue
        check(cancel)
        lacking = missing_atoms(s, r, library, ends, cancel)
        label = s.residue_label(r)
        if lacking is None:
            found.append(Finding(UNKNOWN, label, "No template: its hydrogens come from RDKit",
                                 warn=True, residue=r))  # fmt: skip
        elif lacking == ["OXT"] and r in ends:
            found.append(Finding(MISSING, label, "C-terminal OXT added", residue=r,
                                 atoms=("OXT",)))  # fmt: skip
        elif lacking:
            backbone = set(BACKBONE[int(s.res_kind[r])]) & set(lacking)
            text = (f"Missing backbone atoms {', '.join(sorted(backbone))}: not rebuilt"
                    if backbone else f"Missing atoms built: {', '.join(lacking)}")  # fmt: skip
            found.append(Finding(MISSING, label, text, warn=True, residue=r, atoms=tuple(lacking)))
    return found


def _links(s: Structure) -> list[Finding]:
    found = []
    pairs = set()
    for a, b in s.bonds.tolist():
        ra, rb = int(s.residue_of[a]), int(s.residue_of[b])
        if ra == rb or s.element[a] == HYDROGEN or s.element[b] == HYDROGEN:
            continue
        ka, kb = int(s.res_kind[ra]), int(s.res_kind[rb])
        names = (str(s.name[a]), str(s.name[b]))
        if names == ("SG", "SG"):
            pairs.add((DISULFIDE, min(ra, rb), max(ra, rb), "Disulfide"))
        elif (ka in POLYMER) != (kb in POLYMER) and ION not in (ka, kb):
            r, other = (ra, rb) if ka in POLYMER else (rb, ra)
            kind = KIND_NAMES.get(int(s.res_kind[other]), "residue").lower()
            pairs.add((LINK, r, other, f"Covalent link to a {kind}"))
    for kind, ra, rb, text in sorted(pairs):
        found.append(Finding(kind, f"{s.residue_label(ra)} – {s.residue_label(rb)}", text,
                             residue=ra))  # fmt: skip
    return found


def _metals(s: Structure) -> list[Finding]:
    metals = np.flatnonzero((s.atom_kind == ION) & np.isin(s.element, list(SOFT_METALS)))
    if not len(metals):
        return []
    x = s.xyz()
    polar = np.flatnonzero(np.isin(s.element, [7, 8, 16]) & (s.atom_kind != ION))
    found = []
    for m in metals.tolist():
        d = np.linalg.norm(x[polar] - x[m], axis=1)
        near = polar[d < METAL_CONTACT]
        partners = sorted({s.residue_label(int(s.residue_of[i])) for i in near})
        r = int(s.residue_of[m])
        found.append(Finding(METAL, s.residue_label(r),
                             f"Bound to {', '.join(partners)}" if partners else "No partners",
                             residue=r))  # fmt: skip
    return found
