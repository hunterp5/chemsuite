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
  from the deposited sequence (no loops are built: the ends stay charged termini);
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
MISSING = "missing"
UNKNOWN = "unknown"
LINK = "link"
DISULFIDE = "disulfide"
METAL = "metal"
POLYMER = (AMINO, NUCLEIC)
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
    labels = [s.residue_label(r) for r in sorted(chosen)]
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
            numbered = s.seqres[r] >= 0 and s.seqres[nxt] >= 0
            missing = int(s.seqres[nxt] - s.seqres[r] - 1) if numbered else 0
            if missing > 0:
                found.append(Finding(GAP, label, f"{missing} residues not modelled: the ends are "
                                     "left as charged termini", warn=True, residue=r))  # fmt: skip
            elif (r, nxt) not in linked:
                found.append(Finding(BREAK, label, "No bond between these residues", warn=True,
                                     residue=r))  # fmt: skip
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


def last_in_chain(s: Structure) -> set[int]:
    """The polymer residues that end their chain (the true C-termini, not those at a gap)."""
    out = set()
    for c in range(s.n_chains):
        span = s.residues_of(c)
        polymer = [r for r in range(span.start, span.stop) if s.res_kind[r] == AMINO]
        if polymer:
            out.add(polymer[-1])
    return out


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
