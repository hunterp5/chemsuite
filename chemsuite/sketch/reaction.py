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

"""A drawing with a reaction arrow, as a reaction: components, RDKit, RXN and reaction SMILES.

**Which fragment is what** (docs/REACTIONS_PLAN.md §2) follows Ketcher's ``getComponents``
(``ketcher-core/src/domain/entities/struct.ts`` at 2b031cc2):

1. Fragments are grouped between the plus signs and the arrow. So fragments with no plus
   between them form one component, such as a salt.
2. A component is a reactant when its centre is on the arrow's starting side, otherwise a
   product.

Ketcher measures along x and assumes a left-to-right arrow. Here positions are measured along
the arrow, from tail to head, which is the same thing for that arrow and also handles arrows
drawn the other way or vertically. A **retrosynthetic** arrow points from the target to its
precursors, so its head side holds the reactants.

**Extension: agents.** Ketcher saves anything drawn over the arrow as a product. Here a
fragment whose centre lies along the arrow's span and off its line (above or below it) is an
**agent**, written between the ``>`` of reaction SMILES and as RXN agents.

RDKit is the oracle for files: components are read by ``convert.to_mol`` from the drawing,
added to an ``rdChemReactions.ChemicalReaction``, and written by RDKit's RXN and SMILES writers.
"""

from __future__ import annotations

import contextlib
import math
import re
from dataclasses import dataclass, field

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from chemsuite.sketch import convert
from chemsuite.sketch.model import Arrow, ArrowKind, Atom, Bond, Group, Plus, Sketch

#: How far off the arrow's line (in bond lengths) a fragment's centre must be to be an agent.
AGENT_OFFSET = 0.5
#: Layout of a reaction that is read in, in bond lengths.
GAP = 0.8  # between a component and a plus sign or the arrow
MIN_ARROW = 2.5
AGENT_LIFT = 0.6  # agents sit this far above the arrow


class ReactionError(ValueError):
    """The drawing is not one reaction (no arrow, several arrows, or an empty side)."""


@dataclass
class Components:
    """Atom ids of each component, in drawing order along the arrow."""

    reactants: list[list[int]] = field(default_factory=list)
    agents: list[list[int]] = field(default_factory=list)
    products: list[list[int]] = field(default_factory=list)


# -- sorting the drawing into components --------------------------------------------------------
def _fragments(sketch: Sketch) -> list[list[int]]:
    seen: set[int] = set()
    out = []
    for aid in sorted(sketch.atoms):
        if aid not in seen:
            part = sketch.fragment(aid)
            seen |= part
            out.append(sorted(part))
    return out


def _centre(sketch: Sketch, atoms: list[int]) -> tuple[float, float]:
    xs = [sketch.atoms[a].x for a in atoms]
    ys = [sketch.atoms[a].y for a in atoms]
    return ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)


def the_arrow(sketch: Sketch) -> Arrow:
    if not sketch.arrows:
        raise ReactionError("The drawing has no reaction arrow.")
    if len(sketch.arrows) > 1:
        raise ReactionError("The drawing has more than one arrow; a reaction has one.")
    return next(iter(sketch.arrows.values()))


def components(sketch: Sketch) -> Components:
    """Reactant, agent and product components of a drawing with one arrow (see the module
    notes). Raises ReactionError without exactly one arrow."""
    arrow = the_arrow(sketch)
    cx, cy = arrow.centre
    length = math.hypot(arrow.x1 - arrow.x0, arrow.y1 - arrow.y0)
    if length < 1e-6:
        raise ReactionError("The reaction arrow has no length.")
    ux, uy = (arrow.x1 - arrow.x0) / length, (arrow.y1 - arrow.y0) / length
    if arrow.kind is ArrowKind.RETROSYNTHETIC:  # target => precursors: read it backwards
        ux, uy = -ux, -uy

    def along(p: tuple[float, float]) -> float:
        return (p[0] - cx) * ux + (p[1] - cy) * uy

    def across(p: tuple[float, float]) -> float:
        return (p[1] - cy) * ux - (p[0] - cx) * uy

    out = Components()
    placed: list[tuple[float, list[int]]] = []
    for part in _fragments(sketch):
        centre = _centre(sketch, part)
        if abs(along(centre)) <= length / 2 and abs(across(centre)) > AGENT_OFFSET:
            out.agents.append(part)
        else:
            placed.append((along(centre), part))
    barriers = sorted([along(p.pos) for p in sketch.pluses.values()] + [0.0])
    grouped: dict[int, list[int]] = {}
    for t, part in placed:
        slot = sum(1 for b in barriers if t > b)  # Ketcher: while (c.x > barriers[j]) ++j
        grouped.setdefault(slot, []).extend(part)
    for slot in sorted(grouped):
        part = sorted(grouped[slot])
        side = out.reactants if along(_centre(sketch, part)) < 0 else out.products
        side.append(part)
    out.agents.sort(key=lambda part: along(_centre(sketch, part)))
    return out


# -- to RDKit, RXN and reaction SMILES --------------------------------------------------------
def _part(sketch: Sketch, atoms: list[int]) -> Sketch:
    keep = set(atoms)
    out = Sketch(next_id=sketch.next_id)
    out.atoms = {aid: a for aid, a in sketch.atoms.items() if aid in keep}
    out.bonds = {bid: b for bid, b in sketch.bonds.items() if b.begin in keep and b.end in keep}
    out.groups = {gid: g for gid, g in sketch.groups.items() if g.atoms <= keep}
    return out


def to_reaction(sketch: Sketch) -> rdChemReactions.ChemicalReaction:
    """The drawing as an RDKit reaction, each component read as the table reads a molecule
    (drawn coordinates and atom maps kept). Raises ReactionError or SketchError."""
    found = components(sketch)
    if not found.reactants and not found.products:
        raise ReactionError("The reaction has no reactants or products.")
    rxn = rdChemReactions.ChemicalReaction()
    for parts, add in (
        (found.reactants, rxn.AddReactantTemplate),
        (found.agents, rxn.AddAgentTemplate),
        (found.products, rxn.AddProductTemplate),
    ):
        for atoms in parts:
            add(convert.to_mol(_part(sketch, atoms)))
    return rxn


def to_smiles(sketch: Sketch) -> str:
    """Reaction SMILES (``reactants>agents>products``)."""
    return rdChemReactions.ReactionToSmiles(to_reaction(sketch))


def to_rxn_block(sketch: Sketch) -> str:
    """An RXN file (V2000; V3000 when a component needs it), agents in their own block.

    Bonds are Kekulé, as in ``convert.to_molblock``: aromatic bond type 4 drops the hydrogen
    of an aromatic [nH], and such a file cannot be read back.
    """
    rxn = _kekule(to_reaction(sketch))
    big = any(
        m.GetNumAtoms() > convert.MAX_V2000_ATOMS
        for m in (*rxn.GetReactants(), *rxn.GetAgents(), *rxn.GetProducts())
    )
    with rdBase.BlockLogs():
        return rdChemReactions.ReactionToRxnBlock(rxn, separateAgents=True, forceV3000=big)


def _kekule(rxn: rdChemReactions.ChemicalReaction) -> rdChemReactions.ChemicalReaction:
    out = rdChemReactions.ChemicalReaction()
    for mols, add in (
        (rxn.GetReactants(), out.AddReactantTemplate),
        (rxn.GetAgents(), out.AddAgentTemplate),
        (rxn.GetProducts(), out.AddProductTemplate),
    ):
        for mol in mols:
            mol = Chem.Mol(mol)
            # keep aromatic bonds that cannot be made alternating
            with rdBase.BlockLogs(), contextlib.suppress(ValueError, RuntimeError):
                Chem.Kekulize(mol, clearAromaticFlags=True)
            add(mol)
    return out


# -- reading ----------------------------------------------------------------------------------
_BRACKETS = re.compile(r"\[[^\]]*\]")


def looks_like_reaction(text: str) -> bool:
    """RXN text, or one line with the two ``>`` of reaction SMILES/SMARTS outside brackets."""
    text = text.strip()
    if text.startswith("$RXN"):
        return True
    if not text or "\n" in text:
        return False
    return _BRACKETS.sub("", text.split()[0]).count(">") == 2


def reaction_from_text(text: str) -> rdChemReactions.ChemicalReaction | None:
    """An RDKit reaction from RXN text, reaction SMILES or reaction SMARTS; None if not one."""
    if not looks_like_reaction(text):
        return None
    text = text.strip()
    with rdBase.BlockLogs():
        try:
            if text.startswith("$RXN"):
                try:
                    return rdChemReactions.ReactionFromRxnBlock(text, sanitize=True)
                except (ValueError, RuntimeError):  # e.g. aromatic bonds RDKit cannot kekulize:
                    # read as drawn; each component is sanitised as far as it can be
                    return rdChemReactions.ReactionFromRxnBlock(text, sanitize=False)
            line = text.split()[0]
            try:
                return rdChemReactions.ReactionFromSmarts(line, useSmiles=True)
            except ValueError:
                return rdChemReactions.ReactionFromSmarts(line)
        except (ValueError, RuntimeError):
            return None


def _template_sketch(mol: Chem.Mol) -> Sketch:
    mol = Chem.Mol(mol)
    with rdBase.BlockLogs():
        try:
            Chem.SanitizeMol(mol)
        except (ValueError, RuntimeError):
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
        else:
            # Reaction templates are read without stereo perception. Read the template again
            # as any molecule is read, so E/Z and centres are perceived the same way (RDKit's
            # legacy AssignStereochemistry drops some bridged-ring centres, e.g. quinuclidine).
            again = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
            if again is not None and again.GetNumAtoms() == mol.GetNumAtoms():
                return convert.from_mol(again)
        # not a molecule RDKit can sanitise: perceive what stereo it can from the template
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return convert.from_mol(mol)


def _width(sketch: Sketch) -> tuple[float, float]:
    x0, y0, x1, y1 = sketch.bounds()
    return x1 - x0, y1 - y0


def _add(target: Sketch, part: Sketch, dx: float, dy: float) -> None:
    """Copy *part* into *target*, moved by (dx, dy), with new ids."""
    ids: dict[int, int] = {}
    for aid, atom in part.atoms.items():
        ids[aid] = target.new_id()
        target.atoms[ids[aid]] = atom.moved(atom.x + dx, atom.y + dy)
    for bond in part.bonds.values():
        target.bonds[target.new_id()] = Bond(
            ids[bond.begin], ids[bond.end], bond.order, bond.stereo, bond.topology
        )
    for group in part.groups.values():
        target.groups[target.new_id()] = Group(
            group.label, frozenset(ids[a] for a in group.atoms), ids[group.attach],
            group.expanded,
        )  # fmt: skip


def from_reaction(rxn: rdChemReactions.ChemicalReaction) -> Sketch:
    """A drawing of *rxn*: reactants, plus signs, the arrow (agents above it), products, left to
    right on one line, each component in its own coordinates or a fresh 2D layout."""
    out = Sketch()
    reactants = [_template_sketch(m) for m in rxn.GetReactants()]
    agents = [_template_sketch(m) for m in rxn.GetAgents()]
    products = [_template_sketch(m) for m in rxn.GetProducts()]
    x = 0.0

    def place(part: Sketch) -> None:
        nonlocal x
        x0, y0, x1, y1 = part.bounds()
        _add(out, part, x - x0, -(y0 + y1) / 2)
        x += x1 - x0

    def side(parts: list[Sketch]) -> None:
        nonlocal x
        for i, part in enumerate(parts):
            if i:
                out.pluses[out.new_id()] = Plus(x + GAP, 0.0)
                x += 2 * GAP
            place(part)

    side(reactants)
    width = max([MIN_ARROW] + [_width(a)[0] + 2 * GAP for a in agents])
    start = x + GAP if reactants else x
    out.arrows[out.new_id()] = Arrow(start, 0.0, start + width, 0.0)
    lift = AGENT_LIFT
    for agent in agents:
        x0, y0, x1, y1 = agent.bounds()
        _add(out, agent, start + width / 2 - (x0 + x1) / 2, lift - y0)
        lift += (y1 - y0) + GAP
    x = start + width + GAP
    side(products)
    return out


def from_text(text: str) -> Sketch | None:
    """A drawing of a reaction given as RXN text or reaction SMILES/SMARTS; None if not one."""
    rxn = reaction_from_text(text)
    return None if rxn is None else from_reaction(rxn)


def mapped_atoms(sketch: Sketch) -> dict[int, list[int]]:
    """Map number -> the atoms carrying it."""
    out: dict[int, list[int]] = {}
    for aid, atom in sketch.atoms.items():
        if atom.mapping:
            out.setdefault(atom.mapping, []).append(aid)
    return out


__all__ = [
    "Atom",
    "Components",
    "ReactionError",
    "components",
    "from_reaction",
    "from_text",
    "looks_like_reaction",
    "mapped_atoms",
    "reaction_from_text",
    "the_arrow",
    "to_reaction",
    "to_rxn_block",
    "to_smiles",
]
