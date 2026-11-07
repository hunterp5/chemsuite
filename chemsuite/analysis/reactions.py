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

"""Using reactions (REACTIONS_PLAN R4): search them, split them, and apply them.

* **Reaction search.** A reaction query (reaction SMILES or SMARTS, or drawn) matches a reaction
  when each query reactant is a substructure of one of its reactants, each query product of one
  of its products, and each query agent of one of its agents (RDKit
  ``HasReactionSubstructMatch``). A structure query matches a reaction with that substructure in
  a chosen role: a reactant, a product, an agent, or any component.
* **Split Reaction.** The reactants, agents and products of each reaction as structures: one
  structure per role (components joined, as reaction SMILES writes them) or one per component.
* **Enumerate Products.** A reaction template (reaction SMARTS with atom maps, drawn or typed)
  applied to building blocks, with RDKit ``RunReactants``: row by row, reactant *i* from column
  *i* of the same row; or every combination of the building blocks (library enumeration). Each
  distinct product set is one result: a template with several products gives them joined
  with ``.``, and a building block that matches at several sites gives several results.

Every per-row function works in chunks (in worker processes with an executor) and, when
cancelled, raises :class:`Cancelled` carrying the results for the rows it finished.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from chemsuite.analysis.structure_search import QueryError, StructureQuery, parse_query
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_text
from chemsuite.core.reactions import looks_like_reaction, reaction_from_blob, reaction_from_text

#: A cell as sent to a worker: RDKit binary, text not parsed yet, or None when blank.
Cell = bytes | str | None
Progress = Callable[[int, int], None]

CHUNK_ROWS = 500
INLINE_CHUNK_ROWS = 100
POOL_MIN_ROWS = 2000

#: query role -> menu title
ROLES = {
    "reaction": "Reaction contains",
    "reactant": "A reactant contains",
    "product": "A product contains",
    "agent": "An agent contains",
    "any": "Any component contains",
}


# -- loading cells ----------------------------------------------------------------------------


def _templates(rxn: rdChemReactions.ChemicalReaction) -> list[Chem.Mol]:
    return [*rxn.GetReactants(), *rxn.GetAgents(), *rxn.GetProducts()]


def _role_templates(rxn: rdChemReactions.ChemicalReaction, role: str) -> list[Chem.Mol]:
    if role == "reactant":
        return list(rxn.GetReactants())
    if role == "product":
        return list(rxn.GetProducts())
    if role == "agent":
        return list(rxn.GetAgents())
    return _templates(rxn)


def load_reaction(cell: Cell) -> rdChemReactions.ChemicalReaction | None:
    """The reaction in a cell, with ring information on every component so substructure queries
    with ring features work (reaction components are read without it)."""
    if not cell:
        return None
    with rdBase.BlockLogs():
        rxn = reaction_from_blob(cell) if isinstance(cell, bytes) else reaction_from_text(cell)
    if rxn is None:
        return None
    for template in _templates(rxn):
        template.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(template)
    return rxn


def load_mol(cell: Cell) -> Chem.Mol | None:
    if not cell:
        return None
    if isinstance(cell, bytes):
        return Chem.Mol(cell)
    return mol_from_text(cell)


def _chunked(total: int, executor: Executor | None) -> Iterator[slice]:
    size = CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS
    for start in range(0, total, size):
        yield slice(start, min(start + size, total))


# -- reaction search ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReactionQuery:
    """A reaction query (``role="reaction"``), or a structure to find in one role."""

    text: str
    role: str = "reaction"

    def describe(self) -> str:
        text = self.text.strip()
        if "\n" in text:
            text = "drawn query"
        elif len(text) > 40:
            text = text[:39] + "…"
        return f"{ROLES.get(self.role, self.role).lower()} {text}"


def parse_reaction_query(
    query: ReactionQuery,
) -> rdChemReactions.ChemicalReaction | Chem.Mol:
    """The query reaction (role "reaction") or query structure; raises :class:`QueryError`."""
    text = query.text
    if not text.strip():
        raise QueryError("Enter a query.")
    if query.role not in ROLES:
        raise QueryError(f"Unknown role {query.role!r}.")
    if query.role == "reaction":
        if not looks_like_reaction(text):
            raise QueryError("Enter a reaction, such as C(=O)O.N>>C(=O)N, or draw one.")
        rxn = _query_reaction(text)
        if rxn is None:
            raise QueryError("Not a valid reaction SMILES or SMARTS.")
        return rxn
    if looks_like_reaction(text):
        raise QueryError("This is a reaction: choose “Reaction contains” to search with it.")
    return parse_query(StructureQuery(text, "substructure")).mol


def _query_reaction(text: str) -> rdChemReactions.ChemicalReaction | None:
    """Reaction SMARTS first, so query features (``[#6]``, ``[R]``) stay queries; reaction
    SMILES (``useSmiles``) for text SMARTS cannot read; an RXN block as drawn."""
    text = text.strip() if not text.lstrip().startswith("$RXN") else text
    with rdBase.BlockLogs():
        if text.lstrip().startswith("$RXN"):
            return reaction_from_text(text)
        try:
            return rdChemReactions.ReactionFromSmarts(text.split()[0])
        except ValueError:
            pass
        try:
            return rdChemReactions.ReactionFromSmarts(text.split()[0], useSmiles=True)
        except ValueError:
            return None


def reaction_matches(
    rxn: rdChemReactions.ChemicalReaction,
    query: rdChemReactions.ChemicalReaction | Chem.Mol,
    role: str,
) -> bool:
    with rdBase.BlockLogs():
        if role == "reaction":
            include_agents = query.GetNumAgentTemplates() > 0
            return rdChemReactions.HasReactionSubstructMatch(
                rxn, query, includeAgents=include_agents
            )
        return any(t.HasSubstructMatch(query) for t in _role_templates(rxn, role))


def search_chunk(args: tuple[str, str, Sequence[Cell]]) -> np.ndarray:
    """Match one chunk of cells (worker-process entry point)."""
    text, role, cells = args
    query = parse_reaction_query(ReactionQuery(text, role))
    out = np.zeros(len(cells), dtype=bool)
    for i, cell in enumerate(cells):
        rxn = load_reaction(cell)
        out[i] = rxn is not None and reaction_matches(rxn, query, role)
    return out


def search_reactions(
    query: ReactionQuery,
    cells: Sequence[Cell],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """``bool`` per cell: whether it matches *query*. Cancelled: ``partial`` holds the matches
    for the first rows."""
    parse_reaction_query(query)  # fail before starting
    total = len(cells)
    chunks = ((query.text, query.role, list(cells[s])) for s in _chunked(total, executor))
    parts: list[np.ndarray] = []
    done = 0
    try:
        for part in ordered_map(search_chunk, chunks, executor, cancel=cancel):
            parts.append(part)
            done += len(part)
            if progress is not None:
                progress(done, total)
    except Cancelled:
        raise Cancelled(np.concatenate(parts) if parts else None) from None
    return np.concatenate(parts) if parts else np.zeros(0, dtype=bool)


# -- split reaction ---------------------------------------------------------------------------

SPLIT_ROLES = ("reactant", "agent", "product")


def _component_smiles(mol: Chem.Mol, keep_maps: bool) -> str:
    mol = Chem.Mol(mol)
    if not keep_maps:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    with rdBase.BlockLogs():
        try:
            Chem.SanitizeMol(mol)
        except (ValueError, RuntimeError):
            mol.UpdatePropertyCache(strict=False)
        try:
            return Chem.MolToSmiles(mol)
        except (ValueError, RuntimeError):
            return ""


def split_reaction(
    rxn: rdChemReactions.ChemicalReaction, keep_maps: bool = False
) -> dict[str, list[str]]:
    """SMILES of each component by role (``reactant``, ``agent``, ``product``), in order. Atom
    maps are removed unless *keep_maps*: in a structure column they would stop identical
    molecules from comparing equal."""
    return {
        role: [_component_smiles(m, keep_maps) for m in _role_templates(rxn, role)]
        for role in SPLIT_ROLES
    }


def split_chunk(args: tuple[bool, Sequence[Cell]]) -> list[dict[str, list[str]] | None]:
    keep_maps, cells = args
    out: list[dict[str, list[str]] | None] = []
    for cell in cells:
        with rdBase.BlockLogs():
            rxn = (
                reaction_from_blob(cell)
                if isinstance(cell, bytes)
                else (reaction_from_text(cell) if cell else None)
            )
        out.append(None if rxn is None else split_reaction(rxn, keep_maps))
    return out


def split_cells(
    cells: Sequence[Cell],
    *,
    keep_maps: bool = False,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> list[dict[str, list[str]] | None]:
    """:func:`split_reaction` per cell (None for a blank or invalid one)."""
    total = len(cells)
    chunks = ((keep_maps, list(cells[s])) for s in _chunked(total, executor))
    out: list[dict[str, list[str]] | None] = []
    try:
        for part in ordered_map(split_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(len(out), total)
    except Cancelled:
        raise Cancelled(out or None) from None
    return out


# -- enumerate products -------------------------------------------------------------------------


class TemplateError(ValueError):
    """The text is not a reaction template that can be applied."""


def parse_template(text: str) -> rdChemReactions.ChemicalReaction:
    """A reaction template from reaction SMARTS (or SMILES, or a drawn RXN block), checked:
    it needs at least one reactant and one product, and atom maps to say which reactant atoms
    become which product atoms."""
    if not text.strip():
        raise TemplateError("Enter a reaction template, or draw one.")
    if not looks_like_reaction(text):
        raise TemplateError("A template is a reaction: reactants >> products.")
    rxn = _query_reaction(text)
    if rxn is None:
        raise TemplateError("Not a valid reaction SMARTS.")
    if not rxn.GetNumReactantTemplates() or not rxn.GetNumProductTemplates():
        raise TemplateError("A template needs at least one reactant and one product.")
    mapped = any(a.GetAtomMapNum() for m in rxn.GetReactants() for a in m.GetAtoms())
    if not mapped:
        raise TemplateError(
            "Map the atoms the reaction keeps (the Atom Map tool, or :1, :2 … in SMARTS), so "
            "the product can be built from the building blocks."
        )
    with rdBase.BlockLogs():
        rxn.Initialize()
        _warnings, errors = rxn.Validate()
    if errors:
        raise TemplateError("RDKit cannot apply this template (check its atom maps).")
    return rxn


def apply_template(
    rxn: rdChemReactions.ChemicalReaction, reactants: Sequence[Chem.Mol]
) -> list[str]:
    """The distinct products of *reactants* (canonical SMILES; several products of one outcome
    joined with ``.``). Outcomes that are not valid molecules are dropped."""
    out: list[str] = []
    seen: set[str] = set()
    with rdBase.BlockLogs():
        try:
            outcomes = rxn.RunReactants(tuple(reactants))
        except (ValueError, RuntimeError):
            return out
        for outcome in outcomes:
            parts = []
            for product in outcome:
                try:
                    Chem.SanitizeMol(product)
                    parts.append(Chem.MolToSmiles(product))
                except (ValueError, RuntimeError):
                    parts = []
                    break
            if not parts:
                continue
            smiles = ".".join(parts)
            canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if smiles else ""
            if canonical and canonical not in seen:
                seen.add(canonical)
                out.append(canonical)
    return out


def enumerate_rows_chunk(args: tuple[str, Sequence[Sequence[Cell]]]) -> list[list[str]]:
    """Products for each row's reactants (worker-process entry point)."""
    text, rows = args
    rxn = parse_template(text)
    out: list[list[str]] = []
    for row in rows:
        mols = [load_mol(cell) for cell in row]
        out.append([] if any(m is None for m in mols) else apply_template(rxn, mols))
    return out


def enumerate_rows(
    template: str,
    columns: Sequence[Sequence[Cell]],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> list[list[str]]:
    """Row by row: reactant *i* of each row from ``columns[i]``. One list of products per row
    (empty when a building block is blank or the template does not apply)."""
    rxn = parse_template(template)
    if len(columns) != rxn.GetNumReactantTemplates():
        raise TemplateError(
            f"The template has {rxn.GetNumReactantTemplates()} reactant(s); "
            f"choose a building block column for each."
        )
    rows = list(zip(*columns, strict=True))
    total = len(rows)
    chunks = ((template, rows[s]) for s in _chunked(total, executor))
    out: list[list[str]] = []
    try:
        for part in ordered_map(enumerate_rows_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(len(out), total)
    except Cancelled:
        raise Cancelled(out or None) from None
    return out


@dataclass
class Library:
    """Products of every combination of building blocks: one entry per product, with the index
    of each building block it came from (``blocks[i][k]`` is the index into list *i*)."""

    blocks: list[list[int]] = field(default_factory=list)
    products: list[str] = field(default_factory=list)
    combinations: int = 0  # combinations tried
    failed: int = 0  # combinations that gave no product

    def extend(self, other: Library) -> None:
        for mine, theirs in zip(self.blocks, other.blocks, strict=True):
            mine.extend(theirs)
        self.products.extend(other.products)
        self.combinations += other.combinations
        self.failed += other.failed


def combinations_count(sizes: Sequence[int]) -> int:
    return math.prod(sizes) if sizes else 0


def enumerate_chunk(
    args: tuple[str, Sequence[Sequence[Cell]], int, int],
) -> Library:
    """Combinations ``start`` to ``stop`` (flat indices, last list fastest)."""
    text, blocks, start, stop = args
    rxn = parse_template(text)
    mols = [[load_mol(c) for c in block] for block in blocks]
    sizes = [len(b) for b in blocks]
    out = Library([[] for _ in blocks])
    for flat in range(start, stop):
        index = np.unravel_index(flat, sizes)
        reactants = [mols[i][int(k)] for i, k in enumerate(index)]
        out.combinations += 1
        products = [] if any(m is None for m in reactants) else apply_template(rxn, reactants)
        if not products:
            out.failed += 1
        for product in products:
            for i, k in enumerate(index):
                out.blocks[i].append(int(k))
            out.products.append(product)
    return out


def enumerate_library(
    template: str,
    blocks: Sequence[Sequence[Cell]],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Library:
    """Every combination of one building block from each list (``blocks[i]`` for reactant
    *i*). Cancelled: ``partial`` holds the products of the combinations finished."""
    rxn = parse_template(template)
    if len(blocks) != rxn.GetNumReactantTemplates():
        raise TemplateError(
            f"The template has {rxn.GetNumReactantTemplates()} reactant(s); "
            f"choose building blocks for each."
        )
    total = combinations_count([len(b) for b in blocks])
    lists = [list(b) for b in blocks]
    chunks = ((template, lists, s.start, s.stop) for s in _chunked(total, executor))
    out = Library([[] for _ in blocks])
    try:
        for part in ordered_map(enumerate_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(out.combinations, total)
    except Cancelled:
        raise Cancelled(out if out.combinations else None) from None
    return out


__all__ = [
    "ROLES",
    "Library",
    "QueryError",
    "ReactionQuery",
    "TemplateError",
    "apply_template",
    "combinations_count",
    "enumerate_library",
    "enumerate_rows",
    "parse_reaction_query",
    "parse_template",
    "reaction_matches",
    "search_reactions",
    "split_cells",
    "split_reaction",
]
