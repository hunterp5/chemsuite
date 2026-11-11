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
"""Selecting atoms with a PyMOL-like language; a selection is a boolean mask over atoms.

Words (case does not matter; several values join with ``+`` or ``,``)::

    chain A+B           resi 10-50+60  resi 52A    resn HIS+ASP    name CA+C*   elem Fe
    protein  nucleic  polymer  ligand  water  ion  additive  glycan  hetero  organic
    backbone  sidechain  hydrogen  helix  sheet  altloc  all  none
    byres S    bychain S    within 5 of S    around 5 of S (within, without S itself)
    S within 5 of T (S near T)    S around 5 of T    S and S    S or S    not S    ( S )
    $name (a saved selection)

``organic`` is ligands, additives and glycans; ``hetero`` is everything not in a polymer.
``name`` and ``resn`` accept ``*`` and ``?`` wildcards.
"""

from __future__ import annotations

import fnmatch
import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .model import (
    ADDITIVE,
    AMINO,
    BACKBONE,
    GLYCAN,
    HELIX,
    ION,
    LIGAND,
    NUCLEIC,
    NUCLEIC_BACKBONE,
    STRAND,
    WATER,
    Structure,
)


class SelectionError(ValueError):
    """The text is not a selection; ``position`` is where reading stopped."""

    def __init__(self, message: str, position: int) -> None:
        super().__init__(message)
        self.position = position


@dataclass(frozen=True)
class _Token:
    text: str
    position: int


_TOKEN = re.compile(r"\s*(\(|\)|[^\s()]+)")
_RESI = re.compile(r"^(-?\d+)([A-Za-z]?)(?:-(-?\d+)([A-Za-z]?))?$")

Mask = np.ndarray
_Eval = Callable[["_Context"], Mask]


@dataclass
class _Context:
    structure: Structure
    xyz: np.ndarray
    named: Mapping[str, Mask]


def select(
    structure: Structure,
    text: str,
    named: Mapping[str, Mask] | None = None,
    model: int = 0,
) -> Mask:
    """The atoms ``text`` selects. Raises :class:`SelectionError` if it cannot be read."""
    return compile_selection(text)(_Context(structure, structure.xyz(model), named or {}))


def compile_selection(text: str) -> _Eval:
    tokens = [_Token(m.group(1), m.start(1)) for m in _TOKEN.finditer(text)]
    if not tokens:
        return lambda c: np.zeros(c.structure.n_atoms, dtype=bool)
    parser = _Parser(tokens, len(text))
    tree = parser.expression()
    if parser.at < len(tokens):
        token = tokens[parser.at]
        raise SelectionError(f"unexpected '{token.text}'", token.position)
    return tree


def _values(text: str) -> list[str]:
    return [v for v in re.split(r"[+,]", text) if v]


def _matches(column: np.ndarray, patterns: list[str], upper: bool = True) -> Mask:
    if upper:
        patterns = [p.upper() for p in patterns]
        column = np.char.upper(column)
    plain = [p for p in patterns if not any(ch in p for ch in "*?[")]
    mask = np.isin(column, plain)
    wild = [p for p in patterns if p not in plain]
    if wild:
        unique = np.unique(column)
        hits = [u for u in unique if any(fnmatch.fnmatchcase(str(u), p) for p in wild)]
        mask |= np.isin(column, hits)
    return mask


def _per_atom(residue_mask: Mask) -> _Eval:
    return lambda c: residue_mask[c.structure.residue_of]


def _kinds(*kinds: int) -> _Eval:
    return lambda c: np.isin(c.structure.atom_kind, kinds)


def _names_in(kind: int, names: tuple[str, ...]) -> Callable[[_Context], Mask]:
    return lambda c: (c.structure.atom_kind == kind) & np.isin(c.structure.name, names)


def _backbone(c: _Context) -> Mask:
    return _names_in(AMINO, BACKBONE)(c) | _names_in(NUCLEIC, NUCLEIC_BACKBONE)(c)


def _sidechain(c: _Context) -> Mask:
    return np.isin(c.structure.atom_kind, (AMINO, NUCLEIC)) & ~_backbone(c)


def _ss(kind: int) -> _Eval:
    return lambda c: c.structure.ss[c.structure.residue_of] == kind


_KEYWORDS: dict[str, _Eval] = {
    "all": lambda c: np.ones(c.structure.n_atoms, dtype=bool),
    "none": lambda c: np.zeros(c.structure.n_atoms, dtype=bool),
    "protein": _kinds(AMINO),
    "nucleic": _kinds(NUCLEIC),
    "polymer": _kinds(AMINO, NUCLEIC),
    "ligand": _kinds(LIGAND),
    "water": _kinds(WATER),
    "solvent": _kinds(WATER),
    "ion": _kinds(ION),
    "additive": _kinds(ADDITIVE),
    "glycan": _kinds(GLYCAN),
    "organic": _kinds(LIGAND, ADDITIVE, GLYCAN),
    "hetero": lambda c: ~np.isin(c.structure.atom_kind, (AMINO, NUCLEIC)),
    "backbone": _backbone,
    "bb": _backbone,
    "sidechain": _sidechain,
    "sc": _sidechain,
    "hydrogen": lambda c: c.structure.element == 1,
    "h": lambda c: c.structure.element == 1,
    "helix": _ss(HELIX),
    "sheet": _ss(STRAND),
    "strand": _ss(STRAND),
    "altloc": lambda c: c.structure.altloc != "",
}


class _Parser:
    def __init__(self, tokens: list[_Token], end: int) -> None:
        self.tokens = tokens
        self.at = 0
        self.end = end

    # helpers
    def peek(self) -> str | None:
        return self.tokens[self.at].text.lower() if self.at < len(self.tokens) else None

    def take(self, what: str = "a value") -> _Token:
        if self.at >= len(self.tokens):
            raise SelectionError(f"expected {what}", self.end)
        token = self.tokens[self.at]
        self.at += 1
        return token

    def expect(self, word: str) -> None:
        token = self.take(f"'{word}'")
        if token.text.lower() != word:
            raise SelectionError(f"expected '{word}', found '{token.text}'", token.position)

    # grammar: or < and < not < primary
    def expression(self) -> _Eval:
        left = self.conjunction()
        while self.peek() in ("or", "|"):
            self.at += 1
            right = self.conjunction()
            left = (lambda a, b: lambda c: a(c) | b(c))(left, right)
        return left

    def conjunction(self) -> _Eval:
        left = self.negation()
        while self.peek() in ("and", "&"):
            self.at += 1
            right = self.negation()
            left = (lambda a, b: lambda c: a(c) & b(c))(left, right)
        return left

    def negation(self) -> _Eval:
        if self.peek() in ("not", "!"):
            self.at += 1
            inner = self.negation()
            return lambda c: ~inner(c)
        return self.near()

    def near(self) -> _Eval:
        """``S within 5 of T``: PyMOL's infix form."""
        left = self.primary()
        while self.peek() in ("within", "around"):
            exclusive = self.take().text.lower() == "around"
            right = self.distance(exclusive)
            left = (lambda a, b: lambda c: a(c) & b(c))(left, right)
        return left

    def primary(self) -> _Eval:
        token = self.take("a selection")
        word = token.text.lower()
        if word == "(":
            inner = self.expression()
            self.expect(")")
            return inner
        if word in _KEYWORDS:
            return _KEYWORDS[word]
        if word == "chain":
            values = _values(self.take("a chain").text)
            return lambda c: np.isin(c.structure.chain_name, values)[c.structure.atom_chain]
        if word == "resn":
            values = _values(self.take("a residue name").text)
            return lambda c: _matches(c.structure.res_name, values)[c.structure.residue_of]
        if word == "name":
            values = _values(self.take("an atom name").text)
            return lambda c: _matches(c.structure.name, values)
        if word in ("elem", "element"):
            values = [v.capitalize() for v in _values(self.take("an element").text)]
            return lambda c: _elements(c, values)
        if word in ("resi", "resid"):
            return self.residue_numbers(self.take("a residue number"))
        if word in ("byres", "bychain"):
            inner = self.primary()
            return _byres(inner) if word == "byres" else _bychain(inner)
        if word in ("within", "around"):
            return self.distance(word == "around")
        if word.startswith("$") and len(word) > 1:
            key = token.text[1:]
            return lambda c: _named(c, key, token.position)
        raise SelectionError(f"unknown word '{token.text}'", token.position)

    def residue_numbers(self, token: _Token) -> _Eval:
        ranges = []
        for value in _values(token.text):
            match = _RESI.match(value)
            if not match:
                raise SelectionError(f"'{value}' is not a residue number", token.position)
            first, first_code, last, _ = match.groups()
            ranges.append((int(first), first_code.upper(), None if last is None else int(last)))

        def evaluate(c: _Context) -> Mask:
            seq, icode = c.structure.res_seq, np.char.upper(c.structure.icode)
            mask = np.zeros(c.structure.n_residues, dtype=bool)
            for first, code, last in ranges:
                if last is None:
                    mask |= (seq == first) & (icode == code)
                else:
                    mask |= (seq >= min(first, last)) & (seq <= max(first, last))
            return mask[c.structure.residue_of]

        return evaluate

    def distance(self, exclusive: bool) -> _Eval:
        token = self.take("a distance")
        try:
            radius = float(token.text)
        except ValueError:
            raise SelectionError(f"'{token.text}' is not a distance", token.position) from None
        if self.peek() == "of":
            self.at += 1
        inner = self.primary()

        def evaluate(c: _Context) -> Mask:
            target = inner(c)
            near = within(c.xyz, target, radius)
            return near & ~target if exclusive else near

        return evaluate


def within(xyz: np.ndarray, target: Mask, radius: float) -> Mask:
    """Atoms closer than ``radius`` to any atom in ``target`` (which are themselves in)."""
    if not target.any():
        return np.zeros(len(xyz), dtype=bool)
    distance, _ = cKDTree(xyz[target]).query(xyz, distance_upper_bound=radius)
    return np.isfinite(distance)


def _byres(inner: _Eval) -> _Eval:
    def evaluate(c: _Context) -> Mask:
        residues = np.zeros(c.structure.n_residues, dtype=bool)
        residues[c.structure.residue_of[inner(c)]] = True
        return residues[c.structure.residue_of]

    return evaluate


def _bychain(inner: _Eval) -> _Eval:
    def evaluate(c: _Context) -> Mask:
        chains = np.zeros(c.structure.n_chains, dtype=bool)
        chains[c.structure.atom_chain[inner(c)]] = True
        return chains[c.structure.atom_chain]

    return evaluate


def _elements(c: _Context, symbols: list[str]) -> Mask:
    numbers = [gemmi.Element(symbol).atomic_number for symbol in symbols]
    return np.isin(c.structure.element, [z for z in numbers if z > 0])


def _named(c: _Context, key: str, position: int) -> Mask:
    mask = c.named.get(key)
    if mask is None:
        raise SelectionError(f"no saved selection '{key}'", position)
    if len(mask) != c.structure.n_atoms:
        raise SelectionError(f"saved selection '{key}' is for other atoms", position)
    return np.asarray(mask, dtype=bool)


def residues_in(structure: Structure, mask: Mask) -> np.ndarray:
    """Indices of the residues with at least one selected atom."""
    return np.unique(structure.residue_of[mask])
