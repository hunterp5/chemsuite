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

"""Text search over a table: quick search, find and find/replace.

Everything matches the cells' displayed text, the same text the user sees and copies.
Structure columns are skipped: their text is SMILES or a molfile, which a text search would
match in surprising ways (substructure search is the tool for those).
"""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from .columns import REACTION, STRUCTURE

if TYPE_CHECKING:
    from .table import Table

SearchMode = Literal["contains", "whole", "regex"]


@dataclass(frozen=True)
class Query:
    """What to look for. ``whole`` matches the entire cell; ``regex`` uses Python syntax."""

    text: str
    mode: SearchMode = "contains"
    case_sensitive: bool = False

    def pattern(self) -> re.Pattern[str]:
        """The compiled pattern; raises :class:`re.error` for an invalid regular expression."""
        flags = 0 if self.case_sensitive else re.IGNORECASE
        if self.mode == "regex":
            return re.compile(self.text, flags)
        escaped = re.escape(self.text)
        return re.compile(f"^{escaped}$" if self.mode == "whole" else escaped, flags)

    def predicate(self) -> Callable[[str], bool]:
        pattern = self.pattern()
        if self.mode == "whole":
            return lambda text: pattern.fullmatch(text) is not None
        return lambda text: pattern.search(text) is not None

    def replace(self, text: str, replacement: str) -> str:
        """*text* with every match replaced.

        In regex mode the replacement may use the groups: ``\\1`` or ``$1``, and
        ``\\g<name>`` or ``${name}`` (``$$`` is a dollar sign). Elsewhere it is literal.
        """
        pattern = self.pattern()
        if self.mode == "whole":
            return replacement if pattern.fullmatch(text) else text
        if self.mode == "regex":
            return pattern.sub(regex_template(replacement), text)
        return pattern.sub(lambda _m: replacement, text)  # literal: no escapes in replacement

    def replacement_problem(self, replacement: str) -> str:
        """Why *replacement* cannot be used with this query, or ``""``.

        Checked before anything is replaced, so a bad group reference never fails halfway: the
        template is tried on a pattern with the same groups that matches the empty string.
        """
        if self.mode != "regex":
            return ""
        pattern = self.pattern()
        names = {index: name for name, index in pattern.groupindex.items()}
        groups = "".join(
            f"(?P<{names[i]}>)" if i in names else "()" for i in range(1, pattern.groups + 1)
        )
        try:
            re.compile(groups).sub(regex_template(replacement), "")
        except (re.error, IndexError) as err:
            reason = str(err).split(" at position")[0]  # positions are in the rewritten form
            return f"Invalid replacement: {reason}"
        return ""


_DOLLAR = re.compile(r"\$(?:(?P<dollar>\$)|(?P<number>\d+)|\{(?P<name>[A-Za-z_]\w*|\d+)\})")


def regex_template(replacement: str) -> str:
    """*replacement* with ``$1`` and ``${name}`` written as Python's ``\\g<1>``, ``\\g<name>``."""

    def swap(match: re.Match[str]) -> str:
        if match["dollar"]:
            return "$"
        return "\\g<" + (match["number"] or match["name"]) + ">"

    return _DOLLAR.sub(swap, replacement)


# -- field syntax (the quick search bar) --------------------------------------------------------
Op = Literal["contains", "equals", "<", "<=", ">", ">=", "!="]

#: ``field`` then an operator then the value; the field may be quoted to hold spaces
_FIELD = re.compile(r'^(?P<field>"[^"]*"|[^"<>=!:\s]+)(?P<op>:|<=|>=|!=|<|>|=)(?P<value>.*)$')


@dataclass(frozen=True)
class Term:
    """One condition of a quick search: text in a cell, or a number compared.

    *column* None means the bar's own choice (all columns, or the one it is set to).
    """

    text: str
    column: str | None = None
    op: Op = "contains"
    number: float | None = None
    negate: bool = False


def tokens(text: str) -> list[str]:
    """*text* split on spaces outside double quotes (an unclosed quote runs to the end)."""
    out: list[str] = []
    current: list[str] = []
    quoted = False
    for char in text:
        if char == '"':
            quoted = not quoted
            current.append(char)
        elif char.isspace() and not quoted:
            if current:
                out.append("".join(current))
                current = []
        else:
            current.append(char)
    if current:
        out.append("".join(current))
    return out


def _unquote(text: str) -> str:
    if len(text) >= 2 and text[0] == text[-1] == '"':
        return text[1:-1]
    return text[1:] if text.startswith('"') else text  # a quote still being typed


def parse_terms(text: str, columns: Sequence[str]) -> tuple[Term, ...]:
    """The quick-search *text* as terms, all of which a row must meet.

    - ``word`` – a cell contains it; several words must all be found (anywhere in the row);
    - ``"exact phrase"`` – a cell contains the phrase, spaces included;
    - ``ID:CPD-12`` – that column contains the text (``ID:"two words"`` for spaces);
    - ``pIC50>7``, ``>=``, ``<``, ``<=``, ``=``, ``!=`` – that column's number compared;
      ``Series=A`` with a non-number is the whole cell, ignoring case;
    - a leading ``-`` excludes rows that meet the term (``-Series:B``).

    Field names match a column ignoring case; quote names with spaces (``"Mol Wt"<400``). A
    field that names no column is searched for as plain text, so ``http://`` or ``10:30`` still
    work. Terms with nothing after the operator are skipped while they are being typed.
    """
    by_name = {name.casefold(): name for name in columns}
    terms: list[Term] = []
    for token in tokens(text):
        negate = token.startswith("-") and len(token) > 1
        body = token[1:] if negate else token
        found = _FIELD.match(body)
        column = by_name.get(_unquote(found["field"]).casefold()) if found else None
        if found is None or column is None:
            plain = _unquote(body)
            if plain:
                terms.append(Term(plain, negate=negate))
            continue
        op, value = found["op"], _unquote(found["value"])
        if not value:
            continue
        number = _number(value)
        if op == ":":
            terms.append(Term(value, column, "contains", negate=negate))
        elif op in ("<", "<=", ">", ">=") and number is not None:
            terms.append(Term(value, column, op, number, negate))
        elif op in ("=", "!=") and number is not None:
            terms.append(Term(value, column, op if op == "!=" else "=", number, negate))
        elif op == "=":
            terms.append(Term(value, column, "equals", negate=negate))
        elif op == "!=":
            terms.append(Term(value, column, "equals", negate=not negate))
        else:  # a comparison with no number: look for the text as typed
            terms.append(Term(_unquote(body), negate=negate))
    return tuple(terms)


def terms_mask(
    table: Table,
    terms: Sequence[Term],
    *,
    column: str | None = None,
    case_sensitive: bool = False,
) -> np.ndarray:
    """``bool`` per row (storage order): rows that meet every term.

    *column* is where terms without a field look (None: every searchable column).
    """
    from .columns import NumberColumn, to_number

    keep = np.ones(len(table), dtype=np.bool_)
    for term in terms:
        target = term.column or column
        if target is not None and target not in table:
            hit = np.zeros(len(table), dtype=np.bool_)
        elif term.number is not None:
            values = table.column(target)
            if not isinstance(values, NumberColumn):
                values = to_number(values, force=True)[0]
            numbers = values.values()
            with np.errstate(invalid="ignore"):
                hit = {
                    "<": numbers < term.number,
                    "<=": numbers <= term.number,
                    ">": numbers > term.number,
                    ">=": numbers >= term.number,
                    "=": np.isclose(numbers, term.number),
                    "!=": ~np.isclose(numbers, term.number) & ~np.isnan(numbers),
                }[term.op]
        else:
            mode: SearchMode = "whole" if term.op == "equals" else "contains"
            query = Query(term.text, mode, case_sensitive)
            names = None if target is None else [target]
            hit = match_mask(table, query, names).any(axis=1)
        keep &= ~hit if term.negate else hit
    return keep


def _number(text: str) -> float | None:
    try:
        value = float(text.replace(",", ""))
    except ValueError:
        return None
    return value if np.isfinite(value) else None


def searchable_columns(table: Table, columns: Sequence[str] | None = None) -> list[str]:
    names = table.column_names if columns is None else [c for c in columns if c in table]
    return [name for name in names if table.kind(name) not in (STRUCTURE, REACTION)]


def match_mask(table: Table, query: Query, columns: Sequence[str] | None = None) -> np.ndarray:
    """``bool[rows, columns]`` in storage order for the searchable *columns*."""
    names = searchable_columns(table, columns)
    if not query.text or not names:
        return np.zeros((len(table), len(names)), dtype=np.bool_)
    predicate = query.predicate()
    return np.column_stack([table.column(name).match_text(predicate) for name in names])


@dataclass(frozen=True)
class Matches:
    """Matching cells in display order (row by row, left to right)."""

    oids: np.ndarray  # int64, one per match
    columns: tuple[str, ...]  # the searched columns
    column_index: np.ndarray  # index into *columns*, one per match

    def __len__(self) -> int:
        return len(self.oids)

    def cell(self, i: int) -> tuple[int, str]:
        return int(self.oids[i]), self.columns[int(self.column_index[i])]


#: some cells of a table: a column name to the oids of its cells (a selection)
Cells = Mapping[str, "np.ndarray | Sequence[int]"]


def find(
    table: Table,
    query: Query,
    columns: Sequence[str] | None = None,
    *,
    shown_only: bool = True,
    within: Cells | None = None,
) -> Matches:
    """All matching cells, in the order the table shows them (hidden rows last if included).

    *within* limits the search to those cells; columns it does not name are not searched.
    """
    if within is not None:
        columns = [name for name in (columns or table.column_names) if name in within]
    names = tuple(searchable_columns(table, columns))
    mask = match_mask(table, query, names)
    if within is not None:
        for index, name in enumerate(names):
            mask[:, index] &= np.isin(table.oids, np.asarray(within[name], dtype=np.int64))
    if shown_only:
        slots = table.view()
    else:
        shown = table.view()
        hidden = np.setdiff1d(np.arange(len(table)), shown, assume_unique=True)
        slots = np.concatenate([shown, hidden])
    rows, cols = np.nonzero(mask[slots]) if len(names) else (np.empty(0, int), np.empty(0, int))
    return Matches(table.oids[slots[rows]], names, cols.astype(np.int64))


def replace_all(
    table: Table,
    query: Query,
    replacement: str,
    columns: Sequence[str] | None = None,
    *,
    shown_only: bool = True,
    within: Cells | None = None,
) -> int:
    """Replace in every matching cell as one undoable step; returns how many cells changed."""
    if not query.text:
        return 0
    matches = find(table, query, columns, shown_only=shown_only, within=within)
    updates: dict[str, tuple[np.ndarray, list[str]]] = {}
    changed = 0
    for index, name in enumerate(matches.columns):
        oids = matches.oids[matches.column_index == index]
        if not len(oids):
            continue
        old = table.column(name).displays(table.slots_of(oids))
        new = [query.replace(text, replacement) for text in old]
        keep = np.array([a != b for a, b in zip(old, new, strict=True)], dtype=np.bool_)
        if keep.any():
            updates[name] = (oids[keep], [t for t, k in zip(new, keep, strict=True) if k])
            changed += int(keep.sum())
    if updates:
        table.set_block(updates, label=f"Replace “{query.text}”")
    return changed


def replace_cell(table: Table, oid: int, column: str, query: Query, replacement: str) -> bool:
    """Replace within one cell (Replace button); returns whether the cell changed."""
    old = table.cell(column, oid)
    new = query.replace(old, replacement)
    if new == old:
        return False
    table.set_cells(column, [oid], [new], label=f"Replace “{query.text}”")
    return True
