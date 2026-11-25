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

import dataclasses
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
Op = Literal["contains", "equals", "like", "<", "<=", ">", ">=", "=", "!="]

#: ``field`` then an operator then the value; the field may be quoted to hold spaces
_FIELD = re.compile(r'^(?P<field>"[^"]*"|[^"<>=!:|\s]+)(?P<op>:|<=|>=|!=|<|>|=)(?P<value>.*)$')
#: an operator with no field: the bar's own column
_BARE = re.compile(r"^(?P<op><=|>=|!=|<|>|=)(?P<value>.*)$")


@dataclass(frozen=True)
class Term:
    """One condition of a quick search: text in a cell, or a number compared.

    *column* None means the bar's own choice (all columns, or the one it is set to). ``like``
    matches the whole cell, ``*`` standing for any text.
    """

    text: str
    column: str | None = None
    op: Op = "contains"
    number: float | None = None
    negate: bool = False

    def renamed(self, old: str, new: str) -> Term:
        return dataclasses.replace(self, column=new) if self.column == old else self


@dataclass(frozen=True)
class Either:
    """Alternatives (``a | b``): met when any of *terms* is."""

    terms: tuple[Term | Either, ...]
    negate: bool = False

    def renamed(self, old: str, new: str) -> Either:
        return Either(tuple(term.renamed(old, new) for term in self.terms), self.negate)


Condition = Term | Either


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


def _split_bars(token: str) -> list[str]:
    """*token* split on ``|`` outside double quotes; an empty end marks a ``|`` at that end."""
    parts: list[str] = []
    current: list[str] = []
    quoted = False
    for char in token:
        if char == '"':
            quoted = not quoted
        if char == "|" and not quoted:
            parts.append("".join(current))
            current = []
        else:
            current.append(char)
    parts.append("".join(current))
    return parts


def _unquote(text: str) -> str:
    if len(text) >= 2 and text[0] == text[-1] == '"':
        return text[1:-1]
    return text[1:] if text.startswith('"') else text  # a quote still being typed


def _is_wildcard(raw: str) -> bool:
    """Unquoted text with a ``*``: a whole-cell pattern (quoted text is literal)."""
    return "*" in raw and not raw.startswith('"')


def parse_terms(text: str, columns: Sequence[str]) -> tuple[Condition, ...]:
    """The quick-search *text* as conditions, all of which a row must meet.

    - ``word`` – a cell contains it; several words must all be found (anywhere in the row);
    - ``"exact phrase"`` – a cell contains the phrase, spaces included;
    - ``ID:CPD-12`` – that column contains the text (``ID:"two words"`` for spaces);
    - ``pIC50>7``, ``>=``, ``<``, ``<=``, ``=``, ``!=`` – that column's number compared;
      ``Series=A`` with a non-number is the whole cell, ignoring case;
    - ``>7`` or ``=A`` with no column – the bar's column (any column when it searches all);
    - ``*`` in unquoted text – the whole cell, ``*`` standing for any text (``ID:CPD*``,
      ``*ase``); ``ID:*`` is any value, ``-ID:*`` a blank;
    - ``a | b`` – either; in ``Series=A|B`` a bare value after ``|`` keeps the column and
      operator before it (``Series!=A|B`` is neither);
    - a leading ``-`` or ``!`` excludes rows that meet the term (``-Series:B``, ``!kinase``).

    Field names match a column ignoring case; quote names with spaces (``"Mol Wt"<400``). A
    field that names no column is searched for as plain text, so ``http://`` or ``10:30`` still
    work. Terms with nothing after the operator are skipped while they are being typed.
    """
    by_name = {name.casefold(): name for name in columns}
    out: list[Condition] = []
    join = False  # a "|" waits for the next condition
    for token in tokens(text):
        negate = len(token) > 1 and (
            token.startswith("-") or (token.startswith("!") and not token.startswith("!="))
        )
        parts = _split_bars(token[1:] if negate else token)
        join = join or parts[0] == ""
        condition = _token_condition([p for p in parts if p], by_name, negate)
        if condition is not None:
            if join and out:
                out[-1] = _either(out[-1], condition)
            else:
                out.append(condition)
            join = False
        join = join or (len(parts) > 1 and parts[-1] == "")
    return tuple(out)


def _either(first: Condition, second: Condition) -> Either:
    def flat(condition: Condition) -> tuple[Condition, ...]:
        if isinstance(condition, Either) and not condition.negate:
            return condition.terms
        return (condition,)

    return Either(flat(first) + flat(second))


def _operator(part: str, by_name: Mapping[str, str]) -> re.Match[str] | None:
    """*part*'s operator (with the column it names, if any), or None for plain text."""
    bare = _BARE.match(part)
    if bare is not None:
        return bare
    found = _FIELD.match(part)
    if found is None or by_name.get(_unquote(found["field"]).casefold()) is None:
        return None
    return found


def _token_condition(
    parts: list[str], by_name: Mapping[str, str], negate: bool
) -> Condition | None:
    """One token's parts (split on ``|``) as a condition, or None while it is being typed.

    Later parts with no operator of their own reuse the first's column and operator.
    """
    if not parts:
        return None
    head = _operator(parts[0], by_name)
    prefix = "" if head is None else head[0][: len(head[0]) - len(head["value"])]
    if len(parts) > 1 and prefix.endswith("!=") and head is not None:  # Series!=A|B: neither
        prefix = prefix[:-2] + "="
        parts = [prefix + head["value"], *parts[1:]]
        negate = not negate
    field = head["field"] if head is not None and "field" in head.re.groupindex else ""
    terms: list[Term] = []
    for index, part in enumerate(parts):
        own = _operator(part, by_name) if index else None
        if index and own is None:
            part = prefix + part
        elif field and own is not None and "field" not in own.re.groupindex:
            part = field + part  # pIC50<5|>8: the column again
        term = _part_term(part, by_name)
        if term is not None:
            terms.append(term)
    if not terms:
        return None
    if len(terms) == 1:
        return dataclasses.replace(terms[0], negate=terms[0].negate != negate)
    return Either(tuple(terms), negate)


def _part_term(body: str, by_name: Mapping[str, str]) -> Term | None:
    """One condition with no ``|`` or exclusion prefix, or None while it is being typed."""
    match = _operator(body, by_name)
    if match is None:
        if _is_wildcard(body):
            return Term(body, op="like")
        plain = _unquote(body)
        return Term(plain) if plain else None
    column = (
        by_name[_unquote(match["field"]).casefold()] if "field" in match.re.groupindex else None
    )
    op, raw = match["op"], match["value"]
    value = _unquote(raw)
    if not value:
        return None
    number = _number(value)
    if op == ":":
        return Term(value, column, "like" if _is_wildcard(raw) else "contains")
    if op != ":" and number is not None:
        return Term(value, column, op, number)
    if op in ("=", "!="):
        return Term(value, column, "like" if _is_wildcard(raw) else "equals", negate=op == "!=")
    return Term(_unquote(body))  # a comparison with no number: look for the text as typed


def like_pattern(text: str) -> str:
    """A regular expression for the whole cell, ``*`` standing for any text; ``*`` alone is
    any value (not a blank)."""
    if not text.strip("*"):
        return r"(?s)."
    return r"(?s)\A" + ".*".join(re.escape(piece) for piece in text.split("*")) + r"\Z"


def terms_mask(
    table: Table,
    terms: Sequence[Condition],
    *,
    column: str | None = None,
    case_sensitive: bool = False,
) -> np.ndarray:
    """``bool`` per row (storage order): rows that meet every condition.

    *column* is where terms without a field look (None: every searchable column).
    """
    keep = np.ones(len(table), dtype=np.bool_)
    for term in terms:
        keep &= _condition_mask(table, term, column, case_sensitive)
    return keep


def _condition_mask(
    table: Table, condition: Condition, column: str | None, case_sensitive: bool
) -> np.ndarray:
    if isinstance(condition, Either):
        hit = np.zeros(len(table), dtype=np.bool_)
        for term in condition.terms:
            hit |= _condition_mask(table, term, column, case_sensitive)
    else:
        hit = _term_mask(table, condition, column, case_sensitive)
    return ~hit if condition.negate else hit


def _term_mask(table: Table, term: Term, column: str | None, case_sensitive: bool) -> np.ndarray:
    target = term.column or column
    if target is not None and target not in table:
        return np.zeros(len(table), dtype=np.bool_)
    names = searchable_columns(table, None if target is None else [target])
    if term.number is not None:
        hit = np.zeros(len(table), dtype=np.bool_)
        for name in names:
            hit |= _compare(table, name, term.op, term.number)
        return hit
    if term.op == "like":
        query = Query(like_pattern(term.text), "regex", case_sensitive)
    else:
        query = Query(term.text, "whole" if term.op == "equals" else "contains", case_sensitive)
    return match_mask(table, query, names).any(axis=1)


def _compare(table: Table, name: str, op: Op, number: float) -> np.ndarray:
    """``bool`` per row: the number in column *name* compared with *number*."""
    from .columns import NumberColumn, to_number

    values = table.column(name)
    if not isinstance(values, NumberColumn):
        values = to_number(values, force=True)[0]
    numbers = values.values()
    with np.errstate(invalid="ignore"):
        if op == "<":
            return numbers < number
        if op == "<=":
            return numbers <= number
        if op == ">":
            return numbers > number
        if op == ">=":
            return numbers >= number
        if op == "=":
            return np.isclose(numbers, number)
        return ~np.isclose(numbers, number) & ~np.isnan(numbers)


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
