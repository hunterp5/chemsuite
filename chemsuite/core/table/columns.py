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

"""Typed, columnar cell storage.

Every column stores one value per *slot* (a row's physical position in the table) and supports
the same small set of vectorised operations, so the table, its undo history and the filters never
need to know which kind of column they are handling.

* :class:`NumberColumn`: ``float64`` values, a validity mask, and an ``int8`` display code per
  cell (see :mod:`chemsuite.core.table.numbers`), about 10 bytes per row.
* :class:`TextColumn`: dictionary-encoded; ``uint32`` codes into a shared, append-only list of
  unique strings. Repeated values (series, assay, project) cost 4 bytes per row, and per-value
  work such as text filters or sort ranks is done once per unique string, not once per row.
* :class:`StructureColumn`: one molecule per cell, as RDKit binary plus the text the user sees.
  Structures live in the table like any other cell, so undo, copy and saving need no special
  case, and a table may hold several structure columns (e.g. compound, core, R-group).
* :class:`ReactionColumn`: one reaction per cell, stored the same way (both are
  :class:`ParsedColumn`).

Blank is a first-class state in both kinds: a missing number is ``valid == False``, and a missing
text is code ``0`` (the empty string).
"""

from __future__ import annotations

import math
import re
import sys
from abc import ABC, abstractmethod
from collections import OrderedDict
from collections.abc import Callable, Sequence
from typing import ClassVar

import numpy as np

from .numbers import GENERAL, classify_number, decimals_for_value, format_number
from .style import DEFAULT_STYLE, ColumnStyle, format_fixed

NUMBER = "number"
TEXT = "text"
STRUCTURE = "structure"
REACTION = "reaction"
KINDS = (NUMBER, TEXT, STRUCTURE, REACTION)


def as_slots(slots: Sequence[int] | np.ndarray) -> np.ndarray:
    """Normalise a slot list to an ``int64`` array."""
    return np.asarray(slots, dtype=np.int64)


class Column(ABC):
    """Interface shared by all column kinds. *slots* are always ``int64`` arrays."""

    kind: ClassVar[str]
    #: Whether :meth:`sort_key` is meaningful (the table refuses to sort otherwise).
    sortable: ClassVar[bool] = True
    #: How the column is shown (precision, log10 marker, colours). Change it only through
    #: ``Table.set_style`` so the change is undoable; copies made with :meth:`styled_take`
    #: keep it.
    style: ColumnStyle = DEFAULT_STYLE

    def styled_take(self, slots: np.ndarray) -> Column:
        """:meth:`take` that keeps the column's style (snapshots, duplicates)."""
        copy = self.take(slots)
        copy.style = self.style
        return copy

    @classmethod
    @abstractmethod
    def blanks(cls, n: int) -> Column:
        """A column of *n* blank cells."""

    @abstractmethod
    def __len__(self) -> int: ...

    @abstractmethod
    def display(self, slot: int) -> str:
        """The cell exactly as the user sees it ("" when blank)."""

    @abstractmethod
    def blank_mask(self) -> np.ndarray:
        """``bool`` per row, True where the cell is blank."""

    @abstractmethod
    def values(self, slots: np.ndarray | None = None) -> np.ndarray:
        """Typed values for computation: ``float64`` with NaN for numbers; ``object`` with None
        for text."""

    @abstractmethod
    def sort_key(self) -> np.ndarray:
        """``float64`` per row whose order is the column's sort order (blanks: any value)."""

    @abstractmethod
    def match_text(self, predicate: Callable[[str], bool]) -> np.ndarray:
        """``bool`` per row: *predicate* applied to each cell's display text."""

    @abstractmethod
    def take(self, slots: np.ndarray) -> Column:
        """A new column holding the cells at *slots* (used for snapshots and copies)."""

    @abstractmethod
    def put(self, slots: np.ndarray, source: Column) -> None:
        """Overwrite the cells at *slots* with *source* (same kind, one cell per slot)."""

    @abstractmethod
    def insert(self, positions: np.ndarray, source: Column) -> None:
        """Insert *source*'s cells before *positions* (``numpy.insert`` semantics)."""

    @abstractmethod
    def delete(self, slots: np.ndarray) -> None:
        """Remove the cells at *slots*."""

    @property
    @abstractmethod
    def nbytes(self) -> int:
        """Approximate memory held by this column."""

    def displays(self, slots: np.ndarray) -> list[str]:
        return [self.display(int(slot)) for slot in slots]

    def is_all_blank(self) -> bool:
        return bool(self.blank_mask().all())


# --------------------------------------------------------------------------------------------
# Numbers
# --------------------------------------------------------------------------------------------


class NumberColumn(Column):
    kind = NUMBER
    __slots__ = ("_decimals", "_valid", "_values")

    def __init__(self, values: np.ndarray, valid: np.ndarray, decimals: np.ndarray) -> None:
        if not (len(values) == len(valid) == len(decimals)):
            raise ValueError("values, valid and decimals must have the same length")
        self._values = np.ascontiguousarray(values, dtype=np.float64)
        self._valid = np.ascontiguousarray(valid, dtype=np.bool_)
        self._decimals = np.ascontiguousarray(decimals, dtype=np.int8)

    @classmethod
    def blanks(cls, n: int) -> NumberColumn:
        return cls(np.zeros(n), np.zeros(n, dtype=np.bool_), np.full(n, GENERAL, dtype=np.int8))

    @classmethod
    def from_floats(cls, values: np.ndarray, decimals: int = GENERAL) -> NumberColumn:
        """Computed values; NaN/inf become blank."""
        values = np.asarray(values, dtype=np.float64)
        valid = np.isfinite(values)
        return cls(np.where(valid, values, 0.0), valid, np.full(len(values), decimals, np.int8))

    def __len__(self) -> int:
        return len(self._values)

    def display(self, slot: int) -> str:
        if not self._valid[slot]:
            return ""
        precision = self.style.precision
        if precision is not None:
            return format_fixed(float(self._values[slot]), precision)
        return format_number(float(self._values[slot]), int(self._decimals[slot]))

    def value(self, slot: int) -> float | None:
        """The cell's number (None when blank)."""
        return float(self._values[slot]) if self._valid[slot] else None

    def blank_mask(self) -> np.ndarray:
        return ~self._valid

    def values(self, slots: np.ndarray | None = None) -> np.ndarray:
        if slots is None:
            return np.where(self._valid, self._values, np.nan)
        return np.where(self._valid[slots], self._values[slots], np.nan)

    def sort_key(self) -> np.ndarray:
        return self._values

    def match_text(self, predicate: Callable[[str], bool]) -> np.ndarray:
        # The text depends only on (value, decimals): test each distinct pair once.
        out = np.zeros(len(self), dtype=np.bool_)
        valid_slots = np.flatnonzero(self._valid)
        if len(valid_slots):
            values, value_index = np.unique(self._values[valid_slots], return_inverse=True)
            decimals = self._decimals[valid_slots].astype(np.int64) + 128  # 0..255
            precision = self.style.precision
            if precision is not None:  # every cell shows the column's precision
                decimals[:] = 128

            def text(value: float, code: int) -> str:
                if precision is not None:
                    return format_fixed(value, precision)
                return format_number(value, code)

            pairs, pair_index = np.unique(value_index * 256 + decimals, return_inverse=True)
            hits = np.fromiter(
                (
                    predicate(text(float(values[p >> 8]), int((p & 255) - 128)))
                    for p in pairs.tolist()
                ),
                dtype=np.bool_,
                count=len(pairs),
            )
            out[valid_slots] = hits[pair_index.ravel()]
        if predicate(""):
            out[~self._valid] = True
        return out

    def take(self, slots: np.ndarray) -> NumberColumn:
        return NumberColumn(self._values[slots], self._valid[slots], self._decimals[slots])

    def put(self, slots: np.ndarray, source: Column) -> None:
        src = _require(source, NumberColumn)
        self._values[slots] = src._values
        self._valid[slots] = src._valid
        self._decimals[slots] = src._decimals

    def insert(self, positions: np.ndarray, source: Column) -> None:
        src = _require(source, NumberColumn)
        self._values = np.insert(self._values, positions, src._values)
        self._valid = np.insert(self._valid, positions, src._valid)
        self._decimals = np.insert(self._decimals, positions, src._decimals)

    def delete(self, slots: np.ndarray) -> None:
        self._values = np.delete(self._values, slots)
        self._valid = np.delete(self._valid, slots)
        self._decimals = np.delete(self._decimals, slots)

    @property
    def nbytes(self) -> int:
        return self._values.nbytes + self._valid.nbytes + self._decimals.nbytes

    def finite_range(self) -> tuple[float, float] | None:
        """``(min, max)`` over non-blank cells, or None when the column is all blank."""
        if not self._valid.any():
            return None
        vals = self._values[self._valid]
        return float(vals.min()), float(vals.max())


# --------------------------------------------------------------------------------------------
# Text
# --------------------------------------------------------------------------------------------

_NATURAL_SPLIT = re.compile(r"(\d+)")


def natural_key(text: str) -> tuple[list, str]:
    """Case-insensitive order in which digit runs compare as numbers: CPD-2 < CPD-10.

    ``re.split`` with a capturing group always alternates text, digits, text, …, so position *i*
    holds the same type in every key and the lists compare without type errors. The original
    text breaks ties ("007" vs "7", "a" vs "A") so the order is total and deterministic.
    """
    parts: list = _NATURAL_SPLIT.split(text.casefold())
    parts[1::2] = map(int, parts[1::2])
    return parts, text


class _Dictionary:
    """Append-only list of unique strings shared by a text column and its snapshots.

    Entries are never changed or removed, so codes held by an undo snapshot stay valid however
    the live column is edited. :meth:`TextColumn.compacted` builds a fresh dictionary when unused
    entries pile up.
    """

    __slots__ = ("_ranks", "index", "strings")

    def __init__(self) -> None:
        self.strings: list[str] = [""]
        self.index: dict[str, int] = {"": 0}
        self._ranks: np.ndarray | None = None

    def encode(self, text: str) -> int:
        code = self.index.get(text)
        if code is None:
            code = len(self.strings)
            self.strings.append(text)
            self.index[text] = code
        return code

    def encode_many(self, texts: Sequence[str]) -> np.ndarray:
        encode = self.encode
        return np.fromiter((encode(t) for t in texts), dtype=np.uint32, count=len(texts))

    def ranks(self) -> np.ndarray:
        """``float64`` rank per code in natural order; recomputed only after new strings."""
        if self._ranks is None or len(self._ranks) != len(self.strings):
            order = sorted(range(len(self.strings)), key=lambda c: natural_key(self.strings[c]))
            ranks = np.empty(len(order), dtype=np.float64)
            ranks[np.asarray(order, dtype=np.int64)] = np.arange(len(order), dtype=np.float64)
            self._ranks = ranks
        return self._ranks

    @property
    def nbytes(self) -> int:
        # str object + dict slot + list slot; an estimate is all memory reporting needs.
        return sum(sys.getsizeof(s) + 24 for s in self.strings)


class TextColumn(Column):
    kind = TEXT
    __slots__ = ("_codes", "_dict")

    def __init__(self, codes: np.ndarray, dictionary: _Dictionary) -> None:
        self._codes = np.ascontiguousarray(codes, dtype=np.uint32)
        self._dict = dictionary

    @classmethod
    def blanks(cls, n: int) -> TextColumn:
        return cls(np.zeros(n, dtype=np.uint32), _Dictionary())

    @classmethod
    def from_texts(cls, texts: Sequence[str]) -> TextColumn:
        dictionary = _Dictionary()
        return cls(dictionary.encode_many(texts), dictionary)

    def __len__(self) -> int:
        return len(self._codes)

    def display(self, slot: int) -> str:
        return self._dict.strings[self._codes[slot]]

    def displays(self, slots: np.ndarray) -> list[str]:
        strings = self._dict.strings
        return [strings[c] for c in self._codes[slots].tolist()]

    def blank_mask(self) -> np.ndarray:
        return self._codes == 0

    def coded(self, slots: np.ndarray) -> tuple[np.ndarray, list[str]]:
        """``(code per slot, strings)``: the dictionary encoding, code 0 being blank.

        Lets callers group or count by value without building a string per row.
        """
        return self._codes[slots].astype(np.int64), list(self._dict.strings)

    def values(self, slots: np.ndarray | None = None) -> np.ndarray:
        lookup = np.array(self._dict.strings, dtype=object)
        lookup[0] = None
        return lookup[self._codes if slots is None else self._codes[slots]]

    def sort_key(self) -> np.ndarray:
        return self._dict.ranks()[self._codes]

    def value_counts(self) -> list[tuple[str, int]]:
        """``(text, rows)`` for each text in use (one bincount over the codes)."""
        counts = np.bincount(self._codes, minlength=len(self._dict.strings))
        strings = self._dict.strings
        return [(strings[i], int(counts[i])) for i in np.flatnonzero(counts).tolist()]

    def match_text(self, predicate: Callable[[str], bool]) -> np.ndarray:
        strings = self._dict.strings
        per_code = np.fromiter((predicate(s) for s in strings), dtype=np.bool_, count=len(strings))
        return per_code[self._codes]

    def unique_values(self) -> list[str]:
        """Distinct non-blank values present in the column, in natural order."""
        used = np.unique(self._codes)
        return sorted((self._dict.strings[c] for c in used.tolist() if c), key=natural_key)

    def take(self, slots: np.ndarray) -> TextColumn:
        return TextColumn(self._codes[slots], self._dict)

    def put(self, slots: np.ndarray, source: Column) -> None:
        self._codes[slots] = self._translate(_require(source, TextColumn))

    def insert(self, positions: np.ndarray, source: Column) -> None:
        self._codes = np.insert(
            self._codes, positions, self._translate(_require(source, TextColumn))
        )

    def delete(self, slots: np.ndarray) -> None:
        self._codes = np.delete(self._codes, slots)

    def _translate(self, source: TextColumn) -> np.ndarray:
        """*source*'s codes expressed in this column's dictionary."""
        if source._dict is self._dict:
            return source._codes
        mapping = np.zeros(len(source._dict.strings), dtype=np.uint32)
        for code in np.unique(source._codes).tolist():
            mapping[code] = self._dict.encode(source._dict.strings[code])
        return mapping[source._codes]

    def compacted(self) -> TextColumn:
        """Same cells, dictionary trimmed to the values actually used (blank stays code 0)."""
        used, inverse = np.unique(self._codes, return_inverse=True)
        if not len(used) or used[0] != 0:
            used = np.concatenate([[0], used]).astype(self._codes.dtype)
            inverse = inverse + 1
        dictionary = _Dictionary()
        strings = self._dict.strings
        for code in used[1:].tolist():
            dictionary.encode(strings[code])
        return TextColumn(inverse.astype(np.uint32), dictionary)

    @property
    def unused_entries(self) -> int:
        return len(self._dict.strings) - len(np.unique(self._codes))

    @property
    def nbytes(self) -> int:
        return self._codes.nbytes + self._dict.nbytes


# --------------------------------------------------------------------------------------------
# Structures
# --------------------------------------------------------------------------------------------


class ParsedColumn(Column):
    """Cells holding text that parses to a chemical object (a molecule or a reaction).

    Each cell has a display *text* and a *blob* (RDKit binary). The blob's state tells the
    cell's status:

    * text ``""``                     blank
    * blob ``None`` (text present)    not parsed yet; files open text-first and parse later
    * blob ``b""``                    the text does not parse (kept, flagged)
    * blob non-empty                  parsed

    Parsing fills blobs in place (:meth:`fill`) because it changes no visible value.
    Subclasses say how to parse a text and load a blob.
    """

    sortable = False
    __slots__ = ("_blobs", "_texts")

    #: Loaded objects by blob, shared by the columns of one subclass.
    _cache: ClassVar[OrderedDict[bytes, object]]
    _CACHE_SIZE: ClassVar[int] = 1024

    def __init__(self, blobs: np.ndarray, texts: np.ndarray) -> None:
        if len(blobs) != len(texts):
            raise ValueError("blobs and texts must have the same length")
        self._blobs = object_array(blobs)
        self._texts = object_array(texts)

    @staticmethod
    @abstractmethod
    def parse(text: str) -> bytes:
        """The blob for one cell text (``b""`` when it does not parse)."""

    @staticmethod
    @abstractmethod
    def load(blob: bytes) -> object:
        """The object a non-empty blob holds."""

    @classmethod
    def blanks(cls, n: int):
        return cls(np.full(n, None, dtype=object), np.full(n, "", dtype=object))

    @classmethod
    def from_texts(cls, texts: Sequence[str]):
        """Unparsed cells (parse later, see :mod:`chemsuite.core.io.load`)."""
        clean = [text_of(t) for t in texts]
        return cls(np.full(len(clean), None, dtype=object), object_array(clean))

    @classmethod
    def from_parsed(cls, blobs: Sequence[bytes], texts: Sequence[str]):
        """Already-parsed cells: ``b""`` blobs mark invalid ones."""
        blob_list = [b if t else None for b, t in zip(blobs, texts, strict=True)]
        return cls(object_array(blob_list), object_array(list(texts)))

    def __len__(self) -> int:
        return len(self._texts)

    def display(self, slot: int) -> str:
        return self._texts[slot]

    def displays(self, slots: np.ndarray) -> list[str]:
        return self._texts[slots].tolist()

    def blank_mask(self) -> np.ndarray:
        return self._texts == ""

    def pending_mask(self) -> np.ndarray:
        """Cells with text that has not been parsed yet."""
        return np.equal(self._blobs, None) & ~self.blank_mask()

    def invalid_mask(self) -> np.ndarray:
        return self._blobs == b""

    def values(self, slots: np.ndarray | None = None) -> np.ndarray:
        texts = self._texts if slots is None else self._texts[slots]
        return np.where(texts == "", None, texts)

    def sort_key(self) -> np.ndarray:
        raise TypeError(f"{self.kind} columns cannot be sorted")

    def match_text(self, predicate: Callable[[str], bool]) -> np.ndarray:
        return np.fromiter((predicate(t) for t in self._texts), dtype=np.bool_, count=len(self))

    def blobs(self, slots: np.ndarray) -> list[bytes | None]:
        """Raw blobs, e.g. to send to a worker process without building objects."""
        return self._blobs[slots].tolist()

    def obj(self, slot: int):
        """The parsed object (a cached, shared object: copy it before modifying), or None.

        A cell that was not parsed yet is parsed now and its blob kept.
        """
        blob = self._blobs[slot]
        if blob is None:
            text = self._texts[slot]
            if not text:
                return None
            blob = self.parse(text)
            self._blobs[slot] = blob
        if not blob:
            return None
        cache = type(self)._cache
        value = cache.get(blob)
        if value is None:
            value = self.load(blob)
            cache[blob] = value
            if len(cache) > self._CACHE_SIZE:
                cache.popitem(last=False)
        else:
            cache.move_to_end(blob)
        return value

    def fill(self, slots: np.ndarray, blobs: Sequence[bytes]) -> None:
        """Store parse results for *slots* (``b""`` = invalid). Leaves the texts unchanged."""
        values = object_array(list(blobs))
        values[self._texts[slots] == ""] = None
        self._blobs[slots] = values

    def take(self, slots: np.ndarray):
        return type(self)(self._blobs[slots], self._texts[slots])

    def put(self, slots: np.ndarray, source: Column) -> None:
        src = _require(source, type(self))
        self._blobs[slots] = src._blobs
        self._texts[slots] = src._texts

    def insert(self, positions: np.ndarray, source: Column) -> None:
        src = _require(source, type(self))
        self._blobs = np.insert(self._blobs, positions, src._blobs)
        self._texts = np.insert(self._texts, positions, src._texts)

    def delete(self, slots: np.ndarray) -> None:
        self._blobs = np.delete(self._blobs, slots)
        self._texts = np.delete(self._texts, slots)

    @property
    def nbytes(self) -> int:
        blob_bytes = sum(len(b) + 33 for b in self._blobs.tolist() if b is not None)
        text_bytes = sum(len(t) + 49 for t in self._texts.tolist())
        return self._blobs.nbytes + self._texts.nbytes + blob_bytes + text_bytes


class StructureColumn(ParsedColumn):
    """Molecules, one per cell: the display text is SMILES as typed or read (canonical SMILES
    for molecules read from a molfile) and the blob is RDKit's molecule binary."""

    kind = STRUCTURE
    __slots__ = ()
    _cache: ClassVar[OrderedDict[bytes, object]] = OrderedDict()

    @staticmethod
    def parse(text: str) -> bytes:
        from chemsuite.core.mols import parse_structure

        return parse_structure(text)[0]

    @staticmethod
    def load(blob: bytes) -> object:
        from chemsuite.core.mols import mol_from_blob

        return mol_from_blob(blob)

    def mol(self, slot: int):
        """The RDKit molecule (cached and shared: copy it before modifying), or None."""
        return self.obj(slot)


class ReactionColumn(ParsedColumn):
    """Reactions, one per cell: the display text is reaction SMILES (as typed or read, or
    written from an RXN file) and the blob is RDKit's reaction binary."""

    kind = REACTION
    __slots__ = ()
    _cache: ClassVar[OrderedDict[bytes, object]] = OrderedDict()

    @staticmethod
    def parse(text: str) -> bytes:
        from chemsuite.core.reactions import parse_reaction

        return parse_reaction(text)[0]

    @staticmethod
    def load(blob: bytes) -> object:
        from chemsuite.core.reactions import reaction_from_blob

        return reaction_from_blob(blob)

    def reaction(self, slot: int):
        """The RDKit ``ChemicalReaction`` (cached and shared: do not modify), or None."""
        return self.obj(slot)


def object_array(values: Sequence[object] | np.ndarray) -> np.ndarray:
    """A 1-D object array, even when *values* are equal-length bytes or strings."""
    if isinstance(values, np.ndarray) and values.dtype == object:
        return values
    out = np.empty(len(values), dtype=object)
    out[:] = list(values)
    return out


def _require(source: Column, cls: type) -> object:
    if not isinstance(source, cls):
        raise TypeError(f"expected {cls.__name__}, got {type(source).__name__}")
    return source


# --------------------------------------------------------------------------------------------
# Building and converting columns
# --------------------------------------------------------------------------------------------


def is_blank_value(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value == ""
    return isinstance(value, float) and math.isnan(value)


def text_of(value: object) -> str:
    """Display text for an arbitrary Python value."""
    if is_blank_value(value):
        return ""
    if isinstance(value, str):
        return value
    parsed = decimals_for_value(value)
    if parsed is not None:
        return format_number(*parsed)
    return str(value)


def _number_cell(value: object) -> tuple[float, int] | None:
    if isinstance(value, str):
        return classify_number(value)
    return decimals_for_value(value)


def column_from_values(values: Sequence[object]) -> Column:
    """Infer the column kind: numbers when every non-blank value is a lossless number."""
    n = len(values)
    nums = np.zeros(n, dtype=np.float64)
    valid = np.zeros(n, dtype=np.bool_)
    decimals = np.full(n, GENERAL, dtype=np.int8)
    for i, value in enumerate(values):
        if is_blank_value(value):
            continue
        cell = _number_cell(value)
        if cell is None:
            return TextColumn.from_texts([text_of(v) for v in values])
        nums[i], decimals[i] = cell
        valid[i] = True
    if not valid.any():
        return TextColumn.blanks(n)  # an all-blank (or empty) column starts as text
    return NumberColumn(nums, valid, decimals)


def to_text(column: Column) -> TextColumn:
    if isinstance(column, TextColumn):
        return column
    return TextColumn.from_texts(column.displays(np.arange(len(column))))


def to_number(column: Column, *, force: bool = False) -> tuple[NumberColumn, int] | None:
    """Convert to numbers. Returns ``(column, n_blanked)``.

    Without *force* the conversion must be lossless (otherwise None). With *force*, cells that
    are not plain numbers become blank and are counted in ``n_blanked``.
    """
    if isinstance(column, NumberColumn):
        return column, 0
    n = len(column)
    nums = np.zeros(n, dtype=np.float64)
    valid = np.zeros(n, dtype=np.bool_)
    decimals = np.full(n, GENERAL, dtype=np.int8)
    blanked = 0
    for i, text in enumerate(column.displays(np.arange(n))):
        if not text:
            continue
        cell = classify_number(text.strip()) if force else classify_number(text)
        if cell is None:
            if not force:
                return None
            blanked += 1
            continue
        nums[i], decimals[i] = cell
        valid[i] = True
    return NumberColumn(nums, valid, decimals), blanked


def to_structure(column: Column) -> StructureColumn:
    """Interpret each cell's text as a structure (unparsed until first use)."""
    if isinstance(column, StructureColumn):
        return column
    return StructureColumn.from_texts(column.displays(np.arange(len(column))))


def to_reaction(column: Column) -> ReactionColumn:
    """Interpret each cell's text as a reaction (unparsed until first use)."""
    if isinstance(column, ReactionColumn):
        return column
    return ReactionColumn.from_texts(column.displays(np.arange(len(column))))


def coerce(column: Column, kind: str) -> Column | None:
    """*column* as *kind* when that is lossless, else None."""
    if kind == TEXT:
        return to_text(column)
    if kind == STRUCTURE:
        return to_structure(column)
    if kind == REACTION:
        return to_reaction(column)
    if kind == NUMBER:
        if column.is_all_blank():
            return NumberColumn.blanks(len(column))
        converted = to_number(column)
        return None if converted is None else converted[0]
    raise ValueError(f"unknown column kind {kind!r}")


def blank_column(kind: str, n: int) -> Column:
    if kind == NUMBER:
        return NumberColumn.blanks(n)
    if kind == TEXT:
        return TextColumn.blanks(n)
    if kind == STRUCTURE:
        return StructureColumn.blanks(n)
    if kind == REACTION:
        return ReactionColumn.blanks(n)
    raise ValueError(f"unknown column kind {kind!r}")


def concat_columns(parts: Sequence[Column]) -> Column:
    """Join columns end to end (used to assemble a file read in chunks).

    All-blank parts take the kind of the others. If the non-blank parts disagree, structures (or
    reactions) win over text and numbers, and text wins over numbers, so no cell is lost; a mix
    of structures and reactions becomes text.
    """
    if not parts:
        return TextColumn.blanks(0)
    kinds = {p.kind for p in parts if not p.is_all_blank()} or {parts[0].kind}
    if STRUCTURE in kinds and REACTION not in kinds:
        kind = STRUCTURE
    elif REACTION in kinds and STRUCTURE not in kinds:
        kind = REACTION
    elif kinds == {NUMBER}:
        kind = NUMBER
    else:
        kind = TEXT
    fitted = [blank_column(kind, len(p)) if p.is_all_blank() else coerce(p, kind) for p in parts]
    if kind == NUMBER:
        return NumberColumn(
            np.concatenate([p._values for p in fitted]),  # type: ignore[union-attr]
            np.concatenate([p._valid for p in fitted]),  # type: ignore[union-attr]
            np.concatenate([p._decimals for p in fitted]),  # type: ignore[union-attr]
        )
    if kind in (STRUCTURE, REACTION):
        cls = StructureColumn if kind == STRUCTURE else ReactionColumn
        return cls(
            np.concatenate([p._blobs for p in fitted]),  # type: ignore[union-attr]
            np.concatenate([p._texts for p in fitted]),  # type: ignore[union-attr]
        )
    merged = TextColumn.blanks(0)
    codes = [merged._translate(p) for p in fitted]  # type: ignore[arg-type]
    return TextColumn(np.concatenate(codes), merged._dict)
