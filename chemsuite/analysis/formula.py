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

"""The Calculator's formula language: a safe expression over columns, evaluated per row.

A formula is parsed with Python's ``ast`` and only a fixed set of node types and functions is
evaluated, so nothing can import, call methods, read files or run code. Evaluation is
vectorised with numpy over all rows at once.

Syntax (spreadsheet habits work too):

* columns: ``[IC50 (nM)]`` for any name, or a bare name such as ``MW`` when the name is a
  plain identifier and not a function;
* numbers, ``"text"`` or ``'text'``, ``True`` / ``False``;
* ``+ - * / ^`` (or ``**``), ``%`` (remainder), comparisons ``= == != <> < <= > >=``,
  ``and or not``; ``+`` joins text when either side is text;
* ``if(condition, then, else)``, and Python's ``a if condition else b``.

Blank cells are NaN in numbers and empty in text; arithmetic with a blank gives a blank, and
so does dividing by zero or the log of a negative number. Functions are listed in
:data:`FUNCTIONS` (name, signature, help) and grouped for the dialog's Insert menu.

Statistics (``mean``, ``rank``, ``zscore``…) are taken over the rows the Calculator runs on.
"""

from __future__ import annotations

import ast
import keyword
import math
import re
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

NUM, TEXT, BOOL, MOL = "number", "text", "bool", "mol"


class FormulaError(ValueError):
    """The formula cannot be used; ``message`` says why in plain words."""

    def __init__(self, message: str, column: int | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.column = column  # 0-based position in the typed text, when known


@dataclass
class Value:
    kind: str
    data: object  # a numpy array of length n, or a scalar that broadcasts


@dataclass(frozen=True)
class Function:
    name: str
    group: str
    signature: str
    help: str
    fn: Callable[..., Value]
    args: tuple[int, int]  # (fewest, most); most -1 means any number
    aggregate: bool = False  # needs every row at once (mean, rank…)
    chemistry: bool = False  # needs structures


# -- helpers ------------------------------------------------------------------------------------
def _num(v: Value, n: int) -> np.ndarray:
    """*v* as float64 per row (text that is a number converts; other text is NaN)."""
    if v.kind == NUM:
        data = v.data
    elif v.kind == BOOL:
        data = np.asarray(v.data, dtype=np.float64)
    elif v.kind == TEXT:
        data = _text_to_numbers(v.data, n)
    else:
        raise FormulaError("A structure cannot be used as a number; use a function such as MW().")
    if np.ndim(data) == 0:
        return np.full(n, float(data) if data is not None else np.nan)
    return np.asarray(data, dtype=np.float64)


def _text_to_numbers(data: object, n: int) -> np.ndarray:
    def one(text: object) -> float:
        if text is None or text == "":
            return np.nan
        try:
            return float(str(text).strip().replace(",", ""))
        except ValueError:
            return np.nan

    if not isinstance(data, np.ndarray):
        return np.full(n, one(data))
    return np.fromiter((one(t) for t in data), dtype=np.float64, count=len(data))


def number_text(x: float) -> str:
    """A number as a cell would show it: integers without ``.0``, others to 10 significant
    digits; NaN is blank."""
    if x is None or not math.isfinite(x):
        return ""
    if x == int(x) and abs(x) < 1e15:
        return str(int(x))
    return f"{x:.10g}"


def _text(v: Value, n: int) -> np.ndarray:
    """*v* as text per row ("" for blanks)."""
    if v.kind == TEXT:
        data = v.data
        if not isinstance(data, np.ndarray):
            return np.full(n, "" if data is None else str(data), dtype=object)
        return np.array(["" if t is None else str(t) for t in data], dtype=object)
    if v.kind == NUM:
        return np.array([number_text(float(x)) for x in _num(v, n)], dtype=object)
    if v.kind == BOOL:
        return np.array(["True" if b else "False" for b in _bool(v, n)], dtype=object)
    from chemsuite.core.mols import smiles_of

    return np.array(["" if m is None else smiles_of(m) for m in _mols(v, n)], dtype=object)


def _bool(v: Value, n: int) -> np.ndarray:
    if v.kind == BOOL:
        data = v.data
    elif v.kind == NUM:
        nums = _num(v, n)
        data = np.where(np.isnan(nums), False, nums != 0)
    elif v.kind == TEXT:
        texts = _text(v, n)
        data = np.array([t.strip().lower() not in ("", "0", "false", "no") for t in texts])
    else:
        data = np.array([m is not None for m in _mols(v, n)])
    if np.ndim(data) == 0:
        return np.full(n, bool(data))
    return np.asarray(data, dtype=bool)


def _mols(v: Value, n: int) -> np.ndarray:
    if v.kind != MOL:
        raise FormulaError("This function needs a structure column, e.g. MW([Structure]).")
    return v.data  # type: ignore[return-value]


def _is_text(v: Value) -> bool:
    return v.kind == TEXT


# -- functions ------------------------------------------------------------------------------------
def _math(fn: Callable[[np.ndarray], np.ndarray]) -> Callable[..., Value]:
    def apply(n: int, x: Value) -> Value:
        with np.errstate(all="ignore"):
            out = fn(_num(x, n))
        return Value(NUM, np.where(np.isfinite(out), out, np.nan))

    return apply


def _round(n: int, x: Value, digits: Value | None = None) -> Value:
    places = 0 if digits is None else int(_num(digits, n)[0]) if n else 0
    return Value(NUM, np.round(_num(x, n), places))


def _pairwise(reduce: Callable[..., np.ndarray]) -> Callable[..., Value]:
    """Row by row over several arguments (blanks ignored); with one argument, over the rows."""

    def apply(n: int, *xs: Value) -> Value:
        if len(xs) == 1:
            values = _num(xs[0], n)
            finite = values[np.isfinite(values)]
            return Value(NUM, np.full(n, float(reduce(finite)) if len(finite) else np.nan))
        stack = np.vstack([_num(x, n) for x in xs])
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-blank rows give NaN
            out = reduce(stack, axis=0)
        return Value(NUM, out)

    return apply


def _clip(n: int, x: Value, lo: Value, hi: Value) -> Value:
    return Value(NUM, np.clip(_num(x, n), _num(lo, n), _num(hi, n)))


def _if(n: int, cond: Value, yes: Value, no: Value) -> Value:
    mask = _bool(cond, n)
    if _is_text(yes) or _is_text(no):
        return Value(TEXT, np.where(mask, _text(yes, n), _text(no, n)))
    if yes.kind == BOOL and no.kind == BOOL:
        return Value(BOOL, np.where(mask, _bool(yes, n), _bool(no, n)))
    return Value(NUM, np.where(mask, _num(yes, n), _num(no, n)))


def _isblank(n: int, x: Value) -> Value:
    if x.kind == NUM:
        return Value(BOOL, np.isnan(_num(x, n)))
    if x.kind == MOL:
        return Value(BOOL, np.array([m is None for m in _mols(x, n)]))
    return Value(BOOL, np.array([t == "" for t in _text(x, n)]))


def _ifblank(n: int, x: Value, fallback: Value) -> Value:
    return _if(n, _isblank(n, x), fallback, x)


def _stat(fn: Callable[[np.ndarray], float]) -> Callable[..., Value]:
    def apply(n: int, x: Value) -> Value:
        values = _num(x, n)
        finite = values[np.isfinite(values)]
        result = fn(finite) if len(finite) else np.nan
        return Value(NUM, np.full(n, float(result)))

    return apply


def _count(n: int, x: Value) -> Value:
    blank = _bool(_isblank(n, x), n)
    return Value(NUM, np.full(n, float((~blank).sum())))


def _rank(n: int, x: Value, descending: Value | None = None) -> Value:
    """1 for the smallest (or largest with descending), ties share their average rank."""
    values = _num(x, n)
    down = bool(_bool(descending, n)[0]) if descending is not None and n else False
    out = np.full(n, np.nan)
    ok = np.isfinite(values)
    if ok.any():
        keys = -values[ok] if down else values[ok]
        order = np.argsort(keys, kind="stable")
        ranks = np.empty(len(order))
        sorted_keys = keys[order]
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and sorted_keys[j + 1] == sorted_keys[i]:
                j += 1
            ranks[order[i : j + 1]] = (i + j) / 2 + 1
            i = j + 1
        out[ok] = ranks
    return Value(NUM, out)


def _zscore(n: int, x: Value) -> Value:
    values = _num(x, n)
    finite = values[np.isfinite(values)]
    if len(finite) < 2 or finite.std(ddof=1) == 0:
        return Value(NUM, np.full(n, np.nan))
    return Value(NUM, (values - finite.mean()) / finite.std(ddof=1))


def _normalize(n: int, x: Value) -> Value:
    values = _num(x, n)
    finite = values[np.isfinite(values)]
    if not len(finite) or finite.max() == finite.min():
        return Value(NUM, np.full(n, np.nan))
    return Value(NUM, (values - finite.min()) / (finite.max() - finite.min()))


def _percentile(n: int, x: Value, q: Value) -> Value:
    values = _num(x, n)
    finite = values[np.isfinite(values)]
    p = float(_num(q, n)[0]) if n else 50.0
    if not len(finite) or not 0 <= p <= 100:
        return Value(NUM, np.full(n, np.nan))
    return Value(NUM, np.full(n, float(np.percentile(finite, p))))


def _text_fn(fn: Callable[[str], object], kind: str = TEXT) -> Callable[..., Value]:
    def apply(n: int, x: Value) -> Value:
        texts = _text(x, n)
        if kind == NUM:
            return Value(NUM, np.array([float(fn(t)) for t in texts], dtype=np.float64))
        return Value(kind, np.array([fn(t) for t in texts], dtype=object))

    return apply


def _contains(n: int, x: Value, part: Value) -> Value:
    texts, parts = _text(x, n), _text(part, n)
    found = [p.lower() in t.lower() for t, p in zip(texts, parts, strict=True)]
    return Value(BOOL, np.array(found, dtype=bool))


def _concat(n: int, *xs: Value) -> Value:
    columns = [_text(x, n) for x in xs]
    return Value(TEXT, np.array(["".join(row) for row in zip(*columns, strict=True)], dtype=object))


def _number(n: int, x: Value) -> Value:
    return Value(NUM, _num(x, n))


def _text_of(n: int, x: Value) -> Value:
    return Value(TEXT, _text(x, n))


def _descriptor(key: str) -> Callable[..., Value]:
    def apply(n: int, x: Value) -> Value:
        from rdkit import rdBase

        from chemsuite.descriptors.catalog import CATALOG

        item = CATALOG[key]
        out: list[object] = []
        with rdBase.BlockLogs():
            for mol in _mols(x, n):
                try:
                    out.append(None if mol is None else item.fn(mol))
                except (RuntimeError, ValueError, ArithmeticError):
                    out.append(None)
        if item.is_text:
            return Value(TEXT, np.array(["" if v is None else str(v) for v in out], dtype=object))
        return Value(NUM, np.array([np.nan if v is None else float(v) for v in out]))

    return apply


def _substructure(n: int, x: Value, pattern: Value) -> Value:
    """How many times *pattern* (SMARTS or SMILES) occurs; 0 when absent, blank if no structure."""
    from rdkit import Chem, rdBase

    texts = _text(pattern, n)
    cache: dict[str, Chem.Mol | None] = {}
    out = np.full(n, np.nan)
    with rdBase.BlockLogs():
        for i, (mol, text) in enumerate(zip(_mols(x, n), texts, strict=True)):
            if text not in cache:
                cache[text] = Chem.MolFromSmarts(text) or Chem.MolFromSmiles(text)
            query = cache[text]
            if query is None:
                raise FormulaError(f"“{text}” is not a SMARTS or SMILES pattern.")
            if mol is not None:
                out[i] = len(mol.GetSubstructMatches(query, uniquify=True, maxMatches=10_000))
    return Value(NUM, out)


_CHEM = (
    ("MW", "MolWt", "Average molecular weight"),
    ("ExactMW", "ExactMolWt", "Monoisotopic mass"),
    ("logP", "MolLogP", "Wildman–Crippen cLogP"),
    ("TPSA", "TPSA", "Topological polar surface area"),
    ("HBD", "NumHDonors", "Hydrogen-bond donors"),
    ("HBA", "NumHAcceptors", "Hydrogen-bond acceptors"),
    ("RotB", "NumRotatableBonds", "Rotatable bonds"),
    ("Rings", "RingCount", "Ring count"),
    ("AromaticRings", "NumAromaticRings", "Aromatic rings"),
    ("HeavyAtoms", "HeavyAtomCount", "Heavy atoms"),
    ("Fsp3", "FractionCSP3", "Fraction of sp3 carbons"),
    ("Charge", "NET_CHARGE", "Net formal charge"),
    ("QED", "QED", "Quantitative estimate of drug-likeness"),
    ("SAScore", "SA_SCORE", "Synthetic accessibility, 1 (easy) to 10 (hard)"),
    ("Formula", "FORMULA", "Molecular formula (text)"),
    ("SMILES", "SMILES", "Canonical SMILES (text)"),
    ("InChIKey", "INCHIKEY", "InChIKey (text)"),
    ("Scaffold", "MURCKO", "Bemis–Murcko scaffold SMILES (text)"),
)


def _registry() -> dict[str, Function]:
    items: list[Function] = []

    def add(name, group, signature, help, fn, args, **kw):
        items.append(Function(name, group, signature, help, fn, args, **kw))

    for name, fn, help in (
        ("abs", np.abs, "Absolute value"),
        ("sqrt", np.sqrt, "Square root"),
        ("exp", np.exp, "e to the power x"),
        ("ln", np.log, "Natural logarithm"),
        ("log", np.log10, "Base-10 logarithm"),
        ("log10", np.log10, "Base-10 logarithm"),
        ("log2", np.log2, "Base-2 logarithm"),
        ("floor", np.floor, "Round down"),
        ("ceil", np.ceil, "Round up"),
        ("sin", np.sin, "Sine (radians)"),
        ("cos", np.cos, "Cosine (radians)"),
        ("tan", np.tan, "Tangent (radians)"),
    ):
        add(name, "Math", f"{name}(x)", help, _math(fn), (1, 1))
    add("round", "Math", "round(x, digits)", "Round to *digits* decimals (default 0)", _round,
        (1, 2))  # fmt: skip
    add("min", "Math", "min(a, b, …)", "Smallest of the arguments; min(x) is the column's "
        "smallest value", _pairwise(np.nanmin), (1, -1))  # fmt: skip
    add("max", "Math", "max(a, b, …)", "Largest of the arguments; max(x) is the column's "
        "largest value", _pairwise(np.nanmax), (1, -1))  # fmt: skip
    add("clip", "Math", "clip(x, low, high)", "x kept between low and high", _clip, (3, 3))
    add("if", "Logic", "if(condition, then, else)", "then where the condition holds, else "
        "otherwise", _if, (3, 3))  # fmt: skip
    add("isblank", "Logic", "isblank(x)", "True where the cell is blank", _isblank, (1, 1))
    add("ifblank", "Logic", "ifblank(x, fallback)", "x, or fallback where x is blank",
        _ifblank, (2, 2))  # fmt: skip
    for name, fn, help in (
        ("mean", np.mean, "Mean of the column"),
        ("median", np.median, "Median of the column"),
        ("std", lambda a: a.std(ddof=1) if len(a) > 1 else np.nan, "Sample standard deviation"),
        ("sum", np.sum, "Sum of the column"),
    ):
        add(name, "Statistics", f"{name}(x)", help, _stat(fn), (1, 1), aggregate=True)
    add("count", "Statistics", "count(x)", "Number of non-blank cells", _count, (1, 1),
        aggregate=True)  # fmt: skip
    add("percentile", "Statistics", "percentile(x, p)", "The p-th percentile (0–100)",
        _percentile, (2, 2), aggregate=True)  # fmt: skip
    add("rank", "Statistics", "rank(x, descending)", "1 for the smallest value (largest "
        "with descending = True); ties share the average rank", _rank, (1, 2),
        aggregate=True)  # fmt: skip
    add("zscore", "Statistics", "zscore(x)", "(x − mean) / standard deviation", _zscore, (1, 1),
        aggregate=True)  # fmt: skip
    add("normalize", "Statistics", "normalize(x)", "Scaled to 0–1 between the column's smallest "
        "and largest value", _normalize, (1, 1), aggregate=True)  # fmt: skip
    for name, fn, help, kind in (
        ("len", len, "Number of characters", NUM),
        ("upper", str.upper, "Upper case", TEXT),
        ("lower", str.lower, "Lower case", TEXT),
        ("trim", str.strip, "Without spaces at either end", TEXT),
    ):
        add(name, "Text", f"{name}(t)", help, _text_fn(fn, kind), (1, 1))
    add("contains", "Text", "contains(t, part)", "True where t contains part (any case)",
        _contains, (2, 2))  # fmt: skip
    add("concat", "Text", "concat(a, b, …)", "The values joined as text", _concat, (1, -1))
    add("number", "Text", "number(t)", "Text read as a number (blank if it is not one)",
        _number, (1, 1))  # fmt: skip
    add("text", "Text", "text(x)", "x as text", _text_of, (1, 1))
    for name, key, help in _CHEM:
        add(name, "Chemistry", f"{name}(structure)", help, _descriptor(key), (1, 1),
            chemistry=True)  # fmt: skip
    add("substructure", "Chemistry", 'substructure(structure, "SMARTS")', "How many times the "
        "pattern occurs (0 when absent)", _substructure, (2, 2), chemistry=True)  # fmt: skip
    return {f.name.lower(): f for f in items}


FUNCTIONS: dict[str, Function] = _registry()
CONSTANTS = {"pi": math.pi, "e": math.e}
GROUPS = ("Math", "Logic", "Statistics", "Text", "Chemistry")


# -- parsing ----------------------------------------------------------------------------------
@dataclass
class Formula:
    text: str
    tree: ast.Expression
    columns: list[str]  # referenced columns, in first-use order
    names: dict[str, str] = field(default_factory=dict)  # placeholder -> column
    aggregate: bool = False

    @property
    def chemistry(self) -> bool:
        return any(
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and FUNCTIONS.get(node.func.id.lower()) is not None
            and FUNCTIONS[node.func.id.lower()].chemistry
            for node in ast.walk(self.tree)
        )

    def evaluate(self, data: Mapping[str, Value], n: int) -> Value:
        """Evaluate for *n* rows; *data* maps each referenced column to its values."""
        return _Evaluator(self, data, n).visit(self.tree.body)


_BRACKET = re.compile(r"\[([^\[\]]+)\]")


def parse(text: str, columns: Sequence[str]) -> Formula:  # noqa: C901 - a hand-written parser
    """Check and compile *text*; *columns* are the table's column names."""
    if not text.strip():
        raise FormulaError("Type a formula, e.g. [MW] / [Heavy Atoms] or log10([IC50]).")
    known = set(columns)
    by_lower: dict[str, list[str]] = {}
    for name in columns:
        by_lower.setdefault(name.lower(), []).append(name)
    names: dict[str, str] = {}
    order: list[str] = []

    def placeholder(column: str) -> str:
        for key, value in names.items():
            if value == column:
                return key
        key = f"__col{len(names)}"
        names[key] = column
        return key

    source = _rewrite(text, known, by_lower, placeholder).strip()
    try:
        tree = ast.parse(source, mode="eval")
    except SyntaxError as err:
        raise FormulaError(_syntax_message(err), None) from None
    aggregate = False
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED):
            raise FormulaError(f"{_describe(node)} is not allowed in a formula.")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name):
                raise FormulaError("Only the listed functions can be called.")
            if node.keywords:
                raise FormulaError("Function arguments are given by position, e.g. round(x, 2).")
            function = FUNCTIONS.get(node.func.id.lower().removeprefix("_"))
            if function is None:
                raise FormulaError(f"Unknown function {node.func.id}().")
            fewest, most = function.args
            count = len(node.args)
            if count < fewest or (most >= 0 and count > most):
                raise FormulaError(f"{function.signature} takes {_arg_count(fewest, most)}.")
            aggregate |= function.aggregate
        elif isinstance(node, ast.Name) and node.id not in names:
            if isinstance(getattr(node, "ctx", None), ast.Load) and node.id.lower() in CONSTANTS:
                continue
            if node.id in ("True", "False"):
                continue
            if any(
                isinstance(p, ast.Call) and p.func is node for p in ast.walk(tree)
            ):  # a function name, checked above
                continue
            raise FormulaError(f"No column named “{node.id}”. Put names in brackets: [{node.id}].")
    for key in names:
        order.append(names[key])
    return Formula(text, tree, order, names, aggregate)


def _arg_count(fewest: int, most: int) -> str:
    if most < 0:
        return f"at least {fewest} argument{'s' * (fewest != 1)}"
    if fewest == most:
        return f"{fewest} argument{'s' * (fewest != 1)}"
    return f"{fewest} to {most} arguments"


def _rewrite(text: str, known: set[str], by_lower: dict[str, list[str]], placeholder) -> str:
    """Replace ``[Column]`` and bare column names with placeholders, ``if(`` with ``_if(``,
    and the spreadsheet operators ``= <> ^`` with Python's, leaving text in quotes alone."""
    out: list[str] = []
    pieces = re.split(r"(\"[^\"]*\"|'[^']*')", text)
    for k, piece in enumerate(pieces):
        if k % 2:  # a quoted string
            out.append(piece)
            continue

        def column(match: re.Match) -> str:
            name = match.group(1)
            if name in known:
                return placeholder(name)
            stripped = name.strip()
            if stripped in known:
                return placeholder(stripped)
            matches = by_lower.get(stripped.lower(), [])
            if len(matches) == 1:
                return placeholder(matches[0])
            raise FormulaError(f"No column named “{stripped}”.")

        piece = _BRACKET.sub(lambda m: f" {column(m)} ", piece)
        piece = re.sub(r"<>", "!=", piece)
        piece = re.sub(r"(?<![=!<>])=(?!=)", "==", piece)
        piece = piece.replace("^", "**")
        piece = re.sub(r"\bif\s*\(", "_if(", piece, flags=re.IGNORECASE)
        piece = re.sub(r"\bAND\b", "and", piece)
        piece = re.sub(r"\bOR\b", "or", piece)
        piece = re.sub(r"\bNOT\b", "not", piece)
        out.append(_bare_columns(piece, known, by_lower, placeholder))
    return "".join(out)


def _bare_columns(piece: str, known, by_lower, placeholder) -> str:
    """Bare identifiers that name a column (and not a function or constant) become columns."""

    def swap(match: re.Match) -> str:
        word = match.group(0)
        if word.startswith("__col") or keyword.iskeyword(word) or word in ("True", "False"):
            return word
        follows = piece[match.end() :].lstrip()
        if follows.startswith("(") and word.lower().lstrip("_") in FUNCTIONS:
            return word
        if word.lower() in CONSTANTS and word not in known:
            return word
        if word in known:
            return placeholder(word)
        matches = by_lower.get(word.lower(), [])
        if len(matches) == 1:
            return placeholder(matches[0])
        return word

    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", swap, piece)


def _syntax_message(err: SyntaxError) -> str:
    text = (err.msg or "").lower()
    if "was never closed" in text or "unexpected eof" in text:
        return "A bracket or quote is not closed."
    if "unmatched" in text:
        return "There is a closing bracket without an opening one."
    return "The formula is not complete or has a typing error."


_ALLOWED = (
    ast.Expression,
    ast.BinOp,
    ast.UnaryOp,
    ast.BoolOp,
    ast.Compare,
    ast.Call,
    ast.Name,
    ast.Load,
    ast.Constant,
    ast.IfExp,
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.Pow,
    ast.Mod,
    ast.FloorDiv,
    ast.USub,
    ast.UAdd,
    ast.Not,
    ast.And,
    ast.Or,
    ast.Eq,
    ast.NotEq,
    ast.Lt,
    ast.LtE,
    ast.Gt,
    ast.GtE,
)


def _describe(node: ast.AST) -> str:
    return {
        ast.Attribute: "A dot (.attribute)",
        ast.Subscript: "Indexing with [ ] inside an expression",
        ast.Lambda: "lambda",
        ast.List: "A list",
        ast.Dict: "A dictionary",
        ast.ListComp: "A list comprehension",
        ast.GeneratorExp: "A generator",
        ast.NamedExpr: "Assignment (:=)",
    }.get(type(node), type(node).__name__)


# -- evaluation ---------------------------------------------------------------------------------
class _Evaluator:
    def __init__(self, formula: Formula, data: Mapping[str, Value], n: int) -> None:
        self.formula = formula
        self.data = data
        self.n = n

    def visit(self, node: ast.AST) -> Value:
        method = getattr(self, f"_{type(node).__name__}")
        return method(node)

    def _Constant(self, node: ast.Constant) -> Value:
        value = node.value
        if isinstance(value, bool):
            return Value(BOOL, value)
        if isinstance(value, int | float):
            return Value(NUM, float(value))
        if isinstance(value, str):
            return Value(TEXT, value)
        raise FormulaError(f"{value!r} is not a number or text.")

    def _Name(self, node: ast.Name) -> Value:
        if node.id in self.formula.names:
            return self.data[self.formula.names[node.id]]
        if node.id in ("True", "False"):
            return Value(BOOL, node.id == "True")
        return Value(NUM, CONSTANTS[node.id.lower()])

    def _BinOp(self, node: ast.BinOp) -> Value:
        left, right = self.visit(node.left), self.visit(node.right)
        n = self.n
        if isinstance(node.op, ast.Add) and (_is_text(left) or _is_text(right)):
            a, b = _text(left, n), _text(right, n)
            return Value(TEXT, np.array([x + y for x, y in zip(a, b, strict=True)], dtype=object))
        a, b = _num(left, n), _num(right, n)
        with np.errstate(all="ignore"):
            out = {
                ast.Add: np.add,
                ast.Sub: np.subtract,
                ast.Mult: np.multiply,
                ast.Div: np.divide,
                ast.Pow: np.power,
                ast.Mod: np.mod,
                ast.FloorDiv: np.floor_divide,
            }[type(node.op)](a, b)
        return Value(NUM, np.where(np.isfinite(out), out, np.nan))

    def _UnaryOp(self, node: ast.UnaryOp) -> Value:
        operand = self.visit(node.operand)
        if isinstance(node.op, ast.Not):
            return Value(BOOL, ~_bool(operand, self.n))
        values = _num(operand, self.n)
        return Value(NUM, -values if isinstance(node.op, ast.USub) else values)

    def _BoolOp(self, node: ast.BoolOp) -> Value:
        masks = [_bool(self.visit(v), self.n) for v in node.values]
        combine = np.logical_and if isinstance(node.op, ast.And) else np.logical_or
        out = masks[0]
        for mask in masks[1:]:
            out = combine(out, mask)
        return Value(BOOL, out)

    def _Compare(self, node: ast.Compare) -> Value:
        n = self.n
        result = np.ones(n, dtype=bool)
        left = self.visit(node.left)
        for op, right_node in zip(node.ops, node.comparators, strict=True):
            right = self.visit(right_node)
            if _is_text(left) or _is_text(right):
                a, b = _text(left, n), _text(right, n)
                if isinstance(op, ast.Eq | ast.NotEq):
                    same = np.array([x == y for x, y in zip(a, b, strict=True)])
                    result &= same if isinstance(op, ast.Eq) else ~same
                    left = right
                    continue
                an, bn = _num(left, n), _num(right, n)
                if np.isnan(an).all() or np.isnan(bn).all():
                    raise FormulaError("Text can only be compared with = or !=.")
                a_, b_ = an, bn
            else:
                a_, b_ = _num(left, n), _num(right, n)
            with np.errstate(invalid="ignore"):
                part = {
                    ast.Eq: np.equal,
                    ast.NotEq: np.not_equal,
                    ast.Lt: np.less,
                    ast.LtE: np.less_equal,
                    ast.Gt: np.greater,
                    ast.GtE: np.greater_equal,
                }[type(op)](a_, b_)
            if isinstance(op, ast.NotEq):
                part = part & ~(np.isnan(a_) & np.isnan(b_))
            result &= part
            left = right
        return Value(BOOL, result)

    def _IfExp(self, node: ast.IfExp) -> Value:
        return _if(self.n, self.visit(node.test), self.visit(node.body), self.visit(node.orelse))

    def _Call(self, node: ast.Call) -> Value:
        function = FUNCTIONS[node.func.id.lower().removeprefix("_")]  # type: ignore[union-attr]
        args = [self.visit(a) for a in node.args]
        return function.fn(self.n, *args)


# -- results ------------------------------------------------------------------------------------
def result_values(value: Value, n: int) -> tuple[str, list[object]]:
    """``(column kind, cells)`` for writing: numbers as floats (None for blanks), text as
    strings, True/False as text."""
    if value.kind == NUM:
        nums = _num(value, n)
        return "number", [None if not math.isfinite(x) else float(x) for x in nums]
    if value.kind == BOOL:
        return "text", ["True" if b else "False" for b in _bool(value, n)]
    if value.kind == MOL:
        return "text", list(_text(value, n))
    return "text", [t or None for t in _text(value, n)]
