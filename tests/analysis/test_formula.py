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

"""The Calculator's formula language: syntax, blanks, functions and what it refuses."""

from __future__ import annotations

import math

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.analysis.formula import (
    FUNCTIONS,
    MOL,
    NUM,
    TEXT,
    FormulaError,
    Value,
    parse,
    result_values,
)

COLUMNS = ["MW", "IC50 (nM)", "Series", "Structure", "HBD"]


@pytest.fixture
def data():
    smiles = ["CCO", "c1ccccc1O", None, "CC(=O)O"]
    return {
        "MW": Value(NUM, np.array([300.0, 400.0, np.nan, 250.0])),
        "IC50 (nM)": Value(NUM, np.array([10.0, 1000.0, 5.0, 0.0])),
        "Series": Value(TEXT, np.array(["A", "B", None, "A"], dtype=object)),
        "Structure": Value(
            MOL, np.array([Chem.MolFromSmiles(s) if s else None for s in smiles], dtype=object)
        ),
        "HBD": Value(NUM, np.array([1.0, 1.0, 0.0, 1.0])),
    }


def calc(text, data):
    return result_values(parse(text, COLUMNS).evaluate(data, 4), 4)


@pytest.mark.parametrize(
    ("text", "kind", "cells"),
    [
        ("MW / 2", "number", [150.0, 200.0, None, 125.0]),
        ("9 - log10([IC50 (nM)])", "number", [8.0, 6.0, 9 - math.log10(5), None]),
        ("9 - log10([ic50 (nm)])", "number", [8.0, 6.0, 9 - math.log10(5), None]),
        ("2^3 + 2**3", "number", [16.0] * 4),
        ("[IC50 (nM)] / 0", "number", [None] * 4),
        ('if(MW > 300, "big", "small")', "text", ["small", "big", "small", "small"]),
        ('Series = "A" and MW < 350', "text", ["True", "False", "False", "True"]),
        ('Series <> "A"', "text", ["False", "True", "True", "False"]),
        ('Series + "-" + text(MW)', "text", ["A-300", "B-400", "-", "A-250"]),
        ("max(MW)", "number", [400.0] * 4),
        ("min(MW, 350)", "number", [300.0, 350.0, 350.0, 250.0]),
        ("rank(MW)", "number", [2.0, 3.0, None, 1.0]),
        ("rank(MW, True)", "number", [2.0, 1.0, None, 3.0]),
        ("count(Series)", "number", [3.0] * 4),
        ("ifblank(MW, 0)", "number", [300.0, 400.0, 0.0, 250.0]),
        ("round(pi, 3)", "number", [3.142] * 4),
        ("len(Series)", "number", [1.0, 1.0, 0.0, 1.0]),
        ('contains(Series, "a")', "text", ["True", "False", "False", "True"]),
        ("Formula(Structure)", "text", ["C2H6O", "C6H6O", None, "C2H4O2"]),
        ('substructure(Structure, "[OX2H]")', "number", [1.0, 1.0, None, 1.0]),
        ('substructure(Structure, "c1ccccc1")', "number", [0.0, 1.0, None, 0.0]),
    ],
)
def test_formulas(text, kind, cells, data):
    got_kind, got = calc(text, data)
    assert got_kind == kind
    assert got == pytest.approx(cells, nan_ok=True) if kind == "number" else got == cells


def test_a_column_named_like_a_function_is_both(data):
    columns = [*COLUMNS[:-1], "MW"]
    formula = parse("MW([Structure]) - MW", columns)
    assert formula.columns == ["Structure", "MW"]
    assert formula.chemistry


def test_statistics_mark_the_formula_as_needing_every_row():
    assert parse("zscore(MW)", COLUMNS).aggregate
    assert not parse("MW * 2", COLUMNS).aggregate


def test_zscore_and_normalize(data):
    _, z = calc("zscore(MW)", data)
    finite = np.array([300.0, 400.0, 250.0])
    assert z[0] == pytest.approx((300 - finite.mean()) / finite.std(ddof=1))
    _, scaled = calc("normalize(MW)", data)
    assert scaled == pytest.approx([1 / 3, 1.0, None, 0.0], nan_ok=True)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("MW.real", "is not allowed"),
        ('__import__("os")', "Unknown function"),
        ("open('x')", "Unknown function"),
        ("[Nope] + 1", "No column named “Nope”"),
        ("Nope + 1", "Put names in brackets"),
        ("round(MW, 1, 2)", "takes 1 to 2 arguments"),
        ("(MW", "not closed"),
        ("lambda: 1", "lambda is not allowed"),
        ("", "Type a formula"),
        ("round(x=MW)", "No column named “x”"),  # = means equals, as in a spreadsheet
        ("[c for c in 'ab']", "list comprehension is not allowed"),
    ],
)
def test_refused(text, message):
    with pytest.raises(FormulaError) as err:
        parse(text, COLUMNS)
    assert message in err.value.message


def test_text_cannot_be_ordered_against_numbers(data):
    with pytest.raises(FormulaError, match="compared with"):
        calc("Series < 3", data)


def test_every_function_has_help_and_a_signature():
    for function in FUNCTIONS.values():
        assert function.help and function.signature.startswith(function.name)
