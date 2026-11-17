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

from __future__ import annotations

from chemsuite.core.assets import Asset, PutAsset
from chemsuite.core.table import Table
from chemsuite.core.table.combine import append, merge


def _table(data: dict, structure: str | None = None) -> Table:
    table = Table()
    table.add_rows(data, undoable=False)
    if structure:
        table.convert_column(structure, "structure")
    table.history.mark_clean()
    return table


def cells(table: Table, name: str) -> list[str]:
    return [table.cell(name, int(o)) for o in table.oids]


def test_append_adds_every_row_and_new_columns_in_one_step() -> None:
    table = _table({"ID": ["A", "B"], "pIC50": ["6", "7"]})
    other = _table({"ID": ["C"], "pIC50": ["8.5"], "Series": ["X"]})
    other.edit("p", [PutAsset("p1", None, Asset("p1", "protein", "P1", b""))])
    oids = append(table, other, label="Append b.csv")
    assert len(oids) == 1 and cells(table, "ID") == ["A", "B", "C"]
    assert cells(table, "Series") == ["", "", "X"] and table.kind("pIC50") == "number"
    assert "p1" in table.assets
    assert table.history.undo_label == "Append b.csv"
    table.undo()
    assert len(table) == 2 and "Series" not in table.column_names and not table.assets


def test_merge_updates_matched_rows_and_adds_the_rest() -> None:
    table = _table({"ID": ["A", "B", "C"], "pIC50": ["6", "", "5"]})
    other = _table({"Compound": ["B", "A", "D", "A"], "pIC50": ["7", "", "9", "1"],
                    "LogD": ["2", "3", "", "4"]})  # fmt: skip
    done = merge(table, other, "ID", "Compound")
    assert cells(table, "ID") == ["A", "B", "C", "D", "A"]  # the second A is unmatched
    assert cells(table, "pIC50") == ["6", "7", "5", "9", "1"]  # a blank never overwrites
    assert cells(table, "LogD") == ["3", "2", "", "", "4"]
    assert "Compound" not in table.column_names
    assert done.updated == 2 and len(done.added) == 2 and done.new_columns == ["LogD"]
    assert done.unmatched == 2
    table.undo()
    assert cells(table, "pIC50") == ["6", "", "5"] and "LogD" not in table.column_names


def test_merge_can_keep_values_and_leave_unmatched_rows_out() -> None:
    table = _table({"ID": ["A", "B"], "pIC50": ["6", ""]})
    other = _table({"ID": ["A", "B", "Z"], "pIC50": ["8", "7", "1"]})
    done = merge(table, other, "ID", "ID", replace=False, add_unmatched=False)
    assert cells(table, "pIC50") == ["6", "7"] and len(table) == 2
    assert done.updated == 1 and done.unmatched == 1 and not len(done.added)


def test_merge_matches_structures_however_they_were_written() -> None:
    table = _table({"SMILES": ["OCC", "c1ccccc1"]}, structure="SMILES")
    other = _table({"SMILES": ["CCO", "C1=CC=CC=C1"], "MW": ["46.07", "78.11"]}, "SMILES")
    done = merge(table, other, "SMILES", "SMILES")
    assert cells(table, "MW") == ["46.07", "78.11"] and len(table) == 2
    assert done.updated == 2 and table.kind("SMILES") == "structure"


def test_text_into_a_number_column_makes_it_text_and_undoes() -> None:
    table = _table({"ID": ["A"], "IC50": ["5"]})
    other = _table({"ID": ["A"], "IC50": [">10"]})
    merge(table, other, "ID", "ID")
    assert table.kind("IC50") == "text" and cells(table, "IC50") == [">10"]
    table.undo()
    assert table.kind("IC50") == "number" and cells(table, "IC50") == ["5"]
