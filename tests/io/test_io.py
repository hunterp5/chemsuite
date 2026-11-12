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

import gzip
import io
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from chemsuite.core.jobs.parallel import process_pool
from chemsuite.core.table import NUMBER, TEXT, Table
from chemsuite.core.table.columns import STRUCTURE, StructureColumn
from chemsuite.io import sdf as sdf_io
from chemsuite.io.formats import UnsupportedFormat, detect_encoding, split_suffix
from chemsuite.io.load import TableBuilder, load_table, parse_pending_structures, save_table
from tests.core.conftest import DRUGS


def all_cells(table: Table) -> dict[str, list[str]]:
    slots = np.arange(len(table))
    return {name: table.column(name).displays(slots) for name in table.column_names}


# -- formats ------------------------------------------------------------------------------------


def test_suffix_detection() -> None:
    assert split_suffix("a.SDF") == ("sdf", False)
    assert split_suffix("a.sdf.gz") == ("sdf", True)
    assert split_suffix("x.smi") == ("smiles", False)
    with pytest.raises(UnsupportedFormat):
        split_suffix("a.docx")
    with pytest.raises(UnsupportedFormat):
        split_suffix("a.xlsx.gz")


def test_encoding_detection() -> None:
    assert detect_encoding("µM, Å".encode()) == "utf-8"
    assert detect_encoding(b"\xef\xbb\xbfID") == "utf-8-sig"
    assert detect_encoding("µM, Å".encode("cp1252")) == "cp1252"
    assert detect_encoding("ok µ".encode()[:-1]) == "utf-8"  # µ cut at the sample boundary


# -- SDF ----------------------------------------------------------------------------------------


def test_sdf_reads_fields_titles_and_keeps_broken_records(drug_sdf: Path) -> None:
    table = load_table(drug_sdf)
    assert len(table) == 7
    assert table.column_names[:2] == ["Structure", "Title"]
    assert {"ID", "pIC50", "Name", "Note"} <= set(table.column_names)
    assert table.kind("Structure") == STRUCTURE and table.kind("pIC50") == NUMBER
    cells = all_cells(table)
    assert cells["Title"][:2] == ["", "caffeine"]
    assert cells["Note"][1] == "line one\nline two"
    assert cells["Note"][0] == ""
    assert cells["pIC50"][:2] == ["5.00", "5.25"]
    structures = table.column("Structure")
    assert structures.invalid_mask().tolist() == [False] * 6 + [True]
    assert cells["ID"][6] == "CPD-BAD"  # the broken record's data is kept
    assert "Xx" in cells["Structure"][6]  # ...and so is its molfile
    assert cells["Structure"][2] == "CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1"


def test_sdf_gzip_and_windows_line_endings(drug_sdf: Path, tmp_path: Path) -> None:
    raw = drug_sdf.read_bytes().replace(b"\n", b"\r\n")
    gz = tmp_path / "drugs.sdf.gz"
    gz.write_bytes(gzip.compress(raw))
    assert all_cells(load_table(gz)) == all_cells(load_table(drug_sdf))


def test_sdf_parallel_equals_serial(drug_sdf: Path, monkeypatch) -> None:
    monkeypatch.setattr(sdf_io, "BLOCK_BYTES", 2000)  # force many blocks
    serial = load_table(drug_sdf)
    with process_pool(2) as pool:
        parallel = load_table(drug_sdf, executor=pool)
    assert all_cells(parallel) == all_cells(serial)


@given(
    st.lists(st.integers(1, 5), min_size=1, max_size=12),
    st.integers(1, 400),
)
def test_blocks_always_end_on_record_boundaries(sizes, block_bytes) -> None:
    records = [f"title {i}\n" + "x\n" * n + "M  END\n$$$$\n" for i, n in enumerate(sizes)]
    data = "".join(records).encode()
    blocks = list(sdf_io.iter_blocks(io.BytesIO(data), block_bytes))
    assert b"".join(blocks) == data
    assert all(block.endswith(b"$$$$\n") for block in blocks)


def test_sdf_round_trip(drug_sdf: Path, tmp_path: Path) -> None:
    table = load_table(drug_sdf)
    out = tmp_path / "out.sdf"
    save_table(table, out)
    again = load_table(out)
    assert all_cells(again) == all_cells(table)


def test_sdf_writes_2d_coordinates_for_smiles_input(tmp_path: Path) -> None:
    table = Table()
    table.add_rows({"Structure": StructureColumn.from_texts(["CCO", ""]), "ID": ["a", "b"]})
    out = tmp_path / "x.sdf"
    save_table(table, out)
    again = load_table(out)
    assert all_cells(again)["Structure"] == ["CCO", ""]
    assert again.column("Structure").mol(0).GetNumConformers() == 1


# -- delimited ----------------------------------------------------------------------------------


@pytest.mark.parametrize("delimiter", [",", "\t", ";", "|"])
def test_delimiter_is_sniffed(tmp_path: Path, delimiter: str) -> None:
    path = tmp_path / "data.txt"
    path.write_text(
        delimiter.join(["SMILES", "ID", "pIC50"])
        + "\n"
        + delimiter.join(["CCO", "A-1", "7.50"])
        + "\n",
        encoding="utf-8",
    )
    table = load_table(path)
    assert table.column_names == ["SMILES", "ID", "pIC50"]
    assert table.kind("SMILES") == STRUCTURE
    assert all_cells(table)["pIC50"] == ["7.50"]


def test_csv_quirks(tmp_path: Path) -> None:
    path = tmp_path / "quirks.csv"
    path.write_bytes(
        "ID,,ID,Value,Unit\r\n"
        '"A,1",x,dup,1.50,µM\r\n'
        '"multi\nline",y\r\n'
        "\r\n"
        "C,z,,007,nM,extra\r\n".encode("cp1252")
    )
    table = load_table(path)
    assert table.column_names == ["ID", "Column 2", "ID (2)", "Value", "Unit", "Column 6"]
    cells = all_cells(table)
    assert cells["ID"] == ["A,1", "multi\nline", "C"]
    assert cells["Unit"] == ["µM", "", "nM"]
    assert cells["Value"] == ["1.50", "", "007"]
    assert table.kind("Value") == TEXT  # "007" must not become 7
    assert cells["Column 6"] == ["", "", "extra"]
    assert "SMILES" not in table  # no structure column guessed


def test_csv_round_trip_is_lossless(tmp_path: Path) -> None:
    table = Table()
    table.add_rows(
        {
            "Structure": StructureColumn.from_texts(["CCO", "C1CC", ""]),
            "ID": ["A,1", 'say "hi"', "multi\nline"],
            "pIC50": ["7.50", "", "-3"],
        }
    )
    out = tmp_path / "out.csv"
    save_table(table, out)
    again = load_table(out, structure_column="Structure")
    assert all_cells(again) == all_cells(table)
    assert again.kind("pIC50") == NUMBER


def test_explicit_structure_column_and_none(tmp_path: Path) -> None:
    path = tmp_path / "x.csv"
    path.write_text("Mol,SMILES\nCCO,CCN\n", encoding="utf-8")
    assert load_table(path).kind("Mol") == STRUCTURE  # auto picks the first match
    assert load_table(path, structure_column="SMILES").kind("Mol") == TEXT
    assert load_table(path, structure_column=None).kind("SMILES") == TEXT
    with pytest.raises(ValueError):
        load_table(path, structure_column="nope")


def test_smiles_file(tmp_path: Path) -> None:
    path = tmp_path / "list.smi"
    path.write_text("SMILES Name\nCCO ethanol\nc1ccccc1\n\nCCN ethyl amine\n", encoding="utf-8")
    table = load_table(path)
    assert all_cells(table) == {
        "Structure": ["CCO", "c1ccccc1", "CCN"],
        "Name": ["ethanol", "", "ethyl amine"],
    }
    out = tmp_path / "out.smi"
    save_table(table, out)
    assert all_cells(load_table(out)) == all_cells(table)


def test_xlsx_round_trip(tmp_path: Path) -> None:
    table = Table()
    table.add_rows(
        {
            "SMILES": StructureColumn.from_texts(["CCO", ""]),
            "pIC50": [7.5, None],
            "ID": ["A-1", "B-2"],
        }
    )
    out = tmp_path / "out.xlsx"
    save_table(table, out)
    again = load_table(out)
    assert again.kind("SMILES") == STRUCTURE and again.kind("pIC50") == NUMBER
    np.testing.assert_array_equal(again.read(["pIC50"])["pIC50"], [7.5, np.nan])
    assert all_cells(again)["ID"] == ["A-1", "B-2"]


def test_failed_save_keeps_the_previous_file(tmp_path: Path, monkeypatch) -> None:
    out = tmp_path / "keep.csv"
    out.write_text("old\n", encoding="utf-8")
    table = Table()
    table.add_rows({"A": ["1"]})
    from chemsuite.io import delimited

    def boom(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(delimited.csv, "writer", boom)
    with pytest.raises(OSError):
        save_table(table, out)
    assert out.read_text(encoding="utf-8") == "old\n"
    assert [p.name for p in tmp_path.iterdir()] == ["keep.csv"]  # no temp file left


# -- assembling and parsing ---------------------------------------------------------------------


def test_builder_handles_columns_that_come_and_go() -> None:
    builder = TableBuilder()
    builder.add({"A": ["1", "2"], "B": ["x", "y"]})
    builder.add({"A": ["oops"], "C": ["3"]})
    builder.add({"B": ["z"]})
    table = builder.build()
    assert table.column_names == ["A", "B", "C"]
    cells = all_cells(table)
    assert cells == {
        "A": ["1", "2", "oops", ""],
        "B": ["x", "y", "", "z"],
        "C": ["", "", "3", ""],
    }
    assert table.kind("A") == TEXT and table.kind("C") == NUMBER
    assert table.oids.tolist() == [1, 2, 3, 4]


def test_parse_pending_structures_serial_and_parallel() -> None:
    texts = list(DRUGS.values()) * 3 + ["C1CC"]
    for executor in (None, process_pool(2)):
        table = Table()
        table.add_rows({"S": StructureColumn.from_texts(texts)}, undoable=False)
        try:
            invalid = parse_pending_structures(table, executor=executor)
        finally:
            if executor is not None:
                executor.shutdown()
        column = table.column("S")
        assert invalid == 1
        assert not column.pending_mask().any()
        assert column.invalid_mask().sum() == 1


def test_parse_skips_rows_edited_meanwhile(monkeypatch) -> None:
    from chemsuite.io import load

    table = Table()
    table.add_rows({"S": StructureColumn.from_texts(["CCO", "CCN"])}, undoable=False)
    real = load.parse_structures

    def edit_then_parse(texts):
        table.set_cells("S", [1], ["C1CC"])  # user edits while the job runs
        return real(texts)

    monkeypatch.setattr(load, "parse_structures", edit_then_parse)
    parse_pending_structures(table)
    column = table.column("S")
    assert table.cell("S", 1) == "C1CC"
    assert column.pending_mask().tolist() == [True, False]  # not overwritten with CCO's blob
