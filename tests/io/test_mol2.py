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
"""chemsuite.io.mol2: reading through RDKit, writing with our own SYBYL types."""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from chemsuite.core.jobs.parallel import process_pool
from chemsuite.core.table import Table
from chemsuite.core.table.columns import STRUCTURE
from chemsuite.io import mol2 as mol2_io
from chemsuite.io.formats import split_suffix
from chemsuite.io.load import load_table, save_table

#: what a medicinal chemist's mol2 files hold: charged groups, aromatic N-H rings, amides,
#: sulfonamides, phosphates, nitro groups, salts
ROUND_TRIP = [
    "CCO", "c1ccncc1", "CC(N)=O", "CC(=O)[O-]", "C[NH3+]", "NC(N)=[NH2+]", "C[N+](=O)[O-]",
    "NS(=O)(=O)c1ccccc1", "c1c[nH]cn1", "c1cc[nH]c1", "c1ccsc1", "CC(=O)O", "CS(C)=O",
    "OP(=O)(O)O", "OP(=O)([O-])[O-]", "COP(=O)([O-])OC", "CC#N", "Nc1ccccc1", "C[N+](C)(C)C",
    "Cc1nc2ccccc2[nH]1", "O=c1cc[nH]cc1", "[O-]c1ccccc1", "NC(=[NH2+])c1ccccc1", "CC(=S)N",
    "O=S(=O)([O-])c1ccccc1", "[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O", "C/C=C/C", "F[B-](F)(F)F", "[Na+].[Cl-]",
]  # fmt: skip


def embedded(smiles: str) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    return Chem.RemoveHs(mol)


@pytest.mark.parametrize("smiles", ROUND_TRIP)
def test_written_records_read_back_as_the_same_molecule(smiles) -> None:
    with rdBase.BlockLogs():
        back = Chem.MolFromMol2Block(mol2_io.mol2_block(embedded(smiles), smiles))
    assert back is not None and Chem.MolToSmiles(back) == Chem.CanonSmiles(smiles)


def test_sybyl_types_of_the_usual_groups() -> None:
    def types(smiles):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        return [t for t in mol2_io.sybyl_types(mol) if t != "H"]

    assert types("CC(=O)[O-]") == ["C.3", "C.2", "O.co2", "O.co2"]
    assert types("CC(N)=O") == ["C.3", "C.2", "N.am", "O.2"]
    assert types("C[NH3+]") == ["C.3", "N.4"]
    assert types("NC(N)=[NH2+]") == ["N.pl3", "C.cat", "N.pl3", "N.pl3"]
    assert types("c1ccncc1").count("N.ar") == 1
    assert types("CS(C)=O")[1] == "S.O" and types("CS(=O)(=O)C")[1] == "S.O2"
    assert types("Nc1ccccc1")[0] == "N.pl3"


def test_records_keep_3d_coordinates_hydrogens_and_charges() -> None:
    mol = embedded("CC(=O)[O-]")
    block = mol2_io.mol2_block(mol, "acetate")
    atoms = block.split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>BOND")[0].splitlines()
    assert len(atoms) == 7  # every hydrogen written
    heavy = np.array([[float(v) for v in line.split()[2:5]] for line in atoms[:4]])
    assert np.allclose(heavy, mol.GetConformer().GetPositions(), atol=1e-3)
    charges = [float(line.split()[-1]) for line in atoms]
    assert sum(charges) == pytest.approx(-1, abs=0.01)  # Gasteiger charges sum to the charge


def test_table_round_trip_with_titles_and_gzip(tmp_path: Path) -> None:
    table = Table()
    table.add_rows(
        {"Structure": [Chem.MolToMolBlock(embedded(s)) for s in ("CCO", "c1ccncc1")],
         "Title": ["ethanol", "pyridine"]},
        undoable=False,
    )  # fmt: skip
    table.convert_column("Structure", STRUCTURE)
    out = tmp_path / "out.mol2"
    save_table(table, out)
    back = load_table(out)
    assert back.column_names[:2] == ["Structure", "Title"]
    slots = np.arange(len(back))
    assert back.column("Structure").displays(slots) == ["CCO", "c1ccncc1"]
    assert back.column("Title").displays(slots) == ["ethanol", "pyridine"]
    assert back.column("Structure").mol(0).GetConformer().Is3D()
    gz = tmp_path / "out.mol2.gz"
    gz.write_bytes(gzip.compress(out.read_bytes().replace(b"\n", b"\r\n")))
    assert split_suffix(gz) == ("mol2", True)
    assert load_table(gz).column("Title").displays(slots) == ["ethanol", "pyridine"]


def test_a_record_rdkit_cannot_read_is_kept_and_marked_invalid(tmp_path: Path) -> None:
    good = mol2_io.mol2_block(embedded("CCO"), "ok")
    broken = "@<TRIPOS>MOLECULE\nbroken\n1 0 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n"
    broken += "      1 X1   0.0 0.0 0.0 Xx 1 UNL1 0.0\n"
    path = tmp_path / "mixed.mol2"
    path.write_text("# a comment before the first record\n" + good + broken + good)
    table = load_table(path)
    assert len(table) == 3
    assert table.column("Structure").invalid_mask().tolist() == [False, True, False]
    assert table.column("Title").displays(np.arange(3)) == ["ok", "broken", "ok"]
    assert "Xx" in table.column("Structure").display(1)  # its text is kept


def test_parallel_reading_equals_serial(tmp_path: Path, monkeypatch) -> None:
    path = tmp_path / "many.mol2"
    path.write_text("".join(mol2_io.mol2_block(embedded(s), s) for s in ROUND_TRIP[:12]))
    monkeypatch.setattr(mol2_io, "BLOCK_BYTES", 900)  # force many blocks
    serial = load_table(path)
    with process_pool(2) as pool:
        parallel = load_table(path, executor=pool)
    slots = np.arange(len(serial))
    assert len(serial) == 12
    for name in serial.column_names:
        assert parallel.column(name).displays(slots) == serial.column(name).displays(slots)


@given(st.lists(st.integers(1, 6), min_size=1, max_size=10), st.integers(1, 300))
def test_blocks_always_end_on_record_boundaries(sizes, block_bytes) -> None:
    records = [f"@<TRIPOS>MOLECULE\nm{i}\n" + "x\n" * n for i, n in enumerate(sizes)]
    data = "".join(records).encode()
    blocks = list(mol2_io.iter_blocks(io.BytesIO(data), block_bytes))
    assert b"".join(blocks) == data
    assert all(block.startswith(b"@<TRIPOS>MOLECULE") for block in blocks)
