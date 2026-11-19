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

"""ChemDraw, Marvin, Maestro, XYZ and InChI files open as tables; so does a lone molfile."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.io.formats import UnsupportedFormat, read_filter, split_suffix
from chemsuite.io.load import load_table, save_table


def smiles(table) -> list[str]:
    column = table.column("Structure")
    return [Chem.MolToSmiles(column.mol(i)) if column.mol(i) else "" for i in range(len(table))]


def named(smi: str, name: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi)
    mol.SetProp("_Name", name)
    return mol


def test_suffixes() -> None:
    assert split_suffix("a.cdxml") == ("cdxml", False)
    assert split_suffix("a.MRV") == ("mrv", False)
    assert split_suffix("a.maegz") == ("mae", True)
    assert split_suffix("a.xyz.gz") == ("xyz", True)
    with pytest.raises(UnsupportedFormat):
        split_suffix("a.cdx.gz")


def test_open_filter_lists_every_type() -> None:
    text = read_filter()
    for pattern in ("*.mol", "*.cdxml", "*.cdx", "*.mrv", "*.mae", "*.maegz", "*.xyz", "*.inchi"):
        assert pattern in text.split(";;")[0]
    assert "Reactions" not in read_filter(reactions=False)


def test_molfile(tmp_path) -> None:
    path = tmp_path / "phenol.mol"
    path.write_text(Chem.MolToMolBlock(named("Oc1ccccc1", "phenol")))
    table = load_table(path)
    assert smiles(table) == ["Oc1ccccc1"]
    assert table.cell("Title", table.oids[0]) == "phenol"


def test_chemdraw(tmp_path) -> None:
    path = tmp_path / "scheme.cdxml"
    path.write_text(Chem.MolToCDXMLBlock(Chem.MolFromSmiles("OC(=O)c1ccccc1")))
    assert smiles(load_table(path)) == ["O=C(O)c1ccccc1"]


def test_marvin(tmp_path) -> None:
    path = tmp_path / "m.mrv"
    path.write_text(Chem.MolToMrvBlock(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")))
    assert smiles(load_table(path)) == ["CC(=O)Nc1ccc(O)cc1"]


def test_maestro_properties_become_columns(tmp_path) -> None:
    path = tmp_path / "ligands.mae"
    writer = Chem.MaeWriter(str(path))
    for smi, name, value in (("CCN", "ethylamine", 6.5), ("CCO", "ethanol", 5.0)):
        mol = named(smi, name)
        mol.SetDoubleProp("r_user_pIC50", value)
        writer.write(mol)
    writer.close()
    zipped = tmp_path / "ligands.maegz"
    zipped.write_bytes(gzip.compress(path.read_bytes()))
    for source in (path, zipped):
        table = load_table(source)
        assert smiles(table) == ["CCN", "CCO"]
        assert [table.cell("Title", o) for o in table.oids] == ["ethylamine", "ethanol"]
        assert table.column("r_user_pIC50").kind == "number"
        assert [float(table.cell("r_user_pIC50", o)) for o in table.oids] == [6.5, 5.0]


def test_xyz_frames_get_bonds(tmp_path) -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    frame = Chem.MolToXYZBlock(mol)
    path = tmp_path / "traj.xyz"
    path.write_text(frame + frame)
    table = load_table(path)
    assert [
        Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(table.column("Structure").mol(i))))
        for i in range(2)
    ] == ["CCO", "CCO"]


def test_inchi_keeps_names_and_bad_lines(tmp_path) -> None:
    aspirin = Chem.MolToInchi(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    path = tmp_path / "list.inchi"
    lines = [
        f"{aspirin} aspirin",
        f"{aspirin}\tacetyl salicylic acid",
        f"ASA-1 {aspirin}",
        "AuxInfo=1/1/N:1",
        "InChI=1S/broken",
        "no identifier here",
    ]
    path.write_text("\n".join(lines) + "\n")
    table = load_table(path)
    assert len(table) == 5
    assert smiles(table)[:3] == ["CC(=O)Oc1ccccc1C(=O)O"] * 3
    titles = [table.cell("Title", o) for o in table.oids]
    assert titles == ["aspirin", "acetyl salicylic acid", "ASA-1", "", ""]
    structures = table.column("Structure")
    assert structures.mol(3) is None and structures.mol(4) is None  # kept, flagged invalid
    assert structures.display(4) == "no identifier here"


@pytest.mark.parametrize("suffix", [".cdxml", ".cdx", ".mrv", ".mae", ".xyz"])
def test_a_malformed_file_opens_as_an_invalid_row(tmp_path, suffix) -> None:
    path = tmp_path / f"broken{suffix}"
    path.write_bytes(b"garbage <not a chemistry file")
    table = load_table(path)
    assert len(table) == 1
    assert table.column("Structure").mol(0) is None
    assert table.column("Structure").display(0)


# -- writing ------------------------------------------------------------------------------------

ROWS = (("CCO", "ethanol"), ("Oc1ccccc1", "phenol"))


def sd_table(tmp_path, rows=ROWS, three_d: bool = False):
    """A table read from an SD file: Structure, Title, a number and a text column."""
    from tests.core.conftest import molblock, sdf_record

    path = tmp_path / "in.sdf"
    records = [
        sdf_record(molblock(smi, name, three_d), {"pIC50": f"{6 + i}.5", "Series": "A"})
        for i, (smi, name) in enumerate(rows)
    ]
    path.write_text("".join(records))
    return load_table(path)


def heavy_smiles(table) -> list[str]:
    column = table.column("Structure")
    return [Chem.MolToSmiles(Chem.RemoveHs(column.mol(i))) for i in range(len(table))]


def test_maestro_round_trip_keeps_titles_and_data(tmp_path) -> None:
    out = tmp_path / "out.mae"
    save_table(sd_table(tmp_path), out)
    again = load_table(out)
    assert heavy_smiles(again) == ["CCO", "Oc1ccccc1"]
    assert [again.cell("Title", o) for o in again.oids] == ["ethanol", "phenol"]
    assert again.column("pIC50").kind == "number"
    assert [float(again.cell("pIC50", o)) for o in again.oids] == [6.5, 7.5]
    assert [again.cell("Series", o) for o in again.oids] == ["A", "A"]


def test_maestro_3d_structures_keep_their_hydrogens(tmp_path) -> None:
    pose = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccccc1"))
    AllChem.EmbedMolecule(pose, randomSeed=7)
    path = tmp_path / "pose.mae"
    writer = Chem.MaeWriter(str(path))
    writer.write(pose)
    writer.close()
    structures = load_table(path).column("Structure")
    assert structures.mol(0).GetNumAtoms() == pose.GetNumAtoms()
    assert structures.display(0) == "CC(=O)Nc1ccccc1"


def test_inchi_round_trip(tmp_path) -> None:
    out = tmp_path / "out.inchi"
    save_table(sd_table(tmp_path), out)
    assert out.read_text().splitlines()[0].endswith("\tethanol")
    again = load_table(out)
    assert heavy_smiles(again) == ["CCO", "Oc1ccccc1"]
    assert [again.cell("Title", o) for o in again.oids] == ["ethanol", "phenol"]


def test_xyz_needs_3d_coordinates(tmp_path) -> None:
    out = tmp_path / "out.xyz"
    with pytest.raises(ValueError, match="3D"):
        save_table(sd_table(tmp_path), out)
    assert not out.exists()
    save_table(sd_table(tmp_path, three_d=True), out)
    again = load_table(out)
    assert heavy_smiles(again) == ["CCO", "Oc1ccccc1"]
    assert [again.cell("Title", o) for o in again.oids] == ["ethanol", "phenol"]


@pytest.mark.parametrize("suffix", [".cdxml", ".mrv"])
def test_one_molecule_formats(tmp_path, suffix) -> None:
    table = sd_table(tmp_path)
    out = tmp_path / f"out{suffix}"
    with pytest.raises(ValueError, match="one molecule"):
        save_table(table, out)
    save_table(table, out, oids=table.oids[1:])
    assert heavy_smiles(load_table(out)) == ["Oc1ccccc1"]


def test_export_filter_names_every_written_format() -> None:
    from chemsuite.io.molecules import WRITE_FORMATS

    source = Path(__file__).parents[2] / "chemsuite" / "app" / "document" / "main_window.py"
    text = source.read_text(encoding="utf-8")
    export = text[text.index("EXPORT_FILTER") : text.index("_windows")]
    for fmt in WRITE_FORMATS:
        assert f"*.{fmt})" in export
