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
"""A docking run saved as files for the command line."""

from __future__ import annotations

import subprocess

from rdkit import Chem

from chemsuite.macromolecule.dock import export, vina
from chemsuite.macromolecule.dock.box import Box

BOX = Box((1.0, -2.5, 3.25), (20.0, 18.0, 16.0))


def _config(text: str) -> dict[str, str]:
    """A Vina config file as its settings (comments and blank lines left out)."""
    out = {}
    for line in text.splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            key, _, value = line.partition("=")
            out[key.strip()] = value.strip()
    return out


def test_the_config_carries_the_receptor_box_and_conditions() -> None:
    options = vina.Options(scoring=vina.VINARDO, exhaustiveness=16, poses=5, energy_range=2.5,
                           seed=7)  # fmt: skip
    found = _config(export.config_text(BOX, options))
    assert found["receptor"] == export.RECEPTOR
    assert [found[k] for k in ("center_x", "center_y", "center_z")] == ["1.000", "-2.500", "3.250"]
    assert [found[k] for k in ("size_x", "size_y", "size_z")] == ["20.000", "18.000", "16.000"]
    assert found["scoring"] == "vinardo" and found["exhaustiveness"] == "16"
    assert found["num_modes"] == "5" and found["energy_range"] == "2.5" and found["seed"] == "7"
    # the ligands are the command's argument, not a setting; cpu is left to Vina
    assert "ligand" not in found and "batch" not in found and "cpu" not in found


def test_a_config_without_a_receptor_names_none() -> None:
    """Saving the box alone (no receptor file) must not leave a line pointing at nothing."""
    found = _config(export.config_text(BOX, vina.Options(), receptor=None))
    assert "receptor" not in found and found["center_x"] == "1.000"


def _ligand(name: str, smiles: str, row: int, conformer: int = 1, of: int = 1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    from rdkit.Chem import AllChem

    AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
    return export.LigandInput(mol, name, row, conformer, of)


def test_file_names_are_numbered_safe_and_per_conformer() -> None:
    plain = export.file_stem(7, _ligand("imatinib", "c1ccccc1", 11))
    assert plain == "0007_imatinib"
    messy = export.file_stem(1, _ligand("N-(2 chloro)/amide, 3", "c1ccccc1", 11))
    assert messy == "0001_N-_2_chloro_amide_3"  # no separators a file system would mind
    many = export.file_stem(2, _ligand("x", "c1ccccc1", 11, conformer=4, of=9))
    assert many == "0002_x_c04"
    blank = export.file_stem(3, _ligand("", "c1ccccc1", 11))
    assert blank == "0003"


def test_saving_everything_writes_what_vina_needs(tmp_path) -> None:
    ligands = [
        _ligand("phenol", "Oc1ccccc1", 11),
        _ligand("aniline", "Nc1ccccc1", 12, conformer=1, of=2),
        _ligand("aniline", "Nc1ccccc1", 12, conformer=2, of=2),
    ]
    saved = export.save(tmp_path, BOX, vina.Options(), "REMARK receptor\n", ligands)
    assert saved.ligands == 3 and saved.rows == 2 and not saved.failed
    assert (tmp_path / export.CONFIG).exists()
    assert (tmp_path / export.RECEPTOR).read_text(encoding="utf-8") == "REMARK receptor\n"
    written = sorted(p.name for p in (tmp_path / export.LIGANDS).iterdir())
    assert written == ["0001_phenol.pdbqt", "0002_aniline_c01.pdbqt", "0003_aniline_c02.pdbqt"]
    assert (tmp_path / export.OUT).is_dir()  # --dir wants it there
    rows = export.read_manifest(tmp_path / export.MANIFEST)
    assert [r["file"] for r in rows] == written
    assert [r["row"] for r in rows] == ["11", "12", "12"]
    assert [r["conformer"] for r in rows] == ["1", "1", "2"]
    assert [r["name"] for r in rows] == ["phenol", "aniline", "aniline"]
    assert set(saved.files) >= {export.CONFIG, export.RECEPTOR, export.MANIFEST, export.RUN_SH,
                                export.RUN_BAT}  # fmt: skip
    for name in (export.RUN_SH, export.RUN_BAT):
        assert export.command_line() in (tmp_path / name).read_text(encoding="utf-8")


def test_the_box_alone_writes_only_the_config(tmp_path) -> None:
    saved = export.save(tmp_path, BOX, vina.Options())
    assert saved.files == (export.CONFIG,) and saved.ligands == 0
    assert [p.name for p in tmp_path.iterdir()] == [export.CONFIG]


def test_a_ligand_meeko_cannot_type_is_reported_and_the_rest_written(tmp_path) -> None:
    flat = export.LigandInput(Chem.MolFromSmiles("Oc1ccccc1"), "no 3D", 5)  # no coordinates
    saved = export.save(tmp_path, BOX, vina.Options(), None, [flat, _ligand("ok", "C", 6)])
    assert saved.ligands == 1 and saved.rows == 1
    assert [name for name, _why in saved.failed] == ["no 3D"]
    assert "3D" in saved.failed[0][1]
    assert [r["name"] for r in export.read_manifest(tmp_path / export.MANIFEST)] == ["ok"]


def test_the_saved_folder_docks_at_the_command_line(pocket, tmp_path) -> None:
    """The point of the feature: the folder chemsuite writes runs as it stands, with the
    ligands passed as one --batch argument, and Vina names its output after each input."""
    ligands = [export.LigandInput(pocket.mol, "80U", 1)]
    saved = export.save(tmp_path, pocket.box, vina.Options(exhaustiveness=1, poses=1),
                        pocket.receptor.pdbqt, ligands)  # fmt: skip
    assert not saved.failed
    command = [str(vina.executable()), "--config", export.CONFIG, "--batch", export.LIGANDS,
               "--dir", export.OUT]  # fmt: skip
    done = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=900)
    assert done.returncode == 0, done.stderr or done.stdout
    out = sorted(p.name for p in (tmp_path / export.OUT).iterdir())
    assert out == ["0001_80U_out.pdbqt"]  # the manifest's name, with Vina's suffix
    poses = vina._poses((tmp_path / export.OUT / out[0]).read_text(encoding="utf-8"))
    assert poses and poses[0].affinity < 0
