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
"""Chemical components: from the file's own tables, from the CCD (a stub opener, no network),
and the bond orders they give a ligand."""

from __future__ import annotations

import io as bytes_io
import json
import urllib.error

import gemmi
import numpy as np
import pytest
from rdkit import Chem

from chemsuite.core.settings import POLICY_ENV
from chemsuite.protein import components, ligands, store
from chemsuite.protein.components import CcdCache

from .conftest import DATA
from .test_ligands import only_ligand, without_file_components

ACRYLAMIDE = "C=CC(=O)N1CCC[C@@H](n2nc(-c3ccc(OCc4nccn4C)c(Cl)c3)c3c(N)ncnc32)C1"


class Server:
    """Serves 0G6 as the CCD does; any other component is a 404."""

    def __init__(self) -> None:
        self.asked: list[str] = []

    def __call__(self, request, timeout=None):
        url = request.full_url
        self.asked.append(url)
        if url.endswith("/0G6.cif"):
            return bytes_io.BytesIO(gemmi_text("0g6_ccd.cif.gz").encode())
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)


def gemmi_text(name: str) -> str:
    import gzip

    return gzip.decompress((DATA / name).read_bytes()).decode()


@pytest.fixture(autouse=True)
def no_policy(tmp_path, monkeypatch):
    monkeypatch.setenv(POLICY_ENV, str(tmp_path / "no-policy.json"))


def test_the_files_own_tables_define_its_ligands(egfr) -> None:
    found = components.from_meta(egfr.meta[components.META_KEY])
    assert list(found) == ["80U"]  # ligands only: amino acids have their bundled templates
    eighty = found["80U"]
    assert len(eighty.atoms) == 35 and ("CAO", "CAP", 2, False) in eighty.bonds
    assert Chem.MolToSmiles(eighty.to_mol()) == Chem.CanonSmiles(ACRYLAMIDE.replace("@", ""))


def test_charges_follow_from_the_hydrogens_when_the_table_has_none() -> None:
    block = gemmi.cif.read(str(DATA / "5gnk_pocket.cif.gz")).sole_block()
    found = components.from_block(block, {"LYS", "ARG", "ASP"})
    assert ("NZ", "N", 1) in found["LYS"].atoms
    assert Chem.MolToSmiles(found["ARG"].to_mol()) == "NC(=[NH2+])NCCCC(N)C(=O)O"
    assert all(charge == 0 for _, _, charge in found["ASP"].atoms)  # drawn neutral, as the CCD


def test_the_files_bonds_and_orders_draw_the_ligand(egfr) -> None:
    residue = only_ligand(egfr)
    span = egfr.atoms_of(residue)
    name = {i: str(egfr.name[i]) for i in range(span.start, span.stop)}
    orders = {
        frozenset((name[int(a)], name[int(b)])): int(o)
        for (a, b), o in zip(egfr.bonds, egfr.order, strict=True)
        if int(a) in name and int(b) in name
    }
    assert orders[frozenset(("CAO", "CAP"))] == 2 and orders[frozenset(("OAE", "CBB"))] == 2


def test_a_covalent_ligand_takes_the_files_bond_orders_by_atom_name(egfr) -> None:
    ligand = ligands.to_mol(egfr, only_ligand(egfr))
    assert ligand.source == ligands.FILE and not ligand.guessed
    # the dictionary's form: the acrylamide before it reacted with Cys797, its stereo from the pose
    assert Chem.MolToSmiles(ligand.mol) == Chem.CanonSmiles(ACRYLAMIDE)


def test_the_components_are_saved_with_the_protein(egfr) -> None:
    back = store.unpack(store.pack(egfr))
    assert back.meta[components.META_KEY] == egfr.meta[components.META_KEY]
    assert ligands.to_mol(back, only_ligand(back)).source == ligands.FILE


def test_the_ccd_is_asked_once_and_kept(thrombin, tmp_path) -> None:
    server = Server()
    cache = CcdCache(tmp_path / "ccd", opener=server)
    bare = without_file_components(thrombin)
    ligand = ligands.to_mol(bare, only_ligand(bare), ccd=cache)
    # PPACK's chlorine left when it bonded to His57: the CCD flags it as a leaving atom
    assert ligand.source == ligands.CCD
    assert Chem.MolToSmiles(ligand.mol) == Chem.CanonSmiles(
        "C[C@@H](O)[C@H](CCCNC(N)=[NH2+])NC(=O)[C@@H]1CCCN1C(=O)[C@H](N)Cc1ccccc1"
    )
    assert (tmp_path / "ccd" / "0G6.cif").is_file()
    assert ligands.to_mol(bare, only_ligand(bare), ccd=cache).source == ligands.CCD
    assert server.asked == ["https://files.rcsb.org/ligands/download/0G6.cif"]


def test_a_component_the_ccd_lacks_is_not_asked_for_again(tmp_path) -> None:
    server = Server()
    cache = CcdCache(tmp_path, opener=server)
    assert cache.get("XYZ") is None and cache.get("xyz") is None
    assert len(server.asked) == 1 and (tmp_path / "XYZ.missing").exists()
    assert cache.get("../etc") is None and len(server.asked) == 1  # not a component name


def test_no_download_when_offline_or_turned_off(thrombin, tmp_path, monkeypatch) -> None:
    server = Server()
    assert CcdCache(tmp_path, opener=server, offline=True).get("0G6") is None
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"network_allowed": False}))
    monkeypatch.setenv(POLICY_ENV, str(policy))
    cache = CcdCache(tmp_path, opener=server)
    bare = without_file_components(thrombin)
    ligand = ligands.to_mol(bare, only_ligand(bare), ccd=cache)
    assert ligand.guessed and not server.asked
    assert not list(tmp_path.glob("*.missing"))  # blocked is not the same as missing


def test_cofactors_become_rows_only_when_there_is_nothing_else(kinase) -> None:
    from dataclasses import replace

    (atp,) = ligands.row_residues(kinase)
    assert str(kinase.res_name[atp]) == "ATP"
    kinds = kinase.res_kind.copy()
    manganese = int(np.flatnonzero(kinase.res_name == "MN")[0])
    kinds[manganese] = ligands.LIGAND
    assert ligands.row_residues(replace(kinase, res_kind=kinds)) == [manganese]
