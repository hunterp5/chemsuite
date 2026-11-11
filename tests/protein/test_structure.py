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
"""Structures: reading, classifying, bonding, secondary structure, sequences and storage."""

from __future__ import annotations

from dataclasses import fields

import gemmi
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from chemsuite.protein import io, secondary, store
from chemsuite.protein.model import (
    ADDITIVE,
    AMINO,
    HELIX,
    ION,
    LIGAND,
    STRAND,
    WATER,
    Structure,
    empty,
)
from chemsuite.protein.sequence import chain_sequences, one_letter

from .conftest import DATA, load


def residue(s: Structure, label: str) -> int:
    return next(r for r in range(s.n_residues) if s.residue_label(r) == label)


# reading


@pytest.mark.parametrize(
    "name", ["1crn.pdb.gz", "1crn.cif.gz", "3htb.cif.gz", "1atp.cif.gz", "1ppb.cif.gz"]
)
def test_atoms_residues_and_chains_match_gemmi(name: str) -> None:
    s = load(name)
    model = gemmi.read_structure(str(DATA / name))[0]
    assert s.n_atoms == model.count_atom_sites()
    assert s.n_residues == sum(len(chain) for chain in model)
    assert list(s.chain_name) == [chain.name for chain in model]
    first = model[0][0][0]
    assert s.name[0] == first.name
    assert np.allclose(s.coords[0, 0], first.pos.tolist(), atol=1e-4)


def test_pdb_and_mmcif_give_the_same_structure() -> None:
    pdb, cif = load("1crn.pdb.gz"), load("1crn.cif.gz")
    assert np.array_equal(pdb.name, cif.name)
    assert np.allclose(pdb.coords, cif.coords)
    assert np.array_equal(pdb.bonds, cif.bonds)
    assert np.array_equal(pdb.ss, cif.ss)
    assert np.array_equal(pdb.seqres, cif.seqres)


def test_metadata_is_read(lysozyme: Structure) -> None:
    assert lysozyme.meta["id"] == "3HTB"
    assert lysozyme.meta["method"] == "X-RAY DIFFRACTION"
    assert lysozyme.meta["resolution"] == "1.81"
    assert "T4 lysozyme" in lysozyme.meta["title"]


def test_nmr_models_share_one_set_of_atoms(nmr: Structure) -> None:
    assert nmr.n_models == 2
    assert nmr.coords.shape == (2, nmr.n_atoms, 3)
    assert not np.allclose(nmr.coords[0], nmr.coords[1])
    assert (nmr.element == 1).sum() > 500  # NMR models have their hydrogens


def test_insertion_codes_are_kept(thrombin: Structure) -> None:
    codes = {thrombin.residue_label(r) for r in np.flatnonzero(thrombin.icode != "")}
    assert len(codes) == 49
    assert "L/THR 1H" in codes


def test_a_file_without_atoms_is_refused(tmp_path) -> None:
    path = tmp_path / "empty.pdb"
    path.write_text("HEADER    NOTHING\nEND\n")
    with pytest.raises(io.ProteinFileError):
        io.read(path)


def test_a_file_that_is_not_a_structure_is_refused(tmp_path) -> None:
    path = tmp_path / "junk.cif"
    path.write_text("data_x\nloop_\n_a.b\n'unterminated\n")
    with pytest.raises(io.ProteinFileError):
        io.read(path)


# kinds


def test_residues_are_sorted_into_kinds(lysozyme: Structure, kinase: Structure) -> None:
    kinds = dict(zip(*np.unique(lysozyme.res_kind, return_counts=True), strict=True))
    assert kinds == {AMINO: 163, LIGAND: 1, WATER: 220, ADDITIVE: 3}
    assert lysozyme.res_name[lysozyme.res_kind == LIGAND].tolist() == ["JZ4"]
    assert kinase.res_name[kinase.res_kind == ION].tolist() == ["MN", "MN"]
    assert kinase.res_name[kinase.res_kind == LIGAND].tolist() == ["ATP"]
    # phosphorylated residues stay in the protein
    assert kinase.res_kind[residue(kinase, "E/TPO 197")] == AMINO


# alternate locations


def test_all_alternate_locations_are_kept_and_one_is_shown(lysozyme: Structure) -> None:
    assert (lysozyme.altloc != "").sum() == 128
    shown = lysozyme.active()
    # one location per atom name in each residue
    for r in lysozyme.default_altlocs:
        span = lysozyme.atoms_of(r)
        names = lysozyme.name[span][shown[span]]
        assert len(names) == len(set(names.tolist()))
    r = next(iter(lysozyme.default_altlocs))
    other = "B" if lysozyme.default_altlocs[r] == "A" else "A"
    switched = lysozyme.active({r: other})
    span = lysozyme.atoms_of(r)
    assert set(lysozyme.altloc[span][switched[span]]) <= {"", other}


# bonds


@pytest.mark.parametrize(
    "name", ["1crn.cif.gz", "3htb.cif.gz", "1atp.cif.gz", "1ppb.cif.gz", "1d3z_2models.cif.gz"]
)
def test_every_atom_but_waters_and_ions_is_bonded_with_sensible_lengths(name: str) -> None:
    s = load(name)
    degree = np.bincount(s.bonds.ravel(), minlength=s.n_atoms)
    alone = (degree == 0) & ~np.isin(s.atom_kind, (WATER, ION))
    assert not alone.any(), [s.atom_label(i) for i in np.flatnonzero(alone)]
    length = np.linalg.norm(s.coords[0][s.bonds[:, 0]] - s.coords[0][s.bonds[:, 1]], axis=1)
    assert length.min() > 0.9 and length.max() < 2.1
    assert (s.bonds[:, 0] < s.bonds[:, 1]).all()
    assert degree.max() <= 4 + (s.element == 15).any() + (s.element == 16).any() * 2


def test_bonds_are_never_between_alternate_locations(lysozyme: Structure) -> None:
    a, b = lysozyme.altloc[lysozyme.bonds[:, 0]], lysozyme.altloc[lysozyme.bonds[:, 1]]
    assert not ((a != "") & (b != "") & (a != b)).any()


def test_disulfides_and_covalent_ligands_are_bonded(crambin: Structure, thrombin) -> None:
    labels = {(crambin.atom_label(a), crambin.atom_label(b)) for a, b in crambin.bonds}
    assert ("A/CYS 3/SG", "A/CYS 40/SG") in labels
    assert ("A/CYS 16/SG", "A/CYS 26/SG") in labels
    ppack = residue(thrombin, "H/0G6 1")
    span = thrombin.atoms_of(ppack)
    inside = (thrombin.bonds >= span.start) & (thrombin.bonds < span.stop)
    out = thrombin.bonds[inside.any(axis=1) & ~inside.all(axis=1)]
    links = {tuple(sorted(thrombin.atom_label(i) for i in pair)) for pair in out}
    assert links == {("H/0G6 1/C3", "H/HIS 57/NE2"), ("H/0G6 1/C2", "H/SER 195/OG")}


def test_peptide_bonds_join_neighbours_only(crambin: Structure) -> None:
    peptide = [(a, b) for a, b in crambin.bonds
               if crambin.name[a] == "C" and crambin.name[b] == "N"]  # fmt: skip
    assert len(peptide) == crambin.n_residues - 1
    assert all(crambin.residue_of[b] == crambin.residue_of[a] + 1 for a, b in peptide)


def test_template_bonds_have_orders(crambin: Structure) -> None:
    phe = residue(crambin, "A/PHE 13")
    span = crambin.atoms_of(phe)
    inside = ((crambin.bonds >= span.start) & (crambin.bonds < span.stop)).all(axis=1)
    orders = sorted(crambin.order[inside].tolist())
    assert orders.count(4) == 6  # the aromatic ring
    assert orders.count(2) == 1  # C=O


# secondary structure


def test_secondary_structure_comes_from_the_file(crambin: Structure) -> None:
    assert (crambin.ss == HELIX).sum() == 21
    assert (crambin.ss == STRAND).sum() == 8


@pytest.mark.parametrize(("name", "agreement"), [("1crn.cif.gz", 0.9), ("1atp.cif.gz", 0.9),
                                                  ("1ppb.cif.gz", 0.75)])  # fmt: skip
def test_assigned_secondary_structure_agrees_with_the_file(name: str, agreement: float) -> None:
    s = load(name)
    amino = s.res_kind == AMINO
    assigned = secondary.assign(s)
    assert (assigned[amino] == s.ss[amino]).mean() >= agreement
    assert (assigned[~amino] == 0).all()


def test_a_pdb_file_without_records_gets_its_secondary_structure_assigned(tmp_path) -> None:
    text = io.to_text(load("1atp.cif.gz"), io.PDB)
    lines = [line for line in text.splitlines() if not line.startswith(("HELIX", "SHEET"))]
    s = io.parse("\n".join(lines) + "\n", io.PDB)
    assert (s.ss == HELIX).sum() > 100 and (s.ss == STRAND).sum() > 30


# sequences


def test_sequences_follow_seqres_with_gaps_for_missing_residues(kinase: Structure) -> None:
    chains = chain_sequences(kinase)
    assert [kinase.chain_name[c.chain] for c in chains] == ["E", "I"]
    enzyme = chains[0]
    assert len(enzyme.names) == 350
    assert enzyme.letters.startswith("GNAAAAKKGSEQESVKEFLAKAKEDFLKKWE")
    assert (~enzyme.modelled).sum() == 350 - 336  # the unmodelled N-terminus
    modelled = enzyme.residue_at[enzyme.modelled]
    assert (np.diff(modelled) > 0).all()
    assert kinase.res_name[modelled[0]] == enzyme.names[enzyme.modelled][0]


def test_modified_residues_have_their_parents_letter() -> None:
    assert [one_letter(n) for n in ("SEP", "TPO", "MSE", "HIS", "JZ4")] == list("STMHX")


# storage


def test_a_stored_structure_comes_back_exactly(kinase: Structure) -> None:
    back = store.unpack(kinase.to_bytes())
    _assert_same(back, kinase)


def test_stored_payload_is_not_a_pickle_and_bad_bytes_are_refused() -> None:
    with pytest.raises(store.StoreError):
        store.unpack(b"not a structure")
    with pytest.raises(store.StoreError):
        store.unpack(b"\x80\x04K\x01.")  # a pickled 1


def test_an_empty_structure_round_trips() -> None:
    back = store.unpack(empty().to_bytes())
    assert back.n_atoms == 0 and back.n_residues == 0


@settings(max_examples=25, deadline=None)
@given(st.lists(st.booleans(), min_size=1, max_size=40), st.integers(0, 2**31))
def test_any_subset_round_trips_with_its_uids(pattern: list[bool], seed: int) -> None:
    s = load("3htb.cif.gz")
    keep = np.resize(np.array(pattern), s.n_atoms)
    keep[np.random.default_rng(seed).integers(0, s.n_atoms)] = True
    part = s.subset(keep)
    assert np.array_equal(part.uid, s.uid[keep])
    assert part.n_residues == len(np.unique(s.residue_of[keep]))
    _assert_same(store.unpack(part.to_bytes()), part)


def test_a_subset_keeps_the_bonds_between_kept_atoms(crambin: Structure) -> None:
    keep = crambin.res_kind == AMINO
    keep = keep[crambin.residue_of] & (crambin.residue_of < 10)
    part = crambin.subset(keep)
    old = {tuple(crambin.uid[p]) for p in crambin.bonds if keep[p].all()}
    assert {tuple(part.uid[p]) for p in part.bonds} == old


def test_written_files_read_back_the_same(kinase: Structure, tmp_path) -> None:
    for suffix in (".cif", ".pdb"):
        path = tmp_path / f"out{suffix}"
        io.write(kinase, path)
        back = io.read(path)
        assert np.array_equal(back.name, kinase.name)
        assert np.allclose(back.coords, kinase.coords, atol=1e-3)
        assert np.array_equal(back.res_kind, kinase.res_kind)
        assert np.array_equal(back.bonds, kinase.bonds)


def _assert_same(a: Structure, b: Structure) -> None:
    for f in fields(Structure):
        mine, theirs = a.__dict__[f.name], b.__dict__[f.name]
        if isinstance(theirs, np.ndarray):
            assert mine.dtype == theirs.dtype, f.name
            assert np.array_equal(mine, theirs), f.name
        else:
            assert mine == theirs, f.name
