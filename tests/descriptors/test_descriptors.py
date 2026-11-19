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

import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import Crippen

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import parse_structure
from chemsuite.descriptors import medchem
from chemsuite.descriptors.catalog import CATALOG, COMMON, GROUPS, PKA_KEYS, in_group
from chemsuite.descriptors.compute import compute_chunk, compute_descriptors

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def mol(smiles: str) -> Chem.Mol:
    return Chem.MolFromSmiles(smiles)


def test_catalog_is_consistent() -> None:
    titles = [d.title for d in CATALOG.values()]
    assert len(titles) == len(set(titles)), "column titles must be unique"
    assert set(COMMON) <= set(CATALOG)
    assert sum(len(in_group(g)) for g in GROUPS) == len(CATALOG)


def test_moe_type_vsa_descriptors_match_rdkit_by_name() -> None:
    """The 57 binned VSA columns equal RDKit's descriptors of the same names, bin by bin."""
    from rdkit.Chem import Descriptors

    from chemsuite.descriptors.catalog import VSA_GROUP

    named = dict(Descriptors.descList)
    items = in_group(VSA_GROUP)
    assert len(items) == 12 + 10 + 14 + 11 + 10
    drug = mol("CN1CCC[C@H]1c1cccnc1.OC(=O)c1ccccc1OC(C)=O")
    for item in items:
        assert item.fn(drug) == pytest.approx(named[item.key](drug)), item.key
        assert item.help.startswith("Sum of ")
    assert "< -0.4" in CATALOG["SlogP_VSA1"].help and "≥ 0.6" in CATALOG["SlogP_VSA12"].help


@pytest.fixture
def pka():
    """chemsuite pKa's place taken by a stub with textbook sites."""
    from chemsuite.ionisation import profile
    from tests.ionisation.test_profile import SmilesStub

    profile.use(SmilesStub())
    yield
    profile.use(None)


def test_every_descriptor_runs_on_a_drug(pka) -> None:
    from chemsuite.conformers.etkdg import generate
    from chemsuite.core.mols import mol_from_blob, mol_to_blob

    three_d = generate(mol_from_blob(parse_structure(ASPIRIN)[0])).mol  # for the 3D shapes
    out = compute_chunk((tuple(CATALOG), [mol_to_blob(three_d)]))
    blank = []
    for key, (value,) in out.items():
        if CATALOG[key].is_text:
            ok = isinstance(value, str) and bool(value)
        elif key == "PKA_BASIC":  # aspirin has no basic site
            ok = math.isnan(value)
        else:
            ok = math.isfinite(value)
        if not ok:
            blank.append(f"{key}: {value!r} ({_why(key, three_d)})")
    assert not blank, "\n".join(blank)


@pytest.mark.parametrize("smiles", ["CC(N)C(=O)O", "CC=CC", "OC1CCC(C)CC1"])
def test_inchi_of_3d_conformers_keeps_the_stereo_the_structure_leaves_open(smiles) -> None:
    """Conformers of a structure with an undefined center or double bond each have some
    geometry; the InChI must not take the first one's for the compound's."""
    from rdkit.Chem import AllChem

    from chemsuite.core.mols import mol_to_blob

    three_d = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(three_d, 5, randomSeed=3)
    out = compute_chunk((("INCHI", "INCHIKEY"), [mol_to_blob(three_d), smiles]))
    assert out["INCHI"][0] == out["INCHI"][1] and "/t" not in out["INCHI"][0]
    assert "/b" not in out["INCHI"][0]
    assert out["INCHIKEY"][0] == out["INCHIKEY"][1]


def _why(key: str, three_d: Chem.Mol) -> str:
    """What computing *key* on its own gives or raises, to say why a descriptor came out
    blank."""
    from chemsuite.core.mols import as_drawn as _as_drawn
    from chemsuite.descriptors.compute import SHAPE_GROUP

    spec = CATALOG[key]
    try:
        return repr(spec.fn(three_d if spec.group == SHAPE_GROUP else _as_drawn(three_d)))
    except Exception as error:  # noqa: BLE001 (any error is the answer here)
        return f"{type(error).__name__}: {error}"


def test_the_pka_descriptors_run_on_a_drug(pka) -> None:
    aspirin = mol(ASPIRIN)
    for key in PKA_KEYS:
        value = CATALOG[key].fn(aspirin)
        assert math.isfinite(value) or key == "PKA_BASIC", key  # aspirin has no basic site


def test_known_values(pka) -> None:
    aspirin = mol(ASPIRIN)
    assert medchem.lipinski_violations(aspirin) == 0
    assert medchem.ro5_pass(aspirin) == "Yes"
    assert medchem.lipinski_violations(mol("C" * 40)) == 2  # MW > 500 and cLogP > 5
    assert medchem.ro5_pass(mol("C" * 40)) == "No"
    assert medchem.esol_logs(mol("c1ccccc1")) == pytest.approx(-2.0, abs=0.3)  # exp. -1.64
    assert 0.0 <= medchem.cns_mpo(aspirin) <= 6.0
    assert 0.0 <= medchem.bbb_score(aspirin) <= 6.0
    assert medchem.net_formal_charge(mol("[NH4+].[Cl-]")) == 0
    assert medchem.net_formal_charge(mol("C[N+](C)(C)C")) == 1
    # A base with pKa above 7.4 has cLogD below cLogP.
    amine = mol("CCCCN(C)C")
    assert medchem.clogd74(amine) < Crippen.MolLogP(amine)


def test_sa_score() -> None:
    assert medchem.sa_score(mol(ASPIRIN)) == pytest.approx(1.58, abs=0.01)  # Ertl's paper: 1.58
    assert medchem.sa_score(mol("c1ccccc1")) == pytest.approx(1.0)
    erythromycin = (
        "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)"
        "[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)"
        "[C@@H](O)[C@]1(C)O"
    )
    assert medchem.sa_score(mol(erythromycin)) > 5  # a complex natural product
    assert CATALOG["SA_SCORE"].title == "SAScore" and CATALOG["SA_SCORE"].group == "Drug-likeness"


def test_bbb_score_terms() -> None:
    # Outside every polynomial window only the aromatic-ring term is left.
    assert medchem.bbb_score_from_descriptors(2, 0, 100, 0, 0, 200, 1) == pytest.approx(1.0)


def test_blank_invalid_and_unparsed_cells() -> None:
    cells = [parse_structure("CCO")[0], b"", None, "c1ccccc1", "not smiles"]
    columns = compute_descriptors(cells, ["MolWt", "FORMULA"])
    assert columns["MolWt"].displays(range(5)) == ["46.07", "", "", "78.11", ""]
    assert columns["FORMULA"].displays(range(5)) == ["C2H6O", "", "", "C6H6", ""]
    assert columns["MolWt"].kind == "number" and columns["FORMULA"].kind == "text"


def test_3d_shape_of_a_structure_without_hydrogens_counts_them() -> None:
    """Hydrogens weigh in a shape (without them, PMI1 of hexane falls by half): a 3D
    structure read without them gets them placed, so it scores as the complete one."""
    from rdkit.Chem import AllChem

    from chemsuite.core.mols import mol_to_blob
    from chemsuite.descriptors.catalog import SHAPE_GROUP

    full = Chem.AddHs(mol("CCCCCC"))
    AllChem.EmbedMolecule(full, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(full)
    heavy = Chem.RemoveHs(full)
    keys = tuple(key for key, spec in CATALOG.items() if spec.group == SHAPE_GROUP)
    with_h, without_h = (compute_chunk((keys, [mol_to_blob(m)])) for m in (full, heavy))
    for key in ("PMI1", "PMI3", "NPR1", "RadiusOfGyration", "Asphericity"):
        assert without_h[key][0] == pytest.approx(with_h[key][0], rel=0.01), key


def test_integer_counts_and_no_negative_zero() -> None:
    columns = compute_descriptors(["CCO"], ["NumHDonors", "MolLogP"])
    assert columns["NumHDonors"].display(0) == "1"
    assert columns["MolLogP"].display(0) == "0.00"  # -0.0014 rounds to 0, never "-0.00"


def test_chunks_progress_and_cancel() -> None:
    cells = ["CCO"] * 25
    seen: list[tuple[int, int]] = []
    columns = compute_descriptors(
        cells, ["MolWt"], chunk_rows=10, progress=lambda d, t: seen.append((d, t))
    )
    assert len(columns["MolWt"]) == 25
    assert seen == [(10, 25), (20, 25), (25, 25)]
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        compute_descriptors(cells, ["MolWt"], chunk_rows=10, cancel=token)


def test_cancel_keeps_the_rows_already_done() -> None:
    token = CancelToken()

    def progress(done: int, total: int) -> None:
        if done >= 20:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        compute_descriptors(["CCO"] * 50, ["MolWt", "SMILES"], chunk_rows=10,
                            progress=progress, cancel=token)  # fmt: skip
    partial = stopped.value.partial
    assert len(partial["MolWt"]) == len(partial["SMILES"]) == 20
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled) as stopped:
        compute_descriptors(["CCO"] * 5, ["MolWt"], cancel=token)
    assert stopped.value.partial is None  # nothing finished


def test_unknown_key() -> None:
    with pytest.raises(KeyError, match="nope"):
        compute_descriptors(["C"], ["nope"])


def test_worker_processes_give_the_same_result() -> None:
    from chemsuite.core.jobs.parallel import process_pool

    cells = [parse_structure(s)[0] for s in (ASPIRIN, "CCN", "c1ccncc1")] * 4
    inline = compute_descriptors(cells, ["MolWt", "INCHIKEY"], chunk_rows=5)
    with process_pool(2) as pool:
        pooled = compute_descriptors(cells, ["MolWt", "INCHIKEY"], chunk_rows=5, executor=pool)
    for key in inline:
        rows = np.arange(len(cells))
        assert inline[key].displays(rows) == pooled[key].displays(rows)


def test_every_inchi_key_of_a_3d_cell_is_the_same(tmp_path) -> None:
    """The Descriptors tool, the Calculator, Copy As and InChI export give one 3D cell (its
    center undefined, its first conformer one isomer) one InChIKey."""
    from rdkit.Chem import AllChem

    from chemsuite.analysis.formula import MOL, Value, parse
    from chemsuite.app.table.structure_copy import structure_text
    from chemsuite.core.mols import mol_to_blob
    from chemsuite.core.table import Table
    from chemsuite.core.table.columns import StructureColumn
    from chemsuite.io.molecules import write_molecules

    mol = Chem.AddHs(Chem.MolFromSmiles("CC(N)C(=O)O"))
    AllChem.EmbedMultipleConfs(mol, 3, randomSeed=1)
    flat = Chem.MolToInchiKey(Chem.MolFromSmiles("CC(N)C(=O)O"))
    assert Chem.MolToInchiKey(mol) != flat  # the first conformer would say which isomer

    tool = compute_descriptors([mol_to_blob(mol)], ["INCHIKEY"])["INCHIKEY"].values()[0]
    formula = parse("InChIKey([Structure])", ["Structure"])
    calculated = formula.evaluate({"Structure": Value(MOL, np.array([mol], object))}, 1)
    copied = structure_text(mol, "inchikey")
    blobs, texts = np.array([mol_to_blob(mol)], object), np.array(["CC(N)C(=O)O"], object)
    table = Table.from_columns({"Structure": StructureColumn(blobs, texts)})
    write_molecules(table, tmp_path / "out.inchi", structure="Structure")
    exported = Chem.MolToInchiKey(
        Chem.MolFromInchi((tmp_path / "out.inchi").read_text().split("\t")[0].strip())
    )
    assert [tool, calculated.data[0], copied, exported] == [flat] * 4
