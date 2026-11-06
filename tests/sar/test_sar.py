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

from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.sar.free_wilson import SKIP_RARE, free_wilson
from chemsuite.sar.hansch import HanschSettings, hansch, load_constants, table_key
from chemsuite.sar.matrix import build_matrix
from chemsuite.sar.mmp import find_pairs, fragment, readable, transforms
from chemsuite.sar.regression import least_squares
from chemsuite.sar.rgroups import (
    Decomposition,
    assemble,
    decompose,
    parse_core,
    substituent_label,
    suggest_core,
)
from chemsuite.sar.sali import sali_pairs

R1 = {"H": "", "Cl": "Cl", "C": "C", "OC": "OC"}
R2 = {"H": "", "F": "F", "N": "N"}
EFFECT1 = {"H": 0.0, "Cl": 1.0, "C": 0.5, "OC": -0.3}
EFFECT2 = {"H": 0.0, "F": 0.7, "N": -0.4}


def _series(skip=(("Cl", "F"),)):
    """Para/meta-substituted benzamides: activity = 5 + effect(R1) + effect(R2)."""
    smiles, activity = [], []
    for a, sa in R1.items():
        for b, sb in R2.items():
            if (a, b) in skip:
                continue
            ring = f"c1cc({sa})c({sb})cc1C(=O)N".replace("()", "")
            smiles.append(ring)
            activity.append(5 + EFFECT1[a] + EFFECT2[b])
    return smiles, np.array(activity)


CORE = "[*:1]c1ccc(C(N)=O)cc1[*:2]"


def test_decompose_series() -> None:
    smiles, _ = _series()
    mols = [Chem.MolFromSmiles(s) for s in smiles] + [Chem.MolFromSmiles("CCO"), None]
    d = decompose(mols, parse_core(CORE))
    assert d.positions == ("R1", "R2")
    assert d.matched.tolist() == [True] * len(smiles) + [False, False]
    assert set(d.labels["R1"][: len(smiles)]) == {"H", "*Cl", "*C", "*OC"}
    assert set(d.labels["R2"][: len(smiles)]) == {"H", "*F", "*N"}
    assert d.varying() == ("R1", "R2")
    rebuilt = assemble(d.core, {"R1": d.fragments[("R1", "*Cl")], "R2": d.fragments[("R2", "*F")]})
    assert Chem.MolToSmiles(rebuilt) == Chem.CanonSmiles("NC(=O)c1ccc(Cl)c(F)c1")


def test_labels_and_cores() -> None:
    assert substituent_label("C[*:1]") == "*C"
    assert substituent_label("[H][*:2]") == "H"
    assert substituent_label("[*:1]OC") == "*OC"
    with pytest.raises(ValueError):
        parse_core("")
    with pytest.raises(ValueError):
        parse_core("not((a core")
    core = suggest_core([Chem.MolFromSmiles(s) for s in ("c1ccccc1CO", "c1ccccc1CN")])
    assert Chem.MolFromSmarts(core).GetNumAtoms() >= 6
    with pytest.raises(ValueError):
        suggest_core([Chem.MolFromSmiles("C")])


def test_least_squares_statistics() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(30, 2))
    y = 1 + 2 * x[:, 0] - x[:, 1] + rng.normal(0, 0.1, 30)
    fit = least_squares(x, y, ["a", "b"])
    assert fit.coefficients == pytest.approx([2, -1], abs=0.1) and fit.intercept == pytest.approx(
        1, abs=0.1
    )
    assert 0.98 < fit.q2 < fit.r2 < 1 and fit.parameters == 3
    collinear = least_squares(np.column_stack([x[:, 0], 2 * x[:, 0]]), y, ["a", "2a"])
    assert collinear.estimable.tolist() == [True, False]


def test_free_wilson_recovers_effects() -> None:
    smiles, activity = _series()
    d = decompose([Chem.MolFromSmiles(s) for s in smiles], parse_core(CORE))
    model = free_wilson(d, activity, min_count=1)
    assert model.fit.r2 == pytest.approx(1.0)
    assert model.reference == {"R1": "H", "R2": "H"}
    assert model.term("R1", "*Cl").fujita_ban == pytest.approx(1.0)
    assert model.term("R2", "*N").fujita_ban == pytest.approx(-0.4)
    for p in model.positions:  # Free-Wilson contributions are centred
        terms = [t for t in model.terms if t.position == p]
        assert sum(t.count * t.free_wilson for t in terms) == pytest.approx(0, abs=1e-9)
    best = model.best_unmade(5)
    assert best[0][0] == ("*Cl", "*F") and best[0][1] == pytest.approx(6.7)
    assert len(best) == 1  # every other combination was made
    assert np.allclose(model.predicted, activity)


def test_free_wilson_rare_groups_and_errors() -> None:
    labels = {"R1": ["H", "H", "*Cl", "*Cl", "*Br"], "R2": ["H", "*F", "H", "*F", "H"]}
    d = Decomposition(("R1", "R2"), labels, np.ones(5, bool), "")
    model = free_wilson(d, np.array([5, 6, 6, 7, 9.0]), min_count=2)
    assert model.skipped[SKIP_RARE] == 1 and math.isnan(model.predicted[4])
    with pytest.raises(ValueError):
        free_wilson(d, np.array([5, 6, 6, 7, 9.0]), min_count=5)
    flat = Decomposition(("R1",), {"R1": ["H"] * 5}, np.ones(5, bool), "")
    with pytest.raises(ValueError, match="varies"):
        free_wilson(flat, np.ones(5))


def test_best_unmade_without_enumerating() -> None:
    labels = {
        f"R{p}": [f"*G{(i * (p + 3)) % 40}" for i in range(400)] for p in range(1, 5)
    }  # 40⁴ = 2.56M combinations
    d = Decomposition(tuple(labels), labels, np.ones(400, bool), "")
    rng = np.random.default_rng(1)
    model = free_wilson(d, rng.normal(size=400), min_count=1)
    best = model.best_unmade(10)
    values = [v for _c, v in best]
    assert len(best) == 10 and values == sorted(values, reverse=True)
    assert all(c not in model.observed for c, _v in best)


def test_sar_matrix() -> None:
    rows = ["H", "Cl", "F", "Cl", "H", "F", "C", ""]
    cols = ["acid", "acid", "acid", "amide", "amide", "amide", "amide", "acid"]
    act = np.array([5, 6, 5.5, 7, 6, 6.2, 6.6, 9])
    m = build_matrix(rows, cols, act)
    assert m.row_labels[0] == "H" and m.dropped == 1
    assert m.cell("", "Cl", "amide").value("mean") == 7
    assert m.cell("", "C", "acid") is None
    assert m.predicted("", "C", "acid") == pytest.approx(5.7)
    paged = build_matrix(rows, cols, act, pages=["x"] * 4 + ["y"] * 4, predict=False)
    assert paged.page_labels == ["x", "y"] and paged.model is None
    with pytest.raises(ValueError, match="at most"):
        build_matrix([str(i) for i in range(70)], ["a"] * 70, np.zeros(70))


def test_mmp() -> None:
    pieces = fragment(Chem.MolFromSmiles("Clc1ccccc1C(=O)O"), 4)
    assert ("O=C(O)c1ccccc1[*:1]", "Cl[*:1]") in pieces
    assert any(change == "[*:1][H]" for _c, change in pieces)
    smiles = [
        "c1ccccc1C(=O)O", "Clc1ccccc1C(=O)O", "Fc1ccccc1C(=O)O",
        "Clc1ccc(cc1)C(=O)N", "c1ccc(cc1)C(=O)N", "Fc1ccc(cc1)C(=O)N", None,
    ]  # fmt: skip
    act = np.array([5, 6, 5.5, 7, 6, 6.2, 1])
    pairs = find_pairs(smiles, act, max_variable=3)
    found = {(readable(t.from_group), readable(t.to_group)): t for t in transforms(pairs)}
    h_cl = found[("H", "*Cl")]  # oriented from hydrogen
    assert h_cl.count == 2 and h_cl.mean == pytest.approx(1.0) and h_cl.ups == 2
    assert h_cl.rows() == {0, 1, 3, 4}
    assert ("*Cl", "H") not in found
    assert all(6 not in (p.a, p.b) for p in pairs)


def test_mmp_in_workers() -> None:
    from chemsuite.core.jobs.parallel import process_pool

    smiles = ["c1ccccc1C(=O)O", "Clc1ccccc1C(=O)O", "Fc1ccccc1C(=O)O"] * 3
    act = np.arange(9.0)
    serial = find_pairs(smiles, act)
    with process_pool(2) as pool:
        assert find_pairs(smiles, act, executor=pool) == serial


def test_sali() -> None:
    smiles = ["c1ccccc1CCO", "c1ccccc1CCN", "c1ccccc1CCC", "CCCCCCCCCC"]
    act = np.array([5.0, 8.0, 5.2, np.nan])
    cliffs = sali_pairs(compute_fingerprints(smiles), act, min_similarity=0.3)
    assert (cliffs.a < cliffs.b).all() and 3 not in cliffs.b
    assert (np.diff(cliffs.sali) <= 0).all()
    top = (int(cliffs.a[0]), int(cliffs.b[0]))
    assert 1 in top  # the 8.0 compound makes the cliff
    expected = abs(cliffs.delta[0]) / (1 - cliffs.similarity[0])
    assert cliffs.sali[0] == pytest.approx(expected)
    best, partner = cliffs.per_row(4)
    assert partner[0] == 1 and math.isnan(best[3]) and partner[3] == -1


def test_sali_skips_identical_fingerprints() -> None:
    smiles = ["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "CC(N)C(=O)N"]  # enantiomers: one fingerprint
    act = np.array([5.0, 7.0, 6.0])
    fps = compute_fingerprints(smiles)
    every = sali_pairs(fps, act, min_similarity=0.1)
    assert (int(every.a[0]), int(every.b[0])) == (0, 1) and every.similarity[0] == 1.0
    kept = sali_pairs(fps, act, min_similarity=0.1, skip_identical=True)
    assert len(kept) and (kept.similarity < 1.0).all()


def test_hansch(tmp_path) -> None:
    constants_file = tmp_path / "constants.csv"
    constants_file.write_text(
        "smiles,pi,sigma_m,sigma_p,f,r,es,mr,upsilon,sigma_star\n"
        "# a comment\n"
        "H,0,0,0,0,0,1.24,1.03,,\n"
        "*Cl,0.71,0.37,0.23,,,0.27,6.03,,\n"
        "[*]C,0.56,-0.07,-0.17,,,0,5.65,,\n"
        "CO[*:1],-0.02,0.12,-0.27,,,0.69,7.87,,\n"
        "*F,0.14,0.34,0.06,,,0.78,0.92,,\n"
    )
    constants = load_constants(constants_file)
    assert set(constants) == {"H", "*Cl", "*C", "*OC", "*F"}
    assert table_key("C[*:1]") == "*C"
    smiles, _ = _series(skip=())
    d = decompose([Chem.MolFromSmiles(s) for s in smiles], parse_core(CORE))
    pi = {"H": 0, "*Cl": 0.71, "*C": 0.56, "*OC": -0.02, "*F": 0.14}
    sigma = {"H": 0, "*Cl": 0.23, "*C": -0.17, "*OC": -0.27, "*F": 0.06}
    act = np.array(
        [
            4
            + 1.5 * sum(pi.get(d.labels[p][i], 0) for p in d.positions)
            - 0.8 * sum(sigma.get(d.labels[p][i], 0) for p in d.positions)
            for i in range(len(d))
        ]
    )
    model = hansch(d, act, constants, HanschSettings(), activity_name="pIC50")
    assert model.skipped["substituent missing from the constants"] == 4  # *N has no entry
    assert model.fit.coefficients == pytest.approx(
        [1.5, -0.8]
    ) and model.fit.intercept == pytest.approx(4)
    assert model.equation().startswith("pIC50 = 4 + 1.5·Σπ − 0.8·Σσ")
    per = hansch(d, act, constants, HanschSettings(per_position=True, steric="mr"))
    assert per.fit.names[:3] == ("π(R1)", "σ(R1)", "MR(R1)")
    with pytest.raises(ValueError):
        hansch(d, act, constants, HanschSettings(pi=False, electronic=""))
    bad = tmp_path / "bad.csv"
    bad.write_text("name,pi\nx,1\n")
    with pytest.raises(ValueError, match="header"):
        load_constants(bad)


def test_cancel_keeps_the_rows_decomposed(monkeypatch) -> None:
    from chemsuite.core.jobs import Cancelled, CancelToken
    from chemsuite.sar import rgroups

    monkeypatch.setattr(rgroups, "CHUNK", 4)
    smiles, _ = _series()
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    token = CancelToken()

    def progress(done: int, total: int) -> None:
        if done >= 8:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        decompose(mols, parse_core(CORE), progress=progress, cancel=token)
    d = stopped.value.partial
    assert len(d) == 8 and d.matched.all() and d.positions == ("R1", "R2")
