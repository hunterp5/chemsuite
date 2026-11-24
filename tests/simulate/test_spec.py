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

import json

import numpy as np
import pytest

from chemsuite.simulate import build, spec
from chemsuite.simulate.result import CLASH, Contact, MinimiseResult
from chemsuite.simulate.spec import ForceFieldSpec, MinimiseSpec, SpecError, SystemSpec

from .test_build import _without


@pytest.fixture
def system(lysozyme) -> SystemSpec:
    s, states = _without(lysozyme, ("PO4", "BME"))
    return build.system_spec(s, states)


def test_a_spec_survives_plain_json(system) -> None:
    text = json.dumps(system.to_json())  # nothing but plain values
    again = SystemSpec.from_json(json.loads(text))
    assert again == system and again.digest() == system.digest()


def test_charges_belong_to_the_topology_and_leave_its_digest_alone(system) -> None:
    (molecule,) = system.molecules
    n = system.residues[molecule.residue].count
    charged = system.with_charges(molecule.residue, tuple([0.01] * n), "openff-gnn-am1bcc-1.0.0")
    assert charged.charged and charged.digest() == system.digest()
    assert SystemSpec.from_json(charged.to_json()) == charged
    with pytest.raises(SpecError):
        system.with_charges(molecule.residue, (0.0,), "x")  # one per atom
    with pytest.raises(SpecError):
        system.with_charges(0, tuple([0.0] * system.residues[0].count), "x")  # not a molecule
    moved = system.to_json()
    moved["atoms"]["name"][0] = "X"
    assert SystemSpec.from_json(moved).digest() != system.digest()


def test_damaged_json_is_refused_with_a_reason(system) -> None:
    good = system.to_json()
    for damage in (
        lambda d: d.update(schema=99),
        lambda d: d["atoms"].update(uid=d["atoms"]["uid"][:-1]),
        lambda d: d["atoms"].update(element=["C"] * len(d["atoms"]["element"])),
        lambda d: d["residues"]["first"].__setitem__(1, 0),
        lambda d: d["residues"]["kind"].__setitem__(0, "planet"),
        lambda d: d["bonds"].append([0, 10**9, 1]),
        lambda d: d["molecules"].append({"residue": -1, "molblock": ""}),
        lambda d: d.pop("atoms"),
    ):
        broken = json.loads(json.dumps(good))
        damage(broken)
        with pytest.raises(SpecError):
            SystemSpec.from_json(broken)
    with pytest.raises(SpecError):
        SystemSpec.from_json([])


def test_force_fields_and_protocol_round_trip_and_refuse_nonsense() -> None:
    assert ForceFieldSpec() == ForceFieldSpec(spec.FF14SB, spec.SAGE, spec.TIP3P, spec.GBN2)
    assert ForceFieldSpec.paired(spec.FF19SB).water == spec.OPC
    for ff in (ForceFieldSpec(), ForceFieldSpec.paired(spec.FF19SB)):
        assert ForceFieldSpec.from_json(json.loads(json.dumps(ff.to_json()))) == ff
    with pytest.raises(SpecError):
        ForceFieldSpec(protein="charmm36")
    with pytest.raises(SpecError):
        ForceFieldSpec.from_json({"protein": "ff14SB"})
    protocol = MinimiseSpec()
    assert [s.moving for s in protocol.stages] == [spec.HYDROGENS, spec.HEAVY, spec.HEAVY]
    assert protocol.stages[0].restraints() == (0.0,) * len(spec.RESTRAINED)
    side, whole = protocol.stages[1], protocol.stages[2]
    assert side.backbone > side.side_chains > 0 and whole.ligands < whole.backbone
    assert MinimiseSpec.from_json(json.loads(json.dumps(protocol.to_json()))) == protocol
    # a metal's contacts held in every stage that moves heavy atoms, harder than any atom
    assert protocol.stages[0].contacts == 0
    assert all(s.contacts > max(s.restraints()) for s in (side, whole))
    older = protocol.to_json()
    for stage in older["stages"]:
        del stage["contacts"]  # a protocol saved before metal contacts were held
    assert all(s.contacts == 0 for s in MinimiseSpec.from_json(older).stages)
    halved = spec.MinimiseSettings(restraints=0.5).protocol()
    assert halved.stages[2].contacts == whole.contacts / 2
    with pytest.raises(SpecError):
        spec.Stage("x", spec.HEAVY, ions=-1.0)
    with pytest.raises(SpecError):
        spec.Stage("x", spec.HEAVY, contacts=-1.0)
    with pytest.raises(SpecError):
        MinimiseSpec(tighten=-1)
    with pytest.raises(SpecError):
        MinimiseSpec(stages=())
    with pytest.raises(SpecError):
        spec.Stage("x", "everything")
    with pytest.raises(SpecError):
        MinimiseSpec.from_json({"stages": [{"name": "x"}]})


def _result(system: SystemSpec, s, shift: float = 0.1) -> MinimiseResult:
    xyz = build.positions(s, system) + shift
    return MinimiseResult(
        digest=system.digest(), uid=system.uid, xyz=tuple(map(tuple, xyz.tolist())),
        energy_before=-1000.0, energy_after=-5000.0, rmsd=0.12, rmsd_pocket=0.1,
        rmsd_ligand=0.2, contact_before=Contact("A/LEU 99/CD1", "A/JZ4 167/C1", 1.9),
        contact_after=Contact("A/LEU 99/CD1", "A/JZ4 167/C1", 3.1), max_force=8.5,
        max_force_atom="A/LYS 16/NZ", moved_most=(("A/LYS 16", 0.8),),
        charges=(("A/JZ4 167", "openff-gnn-am1bcc-1.0.0", (0.1, -0.1)),),
        versions=(("openmm", "8.4.0"),), platform="CUDA", seconds=(("Hydrogens", 0.5),),
    )  # fmt: skip


def test_a_result_lands_on_its_atoms_by_uid(system, lysozyme) -> None:
    s, _states = _without(lysozyme, ("PO4", "BME"))
    result = _result(system, s)
    again = MinimiseResult.from_json(json.loads(json.dumps(result.to_json())))
    assert again.uid == result.uid and again.digest == result.digest
    assert np.allclose(again.xyz, result.xyz, atol=1e-4)
    assert again.contact_after == result.contact_after and again.charges == result.charges
    xyz = again.coordinates(s)
    assert np.allclose(xyz, s.xyz() + 0.1, atol=1e-3)  # every atom is in the system here
    # an atom the minimisation did not hold keeps its place
    part = MinimiseResult(result.digest, result.uid[:10], result.xyz[:10], 0.0, 0.0, 0.0)
    moved = part.coordinates(s)
    rest = np.ones(s.n_atoms, dtype=bool)
    rest[s.indices(np.array(result.uid[:10]))] = False
    assert np.allclose(moved[rest], s.xyz()[rest])
    gone = MinimiseResult(result.digest, (10**9,), ((0.0, 0.0, 0.0),), 0.0, 0.0, 0.0)
    with pytest.raises(SpecError):
        gone.coordinates(s)


def test_a_result_says_what_deserves_a_warning(system, lysozyme) -> None:
    s, _states = _without(lysozyme, ("PO4", "BME"))
    result = _result(system, s)
    assert result.warnings(0.3) == []
    assert "moved 0.12" in result.warnings(0.1)[0]
    clash = MinimiseResult(result.digest, (), (), 0.0, 0.0, 0.5,
                           contact_after=Contact("a", "b", CLASH - 0.5))  # fmt: skip
    assert [w.split(":")[0] for w in clash.warnings(1.0)] == ["A clash is left"]


def test_a_damaged_result_is_refused() -> None:
    with pytest.raises(SpecError):
        MinimiseResult.from_json({"uid": [1], "xyz": [0.0, 0.0]})
    with pytest.raises(SpecError):
        MinimiseResult.from_json({"digest": "x", "uid": [1], "xyz": [0, 0, 0]})  # no energies
    with pytest.raises(SpecError):
        MinimiseResult.from_json("result")


def test_minimise_settings_make_the_force_field_and_the_protocol() -> None:
    from chemsuite.simulate.spec import STAGES, MinimiseSettings

    plain = MinimiseSettings()
    assert plain.forcefield() == spec.ForceFieldSpec()
    assert plain.protocol() == spec.MinimiseSpec()  # the defaults are M3's protocol
    firm = MinimiseSettings(protein=spec.FF19SB, solvent=spec.VACUUM, restraints=2.0,
                            rmsd_limit=0.5)  # fmt: skip
    assert firm.forcefield() == spec.ForceFieldSpec(spec.FF19SB, water=spec.OPC,
                                                    solvent=spec.VACUUM)  # fmt: skip
    protocol = firm.protocol()
    assert protocol.rmsd_limit == 0.5
    for stage, base in zip(protocol.stages, STAGES, strict=True):
        assert stage.restraints() == tuple(2.0 * v for v in base.restraints())
    assert firm.problem(has_waters=True) and not firm.problem(has_waters=False)
    assert not plain.problem(has_waters=True)
    with pytest.raises(spec.SpecError):
        MinimiseSettings(solvent=spec.EXPLICIT)


def test_minimise_settings_survive_their_saved_form() -> None:
    from chemsuite.simulate.spec import MinimiseSettings

    settings = MinimiseSettings(spec.FF19SB, spec.VACUUM, 0.5, 0.4, "CUDA")
    assert MinimiseSettings.from_json(settings.to_json()) == settings
    assert MinimiseSettings.from_json(None) == MinimiseSettings()
    junk = {"protein": "ff99", "solvent": 3, "restraints": -1, "rmsd_limit": True}
    assert MinimiseSettings.from_json(junk) == MinimiseSettings()
