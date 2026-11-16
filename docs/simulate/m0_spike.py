"""M0 feasibility gates for docs/SIMULATE_PLAN.md. Not app code: run it in the OpenMM env.

    python -I m0_spike.py [gate ...]      # gates: env charges gb templates massless

Each gate prints one JSON object on a line starting with ``GATE``; a failure prints its
traceback and the run goes on to the next gate.
"""

# ruff: noqa: BLE001  (a spike records whatever a library raises and goes on)

from __future__ import annotations

import json
import sys
import traceback

import numpy as np

LIGANDS = {
    "ibuprofen anion": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)[O-]",
    "benzamidinium": "NC(=[NH2+])c1ccccc1",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1",
}
SAGE = "openff-2.3.0"
PEPTIDE = "AHKDEYCWSR"


def gate_env() -> dict:
    import importlib.metadata as md
    import platform

    import openmm

    out: dict = {"python": sys.version.split()[0], "machine": platform.platform()}
    for name in ("openmm", "openmmforcefields", "openff-toolkit", "openff-interchange",
                 "openff-nagl", "openff-nagl-models", "rdkit", "numpy", "ambertools"):  # fmt: skip
        try:
            out[name] = md.version(name)
        except md.PackageNotFoundError:
            out[name] = None
    out["openmm"] = openmm.version.version
    out["platforms"] = [
        openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())
    ]
    out["plugin_load_failures"] = list(openmm.Platform.getPluginLoadFailures())
    from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

    out["toolkits"] = [type(t).__name__ for t in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits]
    from openff.toolkit import ForceField

    ff = ForceField(f"{SAGE}.offxml")
    out["sage_handlers"] = sorted(ff.registered_parameter_handlers)
    for tag in ("NAGLCharges", "ToolkitAM1BCC"):
        if tag in ff.registered_parameter_handlers:
            out[tag] = ff.get_parameter_handler(tag).to_dict()
    return out


def ligand(smiles: str):
    from openff.toolkit import Molecule

    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    mol.generate_conformers(n_conformers=1)
    return mol


def nagl_charges(mol, model: str):
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    mol.assign_partial_charges(model, toolkit_registry=NAGLToolkitWrapper())
    return np.array(mol.partial_charges.m_as("elementary_charge"))


def system_charges(system) -> np.ndarray:
    import openmm
    from openmm import unit

    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            return np.array(
                [
                    force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                    for i in range(force.getNumParticles())
                ]
            )
    raise RuntimeError("no NonbondedForce")


def nagl_model() -> str:
    from openff.toolkit import ForceField

    ff = ForceField(f"{SAGE}.offxml")
    if "NAGLCharges" in ff.registered_parameter_handlers:
        return str(ff.get_parameter_handler("NAGLCharges").model_file)
    from openff.nagl_models import list_available_nagl_models

    return sorted(p.name for p in list_available_nagl_models())[-1]


def gate_charges() -> dict:
    """Does the template generator use charges it is handed, and what does it do without?"""
    import time

    from openmm import app
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator

    model = nagl_model()
    out: dict = {"model": model, "ligands": {}}
    for name, smiles in LIGANDS.items():
        row: dict = {}
        mol = ligand(smiles)
        t = time.perf_counter()
        q = nagl_charges(mol, model)
        row["nagl_seconds"] = round(time.perf_counter() - t, 2)
        row["net"] = round(float(q.sum()), 4)
        # Marked charges: shifted by a pattern that keeps the total, so "used as given" shows.
        marked = q + 0.01 * (np.arange(len(q)) % 2 * 2 - 1)
        marked -= (marked.sum() - q.sum()) / len(q)
        from openff.units import unit as off_unit

        given = ligand(smiles)
        given.partial_charges = marked * off_unit.elementary_charge
        top = given.to_topology().to_openmm()
        ff = app.ForceField()
        ff.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=[given], forcefield=SAGE).generator
        )
        got = system_charges(ff.createSystem(top))
        row["given_used_max_abs_diff"] = float(np.abs(got - marked).max())
        # No charges handed in: what does openmmforcefields produce on its own?
        bare = ligand(smiles)
        ff = app.ForceField()
        ff.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=[bare], forcefield=SAGE).generator
        )
        try:
            t = time.perf_counter()
            own = system_charges(ff.createSystem(bare.to_topology().to_openmm()))
            row["own_seconds"] = round(time.perf_counter() - t, 2)
            row["own_vs_nagl_max_abs_diff"] = float(np.abs(own - q).max())
        except Exception as error:  # a spike: record whatever it raises
            row["own_error"] = f"{type(error).__name__}: {error}"[:500]
        out["ligands"][name] = row
    return out


def peptide(variants: dict[int, str] | None = None, forcefield: str = "amber14/protein.ff14SB.xml"):
    """A capped-free peptide with hydrogens from Modeller: (topology, positions, ForceField)."""
    import io

    from openmm import app
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSequence(PEPTIDE))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    heavy = Chem.RemoveHs(mol)
    block = Chem.MolToPDBBlock(heavy, flavor=2)  # no CONECT
    pdb = app.PDBFile(io.StringIO(block))
    ff = app.ForceField(forcefield)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    wanted = [None] * modeller.topology.getNumResidues()
    for i, v in (variants or {}).items():
        wanted[i] = v
    modeller.addHydrogens(ff, variants=wanted)
    return modeller.topology, modeller.positions, ff


def without_atoms(topology, positions, drop: set[tuple[int, str]]):
    """A copy lacking the (residue index, atom name) atoms."""
    from openmm import app

    modeller = app.Modeller(topology, positions)
    gone = [a for a in modeller.topology.atoms() if (a.residue.index, a.name) in drop]
    if len(gone) != len(drop):
        raise RuntimeError(f"wanted to drop {sorted(drop)}, found {[a.name for a in gone]}")
    modeller.delete(gone)
    return modeller.topology, modeller.positions


def gate_templates() -> dict:
    """Which of Prepare's states have an Amber template, in ff14SB and ff19SB."""
    from openmm import app

    names = ["HID", "HIE", "HIP", "HIS", "HIN", "HIM", "ASH", "GLH", "LYN", "LYS", "CYS", "CYM",
             "CYX", "TYR", "TYM", "ARG", "ARN", "ACE", "NME", "NHE", "NH2"]  # fmt: skip
    # PEPTIDE residues by index: A0 H1 K2 D3 E4 Y5 C6 W7 S8 R9
    cases = {
        "default (HIE etc.)": ({}, set()),
        "HID": ({1: "HID"}, set()),
        "HIP": ({1: "HIP"}, set()),
        "ASH": ({3: "ASH"}, set()),
        "GLH": ({4: "GLH"}, set()),
        "LYN": ({2: "LYN"}, set()),
        "CYM (CYS minus HG)": ({}, {(6, "HG")}),
        "TYM (TYR minus HH)": ({}, {(5, "HH")}),
        "neutral N-terminus (minus H3)": ({}, {(0, "H3")}),
        "imidazolate (HIE minus HE2)": ({1: "HIE"}, {(1, "HE2")}),
    }
    out: dict = {}
    for label, xml in (("ff14SB", "amber14/protein.ff14SB.xml"),
                       ("ff19SB", "amber19/protein.ff19SB.xml")):  # fmt: skip
        row: dict = {}
        try:
            ff = app.ForceField(xml)
        except Exception as error:  # a spike: record whatever it raises
            out[label] = {"error": f"{type(error).__name__}: {error}"[:300]}
            continue
        have = set(ff._templates)
        row["count"] = len(have)
        row["present"] = {n: [p + n for p in ("", "N", "C") if p + n in have] for n in names}
        row["other"] = sorted(n for n in have if len(n) > 4 or n.startswith(("N", "C")) is False)[
            :80
        ]
        matches: dict = {}
        for case, (variants, drop) in cases.items():
            try:
                top, pos, _ = peptide(variants, xml)
                if drop:
                    top, pos = without_atoms(top, pos, drop)
                bad = ff.getUnmatchedResidues(top)
                matches[case] = "ok" if not bad else f"unmatched: {[r.name for r in bad]}"
                if not bad:
                    templates = ff.getMatchingTemplates(top)
                    idx = {1, 2, 3, 4, 5, 6, 0} & ({k for k in variants} | {r for r, _ in drop})
                    matches[case] += " " + ",".join(templates[i].name for i in sorted(idx))
            except Exception as error:  # a spike: record whatever it raises
                matches[case] = f"{type(error).__name__}: {error}"[:300]
        row["matches"] = matches
        out[label] = row
    return out


def complex_system(solvent_xml: str | None, charged: bool = True):
    """The peptide with a NAGL-charged Sage ligand 2.5 nm away, through SystemGenerator."""
    from openmm import Vec3, app, unit
    from openmmforcefields.generators import SystemGenerator

    top, pos, _ = peptide()
    mol = ligand(LIGANDS["ibuprofen anion"])
    if charged:
        nagl_charges(mol, nagl_model())
    modeller = app.Modeller(top, pos)
    shift = Vec3(2.5, 0, 0) * unit.nanometer
    lig_pos = [p + shift for p in mol.conformers[0].to_openmm()]
    modeller.add(mol.to_topology().to_openmm(), lig_pos)
    forcefields = ["amber/ff14SB.xml"] + ([solvent_xml] if solvent_xml else [])
    generator = SystemGenerator(
        forcefields=forcefields,
        small_molecule_forcefield=SAGE,
        molecules=[mol],
        forcefield_kwargs={"constraints": None, "rigidWater": False},
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
    )
    system = generator.create_system(modeller.topology)
    return system, modeller.topology, modeller.positions, mol


def gate_gb() -> dict:
    """Does implicit/gbn2.xml give GB parameters to a SMIRNOFF ligand through SystemGenerator?"""
    import openmm
    from openmm import unit

    out: dict = {}
    for label, xml in (("vacuum", None), ("gbn2", "implicit/gbn2.xml"),
                       ("obc2", "implicit/obc2.xml")):  # fmt: skip
        row: dict = {}
        try:
            system, top, pos, mol = complex_system(xml)
        except Exception as error:  # a spike: record whatever it raises
            out[label] = {"error": f"{type(error).__name__}: {error}"[:1500]}
            continue
        row["forces"] = [type(f).__name__ for f in system.getForces()]
        n_lig = mol.n_atoms
        n = system.getNumParticles()
        for force in system.getForces():
            if isinstance(force, openmm.CustomGBForce | openmm.GBSAOBCForce):
                row["gb_particles"] = force.getNumParticles()
                row["system_particles"] = n
                params = [force.getParticleParameters(i) for i in range(n - n_lig, n)]
                # the parameters after the charge, per ligand element: one set each means the
                # radii were assigned by element, as GBn2 intends
                symbols = [a.element.symbol for a in top.atoms()][n - n_lig :]
                by_element: dict[str, set] = {}
                for symbol, p in zip(symbols, params, strict=True):
                    by_element.setdefault(symbol, set()).add(
                        tuple(round(float(x), 5) for x in p[1:])
                    )
                row["ligand_gb_params_by_element"] = {
                    k: sorted(v) for k, v in by_element.items()
                }  # fmt: skip
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CPU"))
        context.setPositions(pos)
        e0 = context.getState(getEnergy=True).getPotentialEnergy()
        openmm.LocalEnergyMinimizer.minimize(context, 10.0, 500)
        e1 = context.getState(getEnergy=True).getPotentialEnergy()
        row["energy_kj"] = [round(e0.value_in_unit(unit.kilojoule_per_mole), 1),
                            round(e1.value_in_unit(unit.kilojoule_per_mole), 1)]  # fmt: skip
        out[label] = row
    return out


def gate_massless() -> dict:
    """Does LocalEnergyMinimizer leave zero-mass particles where they are?"""
    import openmm
    from openmm import app, unit

    top, pos, ff = peptide()
    out: dict = {}
    elements = [a.element.symbol for a in top.atoms()]
    heavy = np.array([e != "H" for e in elements])
    start = np.array(pos.value_in_unit(unit.nanometer))
    rng = np.random.default_rng(3)
    start = start + rng.normal(0, 0.01, start.shape)  # strain everything a little
    platforms = [
        openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())
    ]
    for constraints in (None, app.HBonds):
        for name in platforms:
            key = f"{name}, constraints={'HBonds' if constraints else 'None'}"
            try:
                system = ff.createSystem(top, nonbondedMethod=app.NoCutoff,
                                         constraints=constraints)  # fmt: skip
                for i in np.flatnonzero(heavy):
                    system.setParticleMass(int(i), 0.0)
                context = openmm.Context(system, openmm.VerletIntegrator(0.001),
                                         openmm.Platform.getPlatformByName(name))  # fmt: skip
                context.setPositions(start)
                e0 = context.getState(getEnergy=True).getPotentialEnergy()
                openmm.LocalEnergyMinimizer.minimize(context, 1.0, 2000)
                state = context.getState(getEnergy=True, getPositions=True)
                end = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                moved = np.linalg.norm(end - start, axis=1)
                out[key] = {
                    "heavy_max_move_nm": float(moved[heavy].max()),
                    "hydrogen_max_move_nm": float(moved[~heavy].max()),
                    "energy_kj": [
                        round(e0.value_in_unit(unit.kilojoule_per_mole), 1),
                        round(state.getPotentialEnergy().value_in_unit(
                            unit.kilojoule_per_mole), 1),
                    ],
                }  # fmt: skip
            except Exception as error:  # a spike: record whatever it raises
                out[key] = f"{type(error).__name__}: {error}"[:300]
    return out


GATES = {
    "env": gate_env,
    "charges": gate_charges,
    "templates": gate_templates,
    "gb": gate_gb,
    "massless": gate_massless,
}


def main(argv: list[str]) -> int:
    failed = 0
    for name in argv or list(GATES):
        try:
            print("GATE", json.dumps({name: GATES[name]()}, default=str), flush=True)
        except Exception:  # a spike: report and go on
            failed += 1
            print(f"GATE-FAILED {name}", flush=True)
            traceback.print_exc()
    return failed


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
