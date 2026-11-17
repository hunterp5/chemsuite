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

"""The script that runs in OpenMM's own interpreter, and the requests it answers.

It is passed with ``python -I -c`` (through ``wsl.exe -e`` for a WSL interpreter), so nothing is
written to disk and the interpreter sees only its own environment. It speaks the JSON-lines
protocol of :mod:`chemsuite.core.jobs.lines`; systems arrive as :class:`~.spec.SystemSpec` JSON
and positions as flat lists in ångströms.

Requests
    ``hello``
        the Python, OpenMM, OpenFF and force-field versions, and each OpenMM platform with
        whether a context can be made on it (a CUDA platform whose driver is too old is listed by
        OpenMM but fails there).
    ``charges`` ``{"system"}``
        partial charges for each small molecule that has none, from the NAGL model Sage names:
        ``{"charges": [{"residue", "charges", "model"}]}``. Computed once and stored, so every
        later system is built from the same ones.
    ``build`` ``{"system", "forcefield"}``
        the OpenMM system: ``{"particles"}``, or ``ok: false`` with ``problems`` (``label``,
        ``text``). Every protein residue and cap is given its template by name (never matched by
        graph alone), and every residue's charge in the system must be the spec's formal charge.
        Built systems are kept, by the whole request, charges included.
    ``minimise`` ``{"system", "forcefield", "positions", "protocol", "platform", "free"}``
        the staged minimisation of :class:`~.spec.MinimiseSpec`, each class of heavy atom held to
        where it started but those in ``free`` (atoms Prepare built); the last stage run again
        with its restraints doubled while the held atoms end further than the limit. Sends
        ``{"stage", "iteration", "energy"}`` events as it goes and stops between iterations when
        asked; ``stages`` in the reply holds every run, and ``rmsd`` the held atoms' at the end.

Replies go to a duplicate of the original stdout, and file descriptor 1 is pointed at stderr
before anything is imported, so that a library's prints cannot corrupt the stream.
"""

from __future__ import annotations

#: ``python -I -c SCRIPT``: isolated, so the child's imports stay in its own environment
ARGS = ("-I", "-c")

SCRIPT = r'''
"""chemsuite's OpenMM worker. One JSON request per line; see chemsuite/simulate/worker.py."""
import copy
import hashlib
import json
import os
import queue
import sys
import threading
import time
import traceback

_replies = os.fdopen(os.dup(1), "w", encoding="utf-8", newline="\n")
os.dup2(2, 1)
sys.stdout = sys.stderr
_writing = threading.Lock()

NM = 10.0  # angstroms in a nanometre
PROTEIN_FILES = {"ff14SB": "amber14/protein.ff14SB.xml", "ff19SB": "amber19/protein.ff19SB.xml"}
WATER_FILES = {"tip3p": "amber14/tip3p.xml", "opc": "amber19/opc.xml"}
SOLVENT_FILES = {"gbn2": ["implicit/gbn2.xml"], "vacuum": []}
PREFERRED = ("CUDA", "OpenCL", "CPU", "Reference")
BACKBONE = frozenset(["N", "CA", "C", "O", "OXT"])
#: the restraint classes, in the order of spec.RESTRAINED, by residue kind (protein: side chain)
RESTRAINED = ("backbone", "side_chains", "ligands", "waters", "ions")
CLASS_OF_KIND = {"protein": 1, "cap": 1, "molecule": 2, "water": 3, "ion": 4}
EVENT_EVERY = 25  # minimiser iterations between energy events
KEEP = 2  # built systems kept


def _send(payload):
    with _writing:
        _replies.write(json.dumps(payload) + "\n")
        _replies.flush()


class Problems(Exception):
    """What stops a system from being built: a list of (label, text)."""

    def __init__(self, found):
        super().__init__("; ".join(f"{label}: {text}" for label, text in found))
        self.found = found


def _expected():
    """The errors a request can end with and leave the worker up. Anything else ends it, and the
    app reports the traceback."""
    found = [Problems, ValueError, KeyError, IndexError, TypeError, RuntimeError,
             ArithmeticError, OSError, ImportError]
    try:
        from openmm import OpenMMException
        found.append(OpenMMException)
    except ImportError:
        pass
    try:
        from openff.toolkit.utils.exceptions import OpenFFToolkitException
        found.append(OpenFFToolkitException)
    except ImportError:
        pass
    return tuple(found)


# -- reading requests -------------------------------------------------------------------------
class Inbox:
    """Requests from stdin, read on a thread of their own so a stop reaches a running one."""

    def __init__(self):
        self.requests = queue.Queue()
        self.stop = threading.Event()
        self.current = None
        threading.Thread(target=self._read, daemon=True).start()

    def _read(self):
        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue
            try:
                request = json.loads(line)
            except ValueError:
                _send({"id": None, "ok": False, "error": "could not read the request"})
                continue
            if not isinstance(request, dict):
                continue
            if request.get("op") == "stop":
                if request.get("id") == self.current:
                    self.stop.set()
                continue
            self.requests.put(request)
            if request.get("op") == "quit":
                return
        self.stop.set()  # stdin closed: the app is gone, so end what is running and quit
        self.requests.put(None)


# -- hello --------------------------------------------------------------------------------------
_PLATFORMS = None


def platforms():
    """Every OpenMM platform: ``[name, works, detail]``; detail is the device or the error."""
    global _PLATFORMS
    if _PLATFORMS is not None:
        return _PLATFORMS
    import openmm

    found = []
    for i in range(openmm.Platform.getNumPlatforms()):
        platform = openmm.Platform.getPlatform(i)
        name = platform.getName()
        system = openmm.System()
        system.addParticle(1.0)
        try:
            context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
        except openmm.OpenMMException as err:
            found.append([name, False, str(err).splitlines()[0][:200]])
            continue
        detail = ""
        for key in ("DeviceName", "OpenCLDeviceName", "CudaDeviceName"):
            if key in platform.getPropertyNames():
                detail = platform.getPropertyValue(context, key)
                break
        found.append([name, True, detail])
        del context
    _PLATFORMS = found
    return found


def hello():
    import importlib.metadata as metadata

    import openmm

    versions = {"openmm": openmm.version.version}
    for name in ("openmmforcefields", "openff-toolkit", "openff-nagl", "openff-nagl-models",
                 "rdkit", "numpy"):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = ""
    return {"python": sys.version.split()[0], "versions": versions, "platforms": platforms(),
            "plugin_failures": list(openmm.Platform.getPluginLoadFailures())}


# -- charges ------------------------------------------------------------------------------------
def nagl_model(small_molecule):
    from openff.toolkit import ForceField

    handler = ForceField(small_molecule + ".offxml").get_parameter_handler("NAGLCharges")
    return str(handler.model_file)


def off_molecule(molblock):
    """The OpenFF molecule of a molblock, its atoms in the molblock's order."""
    from openff.toolkit import Molecule
    from rdkit import Chem

    mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError("a molblock RDKit cannot read")
    return Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)


def charges(request):
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    system = request["system"]
    model = nagl_model(request.get("small_molecule") or "openff-2.3.0")
    out = []
    for molecule in system["molecules"]:
        if molecule.get("charges"):
            continue
        off = off_molecule(molecule["molblock"])
        off.assign_partial_charges(model, toolkit_registry=NAGLToolkitWrapper())
        values = [float(q) for q in off.partial_charges.m_as("elementary_charge")]
        name = model[:-3] if model.endswith(".pt") else model
        out.append({"residue": molecule["residue"], "charges": values, "model": name})
    return {"charges": out}


# -- build --------------------------------------------------------------------------------------
class Built:
    def __init__(self, system, topology, labels, residue_of, heavy, classes):
        self.system = system
        self.topology = topology
        self.labels = labels  # residue labels
        self.residue_of = residue_of  # residue of each atom
        self.heavy = heavy  # atom indices
        self.classes = classes  # each heavy atom's restraint class (an index into RESTRAINED)


_BUILT = {}


def build(request):
    key = hashlib.sha256(json.dumps([request["system"], request["forcefield"]],
                                    sort_keys=True).encode("utf-8")).hexdigest()
    if key in _BUILT:
        return _BUILT[key]
    made = _build(request["system"], request["forcefield"])
    while len(_BUILT) >= KEEP:
        del _BUILT[next(iter(_BUILT))]
    _BUILT[key] = made
    return made


def _build(spec, ff):
    from openff.units import unit as off_unit
    from openmm import NonbondedForce, app, unit
    from openmm.app.element import Element
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator

    atoms, rows = spec["atoms"], spec["residues"]
    labels = rows["label"]
    problems = []
    if ff["solvent"] not in SOLVENT_FILES:
        problems.append(("system", "explicit solvent is not available yet"))
    files = [PROTEIN_FILES[ff["protein"]], WATER_FILES[ff["water"]],
             *SOLVENT_FILES.get(ff["solvent"], [])]
    forcefield = app.ForceField(*files)
    molecules = []
    for molecule in spec["molecules"]:
        if not molecule.get("charges"):
            problems.append((labels[molecule["residue"]], "has no partial charges yet"))
            continue
        off = off_molecule(molecule["molblock"])
        off.partial_charges = [float(q) for q in molecule["charges"]] * off_unit.elementary_charge
        molecules.append(off)
    if molecules:
        generator = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield=ff["small_molecule"])
        forcefield.registerTemplateGenerator(generator.generator)

    topology = app.Topology()
    chain = None
    residues = []
    created = []
    residue_of = []
    for r, name in enumerate(rows["name"]):
        if chain is None or rows["chain"][r] != chain.id:
            chain = topology.addChain(rows["chain"][r])
        residue = topology.addResidue(name, chain, str(rows["seq"][r]), rows["icode"][r])
        residues.append(residue)
        for i in range(rows["first"][r], rows["first"][r] + rows["count"][r]):
            created.append(topology.addAtom(atoms["name"][i],
                                            Element.getByAtomicNumber(atoms["element"][i]),
                                            residue))
            residue_of.append(r)
    for i, j, _order in spec["bonds"]:
        topology.addBond(created[i], created[j])

    # every residue but a small molecule takes its template by name: an ion's by name too, since
    # FE (Fe3+) and FE2 are the same single atom to a graph match
    # A template is checked against its residue here, atom for atom by element, because OpenMM
    # reports a forced template that does not match with a bare Exception, which would end the
    # worker. A wrong Amber form differs in its hydrogens, so it is caught too.
    templates = {}
    for r, residue in enumerate(residues):
        if rows["kind"][r] not in ("protein", "cap", "water", "ion"):
            continue
        template = forcefield._templates.get(residue.name)
        if template is None:
            problems.append((labels[r], f"{ff['protein']}/{ff['water']} has no "
                                        f"{residue.name} template"))
            continue
        sites = [a for a in template.atoms if a.element is None]
        if sites:
            problems.append((labels[r], f"{ff['water']} water has {len(sites)} more site(s) "
                                        "than its atoms, which chemsuite does not place yet: "
                                        "use ff14SB with TIP3P, or leave the waters out"))
            continue
        wanted = sorted(a.element.atomic_number for a in template.atoms)
        have = sorted(a.element.atomic_number for a in residue.atoms())
        if wanted != have:
            problems.append((labels[r], f"its atoms do not match the {residue.name} template "
                                        f"({len(have)} atoms, the template {len(wanted)})"))
            continue
        templates[residue] = residue.name
    if problems:
        raise Problems(problems)
    try:
        system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None,
                                         rigidWater=False, residueTemplates=templates)
    except ValueError as err:
        raise Problems([("system", str(err).splitlines()[0])]) from None

    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
    totals = [0.0] * len(residues)
    for i in range(system.getNumParticles()):
        charge = nonbonded[0].getParticleParameters(i)[0]
        totals[residue_of[i]] += charge.value_in_unit(unit.elementary_charge)
    formal = [0] * len(residues)
    for i, q in enumerate(atoms["formal"]):
        formal[residue_of[i]] += q
    for r, total in enumerate(totals):
        if abs(total - formal[r]) > 0.01:
            problems.append((labels[r], f"the force field gives it charge {total:+.2f}, "
                                        f"Prepare {formal[r]:+d}"))
    if problems:
        raise Problems(problems)

    kinds = rows["kind"]
    heavy = [i for i, z in enumerate(atoms["element"]) if z != 1]
    classes = [0 if kinds[residue_of[i]] == "protein" and atoms["name"][i] in BACKBONE
               else CLASS_OF_KIND[kinds[residue_of[i]]] for i in heavy]
    return Built(system, topology, labels, residue_of, heavy, classes)


# -- minimise -----------------------------------------------------------------------------------
def _platform(name):
    import openmm

    working = [p[0] for p in platforms() if p[1]]
    if name:
        if name not in working:
            raise ValueError(f"the {name} platform does not work here")
        chosen = name
    else:
        chosen = next(p for p in PREFERRED if p in working)
    platform = openmm.Platform.getPlatformByName(chosen)
    properties = {"Precision": "mixed"} if chosen in ("CUDA", "OpenCL") else {}
    return platform, properties


def minimise(request, inbox, ident):
    import numpy
    import openmm
    from openmm import unit

    made = build(request)
    protocol = request["protocol"]
    platform, properties = _platform(request.get("platform") or "")
    flat = request["positions"]
    n = made.system.getNumParticles()
    if len(flat) != 3 * n:
        raise ValueError(f"{len(flat) // 3} positions for {n} atoms")
    start = numpy.asarray(flat, dtype=float).reshape(n, 3) / NM  # nm
    free = {int(i) for i in request.get("free", [])}
    held = [(i, c) for i, c in zip(made.heavy, made.classes) if i not in free]
    run = Run(made, protocol, platform, properties, start, held, inbox, ident)
    positions = start
    for stage in protocol["stages"]:
        if inbox.stop.is_set():
            break
        positions = run.stage(stage, positions)
    last = protocol["stages"][-1]
    began = run.began_last
    for attempt in range(protocol.get("tighten", 0)):
        if (inbox.stop.is_set() or last["moving"] == "hydrogens"
                or _rmsd(run.done, start, held) <= protocol["rmsd_limit"]):
            break
        scale = 2 ** (attempt + 1)
        kept = run.done
        run.stage(last, began, scale=scale, name=f"{last['name']}, restraints x{scale}")
        if inbox.stop.is_set():  # a restart stopped part-way is worse than the run it repeats
            run.done = kept
            run.stages[-1]["name"] += " (stopped, not used)"
            break
    positions = run.done
    final = openmm.Context(copy.deepcopy(made.system), openmm.VerletIntegrator(0.001), platform,
                           properties)
    final.setPositions(positions)
    state = final.getState(getEnergy=True, getForces=True)
    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    norms = (forces ** 2).sum(axis=1) ** 0.5
    worst = int(norms.argmax())
    xyz = (positions * NM).reshape(-1).tolist()
    return {"positions": xyz, "stages": run.stages, "energy_before": run.first_energy,
            "energy_after": _energy(final), "max_force": float(norms[worst]),
            "max_force_atom": worst, "platform": platform.getName(),
            "stopped": inbox.stop.is_set(), "rmsd": _rmsd(positions, start, held)}


class Run:
    """The stages of one minimisation. :attr:`done` is where the latest stage left the atoms
    (one stopped part-way counts: the minimiser leaves them no worse than it found them);
    :attr:`began_last` where it began."""

    def __init__(self, made, protocol, platform, properties, start, held, inbox, ident):
        self.made, self.protocol = made, protocol
        self.platform, self.properties = platform, properties
        self.start, self.held = start, held
        self.inbox, self.ident = inbox, ident
        self.stages = []
        self.first_energy = None
        self.done = start
        self.began_last = start

    def stage(self, stage, positions, scale=1, name=None):
        import openmm

        name = name or stage["name"]
        system = copy.deepcopy(self.made.system)
        for force in system.getForces():
            force.setForceGroup(0)
        if stage["moving"] == "hydrogens":
            for i in self.made.heavy:
                system.setParticleMass(i, 0.0)
        else:
            restraint = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            for parameter in ("k", "x0", "y0", "z0"):
                restraint.addPerParticleParameter(parameter)
            constants = [float(stage.get(c, 0.0)) * scale for c in RESTRAINED]
            for i, c in self.held:
                if constants[c] > 0:  # always to the deposited start, whatever stage this is
                    x, y, z = (float(v) for v in self.start[i])
                    restraint.addParticle(i, [constants[c], x, y, z])
            restraint.setForceGroup(1)
            system.addForce(restraint)
        context = openmm.Context(system, openmm.VerletIntegrator(0.001), self.platform,
                                 self.properties)
        context.setPositions(positions)
        before = _energy(context)
        if self.first_energy is None:
            self.first_energy = before
        reporter = _reporter(name, self.ident, self.inbox.stop)
        began = time.perf_counter()
        openmm.LocalEnergyMinimizer.minimize(context, self.protocol["tolerance"],
                                             self.protocol["max_iterations"], reporter)
        out = context.getState(getPositions=True).getPositions(asNumpy=True)
        out = out.value_in_unit(openmm.unit.nanometer)
        self.stages.append({"name": name, "energy_before": before,
                            "energy_after": _energy(context), "iterations": reporter.count,
                            "seconds": time.perf_counter() - began})
        del context
        self.began_last, self.done = positions, out
        return out


def _rmsd(positions, start, held):
    """The RMSD (A) of the held heavy atoms from where they started (both in nm)."""
    if not held:
        return 0.0
    at = [i for i, _c in held]
    return float(((positions[at] - start[at]) ** 2).sum(axis=1).mean() ** 0.5 * NM)


def _energy(context):
    """The potential energy without the restraints (kJ/mol)."""
    from openmm import unit

    state = context.getState(getEnergy=True, groups={0})
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def _reporter(stage, ident, stop):
    import openmm

    class Reporter(openmm.MinimizationReporter):
        count = 0

        def report(self, iteration, x, grad, args):
            self.count += 1
            if iteration % EVENT_EVERY == 0:
                _send({"id": ident, "event": {"stage": stage, "iteration": iteration,
                                              "energy": args["system energy"]}})
            return stop.is_set()

    return Reporter()


# -- the loop -----------------------------------------------------------------------------------
def main():
    inbox = Inbox()
    _send({"id": 0, "ok": True, "ready": True})
    expected = None
    while True:
        request = inbox.requests.get()
        if request is None:
            return
        ident = request.get("id")
        op = request.get("op")
        if op == "quit":
            _send({"id": ident, "ok": True})
            return
        inbox.current = ident
        inbox.stop.clear()
        if expected is None:
            expected = _expected()
        try:
            if op == "hello":
                body = hello()
            elif op == "charges":
                body = charges(request)
            elif op == "build":
                body = {"particles": build(request).system.getNumParticles()}
            elif op == "minimise":
                body = minimise(request, inbox, ident)
            else:
                _send({"id": ident, "ok": False, "error": "unknown request " + repr(op)})
                continue
            _send({"id": ident, "ok": True, **body})
        except expected as err:
            reply = {"id": ident, "ok": False, "error": f"{type(err).__name__}: {err}",
                     "traceback": traceback.format_exc(limit=8)}
            if isinstance(err, Problems):
                reply["problems"] = [{"label": label, "text": text} for label, text in err.found]
            _send(reply)
        finally:
            inbox.current = None


main()
'''


def command(python: str, wsl: bool = False, launcher: str = "wsl.exe") -> list[str]:
    """Argv that starts the worker in *python* (a Linux path run through *launcher*, the
    ``wsl.exe`` set in Settings, when *wsl*)."""
    native = [str(python), *ARGS, SCRIPT]
    return [launcher, "-e", *native] if wsl else native
