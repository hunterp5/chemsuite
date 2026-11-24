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
"""Simulation M5a: Prepare + Minimise on the P0 neutron set, measured (docs/SIMULATE_PLAN.md).

Each entry of ``prepare_review``'s set is prepared as the ``chemsuite`` variant (its crystal pH,
the shipped weights and ligand model) and minimised with the app's default settings
(:class:`MinimiseSettings`) on the OpenMM worker that ``CHEMSUITE_OPENMM_PYTHON`` names. One run
writes both ``out/chemsuite-prep/<ID>.pdb`` and ``out/chemsuite-min/<ID>.pdb`` (so before and
after come from the same preparation) and ``out/minimise/<ID>.json``: the result's numbers, the
clash score before and after, and why an entry failed. Then::

    python scripts/review/minimise_review.py run [--redo] [--only ID,ID]  # about 1.5 h
    python scripts/review/minimise_review.py score > samples/prepare/m5a-scores.md

``score`` adds what needs the deposits: the neutron marks of both outputs (``prepare_truth``)
and how far each polar hydrogen is from the deposit's on the same atom.

The clash score is our own, Probe-like: pairs of atoms (hydrogens included) more than three
bonds apart whose van der Waals spheres overlap by :data:`OVERLAP` Å or more, per 1000 atoms. A
polar hydrogen and an N or O acceptor never clash; two N/O atoms clash only closer than
:data:`HBOND_HEAVY`. Ions are left out.
"""

from __future__ import annotations

import json
import os
import sys
import time
import traceback

import gemmi
import numpy as np
import prepare_review as pr
import prepare_truth as pt
from scipy import sparse
from scipy.spatial import cKDTree

OUT = pr.DATA / "out"
PREP, MIN, NUMBERS = "chemsuite-prep", "chemsuite-min", OUT / "minimise"

OVERLAP = 0.4  # Å of overlap of van der Waals spheres that makes a clash
HBOND_HEAVY = 2.5  # Å: N/O pairs closer than this clash even if they could H-bond
#: Bondi radii (Å) by atomic number; anything else 1.8
RADII = {1: 1.1, 6: 1.7, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.8, 16: 1.8, 17: 1.75, 35: 1.85, 53: 1.98}
ACCEPTORS = {7, 8}
#: hydrogens that turn about one bond, scored apart from those the geometry fixes
ROTORS = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("CYS", "SG"), ("LYS", "NZ"),
          ("LYN", "NZ"), ("CYM", "SG")}  # fmt: skip
WATER_NAMES = pt.WATERS
#: the plan's named cases (docs/prepare/p0-review.md §2)
NAMED = ("4PDJ", "7D6G", "4G0C", "2ZYE", "5MO0", "5MOQ", "4JEC", "5T8H")
LIGAND_MOVERS = 8


# -- clash score -------------------------------------------------------------------------------
def clashes(system, xyz: np.ndarray) -> int:
    """How many pairs of *system*'s atoms at *xyz* (Å, the spec's order) clash."""
    from chemsuite.simulate import spec

    n = system.n_atoms
    element = np.array(system.element)
    kinds = np.array([r.kind for r in system.residues])
    ion = kinds[np.array(system.residue_of())] == spec.ION
    radius = np.array([RADII.get(int(e), 1.8) for e in element])
    pairs = cKDTree(xyz).query_pairs(2 * radius.max() - OVERLAP, output_type="ndarray")
    if not len(pairs):
        return 0
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    keep = (radius[i] + radius[j] - d >= OVERLAP) & ~ion[i] & ~ion[j]
    bonds = np.array([(a, b) for a, b, _order in system.bonds], dtype=np.int64).reshape(-1, 2)
    one = sparse.coo_matrix((np.ones(len(bonds)), (bonds[:, 0], bonds[:, 1])), shape=(n, n))
    step = (one + one.T + sparse.identity(n)).tocsr()
    near = (step @ step @ step).tocsr()  # within three bonds
    keep &= np.asarray(near[i, j]).ravel() == 0
    polar_h = np.zeros(n, dtype=bool)
    for a, b, _order in system.bonds:
        if element[a] == 1 and element[b] in ACCEPTORS:
            polar_h[a] = True
        elif element[b] == 1 and element[a] in ACCEPTORS:
            polar_h[b] = True
    acceptor = np.isin(element, list(ACCEPTORS))
    hbond = (polar_h[i] & acceptor[j]) | (polar_h[j] & acceptor[i])
    hbond |= acceptor[i] & acceptor[j] & (d >= HBOND_HEAVY)
    return int((keep & ~hbond).sum())


# -- the run -----------------------------------------------------------------------------------
def _worker():
    from chemsuite.simulate.config import Interpreter
    from chemsuite.simulate.runner import Simulator

    python = os.environ.get("CHEMSUITE_OPENMM_PYTHON", "")
    if not python:
        sys.exit("Set CHEMSUITE_OPENMM_PYTHON to the OpenMM env's python")
    return Simulator(Interpreter(python))


def _one(worker, pdb_id: str, ccd, ligand_model) -> dict:
    from chemsuite.macromolecule.prepare.pipeline import Options, prepare
    from chemsuite.macromolecule.structure import io
    from chemsuite.simulate import build
    from chemsuite.simulate.minimise import minimise
    from chemsuite.simulate.spec import MinimiseSettings

    numbers: dict = {"id": pdb_id}
    began = time.perf_counter()
    try:
        options = Options(ph=pr.crystal_ph(pdb_id), weights=pr.FITTED,
                          **pr.VARIANTS["chemsuite"])  # fmt: skip
        done = prepare(io.read(pr.DATA / "inputs" / f"{pdb_id}.pdb"), options, ccd,
                       ligand_model=ligand_model)  # fmt: skip
    except (ValueError, RuntimeError, KeyError) as error:
        return {**numbers, "failed": "prepare", "reason": str(error)}
    s = done.structure
    io.write(s, OUT / PREP / f"{pdb_id}.pdb")
    numbers["prepare_s"] = time.perf_counter() - began
    try:
        system = build.system_spec(s, done.states, ccd=ccd)
    except build.BuildError as error:
        return {**numbers, "failed": "build", "reason": str(error)}
    settings = MinimiseSettings()
    began = time.perf_counter()
    try:
        relaxed = minimise(worker, s, system, settings.forcefield(), settings.protocol(),
                           built=done.built, platform=settings.platform)  # fmt: skip
    except (build.BuildError, worker.error) as error:
        return {**numbers, "failed": "minimise", "reason": str(error)}
    result = relaxed.result
    io.write(s.with_coords(result.coordinates(s)), OUT / MIN / f"{pdb_id}.pdb")
    start = build.positions(s, relaxed.system)
    end = np.array(result.xyz)
    saved = result.to_json()
    for heavy in ("uid", "xyz", "charges"):
        saved.pop(heavy, None)
    return {
        **numbers,
        "minimise_s": time.perf_counter() - began,
        "atoms": relaxed.system.n_atoms,
        "ligands": len(relaxed.system.molecules),
        "built": len(done.built),
        "clashes_before": clashes(relaxed.system, start),
        "clashes_after": clashes(relaxed.system, end),
        "warnings": result.warnings(settings.rmsd_limit),
        "result": saved,
    }


def run() -> None:
    from chemsuite.macromolecule.structure.components import CcdCache

    for folder in (OUT / PREP, OUT / MIN, NUMBERS):
        folder.mkdir(parents=True, exist_ok=True)
    only = set(filter(None, pr._option("only").split(",")))
    redo = "--redo" in sys.argv
    todo = [i for i in pr.entries() if (not only or i in only)
            and (redo or not (NUMBERS / f"{i}.json").exists())]  # fmt: skip
    ccd = CcdCache(pr.DATA / "ccd", offline=True)
    ligand_model = pr._ligand_model()
    if sys.platform == "win32":  # no sleep while this process runs (ES_CONTINUOUS | SYSTEM)
        import ctypes

        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    worker = _worker()
    print("worker:", worker.hello().versions, flush=True)
    try:
        for k, pdb_id in enumerate(todo, 1):
            try:
                numbers = _one(worker, pdb_id, ccd, ligand_model)
            except worker.error as error:  # the worker died: note it and start another
                numbers = {"id": pdb_id, "failed": "worker", "reason": str(error)}
                worker.close()
                worker = _worker()
            except Exception:  # noqa: BLE001 - a benchmark notes every failure and goes on
                numbers = {"id": pdb_id, "failed": "error", "reason": traceback.format_exc()}
            (NUMBERS / f"{pdb_id}.json").write_text(json.dumps(numbers, indent=1))
            r = numbers.get("result", {})
            print(f"[{k}/{len(todo)}] {pdb_id}", numbers.get("failed", ""),
                  f"rmsd {r['rmsd']:.2f}" if r else "",
                  f"clashes {numbers['clashes_before']}->{numbers['clashes_after']}" if r else "",
                  flush=True)  # fmt: skip
    finally:
        worker.close()


# -- scoring -----------------------------------------------------------------------------------
def load() -> list[dict]:
    """Every entry's numbers that a run has written."""
    return [json.loads(p.read_text()) for p in sorted(NUMBERS.glob("*.json"))]


#: what chemsuite does not simulate yet, by design: (a part of) the problem's text, its name
UNSUPPORTED = (
    ("Glycans cannot be simulated", "glycan"),
    ("Nucleic acids cannot be simulated", "nucleic acid"),
    ("is not one of the 20 amino acids", "nonstandard residue"),
    ("is bonded to", "covalent ligand"),
    ("bond orders are a guess", "ligand with atoms missing"),
    ("forbidden element 26", "heme (iron)"),
)
#: defects: chemsuite should have prepared and minimised these
DEFECTS = (
    ("missing a hydrogen", "ligand hydrogen missing"),
    ("without an OXT or a cap", "chain end not capped"),
    ("Tyrosinate (TYM)", "TYM chosen, no template"),
    ("do not match the", "residue does not match its template"),
    ("has no MN3 template", "Mn(III) ion, no template"),
    ("Particle coordinate is NaN", "NaN in minimisation"),
    ("Topology indices", "ligand topology refused"),
)


def failure(row: dict) -> tuple[str, str]:
    """("unsupported" or "defect", why) for a failed entry: a defect when any of its
    problems is one, unsupported when one is that, and a defect when none is recognised
    (the reason can be a traceback, so lines that match neither are passed over)."""
    problems = [p for p in row["reason"].replace("; ", "\n").splitlines() if p.strip()]
    for problem in problems:
        defect = next((why for text, why in DEFECTS if text in problem), None)
        if defect is not None:
            return "defect", defect
    for problem in problems:
        unsupported = next((why for text, why in UNSUPPORTED if text in problem), None)
        if unsupported is not None:
            return "unsupported", unsupported
    return "defect", "unrecognised: " + (problems[-1][:80] if problems else "no reason")


#: residues whose polar hydrogens are the protein's (any other but water is a ligand)
PROTEIN = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
           "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "ASH",
           "GLH", "LYN", "CYM", "CYX", "TYM", "ACE", "NME"}  # fmt: skip
GROUPS = ("OH/SH rotor", "Lys NH3", "N-terminus", "water", "ligand", "fixed")
TURNED = 60.0  # degrees between parent→H and the deposit's parent→D that count as turned


def _polar_hydrogens(structure: gemmi.Structure):
    """(residue key, parent atom name, residue name, parent position, H position, how many
    H the parent carries) for each hydrogen on N, O or S."""
    for chain in structure[0]:
        for residue in chain:
            key = pt.key_of(chain, residue)
            heavy = [a for a in residue if a.element.name not in pt.HYDROGEN]
            on: dict[str, tuple[gemmi.Atom, list[gemmi.Atom]]] = {}
            for h in pt._hydrogens(residue):
                parent = min(heavy, key=lambda a, h=h: a.pos.dist(h.pos), default=None)
                if (parent is not None and parent.element.name in pt.POLAR
                        and h.pos.dist(parent.pos) <= pt.bond_limit(parent)):  # fmt: skip
                    on.setdefault(parent.name, (parent, []))[1].append(h)
            for name, (parent, hs) in on.items():
                for h in hs:
                    yield (key, name, residue.name, np.array(parent.pos.tolist()),
                           np.array(h.pos.tolist()), len(hs))  # fmt: skip


def _group(residue: str, parent: str, count: int) -> str:
    if residue in WATER_NAMES:
        return "water"
    if residue not in PROTEIN:
        return "ligand"
    if parent == "N" and count >= 3:
        return "N-terminus"
    if (residue, parent) not in ROTORS:
        return "fixed"
    return "Lys NH3" if parent == "NZ" else "OH/SH rotor"


def h_distances(pdb_id: str, method: str) -> dict[str, list[tuple[float, float]]]:
    """For each polar hydrogen of *method*'s output whose atom carries an H/D in the deposit:
    the distance (Å) to the closest of those, and the angle (degrees) between the two
    parent→H directions, by group (:data:`GROUPS`). The angle does not grow when the parent
    drifts, so it tells a turned hydrogen from one carried along with its atom."""
    found: dict[tuple, list[tuple[np.ndarray, np.ndarray]]] = {}
    truth = gemmi.read_structure(str(pr.deposit(pdb_id)))
    truth.remove_alternative_conformations()
    for key, parent, _name, at, xyz, _count in _polar_hydrogens(truth):
        found.setdefault((key, parent), []).append((at, xyz))
    out = gemmi.read_structure(str(OUT / method / f"{pdb_id}.pdb"))
    groups: dict[str, list[tuple[float, float]]] = {g: [] for g in GROUPS}
    for key, parent, name, at, xyz, count in _polar_hydrogens(out):
        there = found.get((key, parent))
        if not there:
            continue
        t_at, t_xyz = min(there, key=lambda t, xyz=xyz: float(np.linalg.norm(xyz - t[1])))
        u, v = xyz - at, t_xyz - t_at
        cos = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
        angle = float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        groups[_group(name, parent, count)].append((float(np.linalg.norm(xyz - t_xyz)), angle))
    return groups


def marks(method: str, ids: list[str]) -> dict[str, pt.Tally]:
    """*method*'s neutron marks on the entries *ids* (both outputs on the same entries)."""
    found = pr.truths()
    pairs = []
    for pdb_id in ids:
        seen = pr.observed(method, pdb_id, found[pdb_id])
        if seen is not None:
            pairs.append((found[pdb_id], seen))
    return pt.score(pairs)


def _line(values) -> str:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    if not len(v):
        return "–"
    p95 = np.percentile(v, 95)
    return f"median {np.median(v):.2f}, p95 {p95:.2f}, max {v.max():.2f} (n {len(v)})"


def score() -> None:
    rows = load()
    done = [r for r in rows if "result" in r]
    failed = [r for r in rows if "failed" in r]
    print(f"# M5a: Prepare + Minimise on {len(rows)} neutron entries\n")
    print(f"Minimised {len(done)}; failed {len(failed)}:\n")
    for kind in ("defect", "unsupported"):
        found: dict[str, list[str]] = {}
        for r in failed:
            k, why = failure(r)
            if k == kind:
                found.setdefault(why, []).append(r["id"])
        print(f"{kind.capitalize()} ({sum(map(len, found.values()))}):\n")
        for why, ids in sorted(found.items(), key=lambda x: -len(x[1])):
            print(f"- {why} ({len(ids)}): {', '.join(ids)}")
        print()
    res = [r["result"] for r in done]
    print("\n## Heavy atoms held (Å, from the prepared start)\n")
    print("- overall:", _line(x["rmsd"] for x in res))
    print("- pocket:", _line(x.get("rmsd_pocket") for x in res))
    print("- ligand:", _line(x.get("rmsd_ligand") for x in res))
    over = [r["id"] for r in done if r["result"]["rmsd"] > 0.3]
    print(f"- over the 0.3 Å limit: {len(over)} {over[:20]}")
    print(f"- with a ⚠: {sum(bool(r['warnings']) for r in done)}")
    print("\n## Clash score (per 1000 atoms)\n")
    before = [1000 * r["clashes_before"] / r["atoms"] for r in done]
    after = [1000 * r["clashes_after"] / r["atoms"] for r in done]
    print("- before:", _line(before))
    print("- after:", _line(after))
    print(f"- worse after: {sum(a > b for a, b in zip(after, before, strict=True))}")
    print("\n## Time (s)\n")
    print("- prepare:", _line(r.get("prepare_s") for r in done))
    print("- minimise:", _line(r.get("minimise_s") for r in done))
    print("\n## Ligands that moved most\n")
    moved = sorted((r for r in done if r["result"].get("rmsd_ligand") is not None),
                   key=lambda r: -r["result"]["rmsd_ligand"])  # fmt: skip
    for r in moved[:LIGAND_MOVERS]:
        top = ", ".join(f"{label} {d:.2f}" for label, d in r["result"]["moved_most"][:3])
        print(f"- {r['id']}: ligand {r['result']['rmsd_ligand']:.2f} Å; most: {top}")
    print("\n## Named cases\n")
    by_id = {r["id"]: r for r in rows}
    for pdb_id in NAMED:
        r = by_id.get(pdb_id)
        if r is None:
            print(f"- {pdb_id}: not run")
        elif "failed" in r:
            print(
                f"- {pdb_id}: failed ({r['failed']}): {r['reason'].strip().splitlines()[-1][:160]}"
            )
        else:
            x = r["result"]
            print(f"- {pdb_id}: RMSD {x['rmsd']:.2f}, pocket {x.get('rmsd_pocket') or 0:.2f}, "
                  f"ligand {x.get('rmsd_ligand') or 0:.2f} Å; clashes "
                  f"{r['clashes_before']} → {r['clashes_after']}")  # fmt: skip
    print("\n## Neutron marks\n")
    print(f"| Question | {PREP} | {MIN} |\n|---|---|---|")
    tallies = {m: marks(m, [r["id"] for r in done]) for m in (PREP, MIN)}
    for label in sorted(set(tallies[PREP]) | set(tallies[MIN])):
        cells = [tallies[m].get(label, pt.Tally()).text() for m in (PREP, MIN)]
        print(f"| {label} | " + " | ".join(cells) + " |")
    print("\n## Polar hydrogens to the deposit's (the same entries)\n")
    print(
        f"| Group | Output | Distance (Å) | Within 0.5 Å | Beyond 1.0 Å | Turned > {TURNED:.0f}° |"
    )
    print("|---|---|---|---|---|---|")
    for group, method, d, a in h_table(done):
        print(f"| {group} | {method} | {_line(d)} | {np.mean(d <= 0.5):.1%} "
              f"| {np.mean(d > 1.0):.1%} | {np.mean(a > TURNED):.1%} |")  # fmt: skip


def h_table(done: list[dict]):
    """(group, output, distances, angles) over the entries *done*, both outputs."""
    for method in (PREP, MIN):
        found: dict[str, list[tuple[float, float]]] = {g: [] for g in GROUPS}
        for r in done:
            for group, values in h_distances(r["id"], method).items():
                found[group] += values
        for group, values in found.items():
            if values:
                pairs = np.array(values)
                yield group, method, pairs[:, 0], pairs[:, 1]


STEPS = {"run": run, "score": score}

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    STEPS[sys.argv[1]]()
