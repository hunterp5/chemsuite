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
"""Loop modeler L0: loops cut out of high-resolution chains, rebuilt, and scored
(docs/LOOPS_PLAN.md).

Everything lives under ``samples/loops/`` (local data, not committed):

- ``ids.txt``: the picked polymer entities (``1ABC_1``), one per 30% sequence cluster;
- ``pdb/<ID>.cif``: the deposits;
- ``ref/<ID>.cif``: each usable deposit as the methods see it (first model, one location,
  polymer chains only, no hydrogens), the truth the loops are scored against;
- ``cuts.json``: the cuts, each with its chain, first residue, length, secondary structure,
  B-factors and half (train or test);
- ``inputs/<cut>.cif``: the reference with the cut's residues deleted (the sequence kept, so
  the gap is there for every method);
- ``out/<method>/<cut>.cif`` and ``.json`` (time, failure): each method's output.

Steps, each skipping work already done::

    python scripts/review/loop_bench.py pick       # RCSB search: ≤ 1.8 Å, 30% clusters
    python scripts/review/loop_bench.py fetch
    python scripts/review/loop_bench.py cut
    python scripts/review/loop_bench.py caps       # chemsuite's Prepare as it is (the floor)
    python scripts/review/loop_bench.py loops      # the loop modeler, then Prepare
    python scripts/review/loop_bench.py pdbfixer   # in the cs-loopbench env (LOOP_BENCH_PYTHON)
    python scripts/review/loop_bench.py minimize --method caps --test  # CHEMSUITE_OPENMM_PYTHON
    python scripts/review/loop_bench.py score [--test|--train] > samples/loops/l0-scores.md

Scored per loop: backbone RMSD (N, CA, C, O) to the deposit with no superposition (the anchors
are fixed: how far the rest of the chain moved is reported too), closure (the peptide bonds and
angles at and inside the loop), heavy-atom clashes between the loop and everything but itself
and its neighbours in the chain, Ramachandran outliers (with ``prepare/data/ramachandran.json``
once it exists) and time.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import time
import traceback
import urllib.request
import zlib
from dataclasses import asdict, dataclass
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "loops"
PDB, REF, INPUTS, OUT = DATA / "pdb", DATA / "ref", DATA / "inputs", DATA / "out"
CUTS = DATA / "cuts.json"
SEED = 20261124
AGENT = "chemsuite-loop-bench/0.1"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD = "https://files.rcsb.org/download/{}.cif.gz"
PICK = 500  # entities sampled from the cluster representatives
LENGTHS = (2, 3, 4, 6, 8, 10, 12)
PER_LENGTH = 80  # cuts per length (about half of them in the test half)
MARGIN = 3  # residues kept between a cut and a chain end
MIN_CHAIN = 60
# fmt: off
AMINO_ACIDS = frozenset(["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"])
# fmt: on
BACKBONE = ("N", "CA", "C", "O")
PEPTIDE, CA_C_N, C_N_CA = 1.329, 116.2, 121.7  # Å, degrees (Engh & Huber)
SEVERE, CONTACT = 2.2, 3.0  # Å, heavy-atom clashes
FIXER_SCRIPT = Path(__file__).resolve().parent / "loop_bench_fixer.py"


# -- pick and fetch ----------------------------------------------------------------------------
def _query(start: int, rows: int) -> dict:
    """Protein entities of X-ray entries at 1.8 Å or better with one protein entity and no
    nucleic acid, 80–400 residues, grouped by 30% sequence identity (best resolution first)."""

    def term(attribute: str, operator: str, value) -> dict:
        return {
            "type": "terminal",
            "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value},
        }

    nodes = [
        term("rcsb_entry_info.resolution_combined", "less_or_equal", 1.8),
        term("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        term("entity_poly.rcsb_entity_polymer_type", "exact_match", "Protein"),
        term(
            "entity_poly.rcsb_sample_sequence_length",
            "range",
            {"from": 80, "to": 400, "include_lower": True, "include_upper": True},
        ),  # fmt: skip
        term("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
        term("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
    ]
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "polymer_entity",
        "request_options": {
            "group_by": {
                "aggregation_method": "sequence_identity",
                "similarity_cutoff": 30,
                "ranking_criteria_type": {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                },
            },  # fmt: skip
            "group_by_return_type": "representatives",
            "paginate": {"start": start, "rows": rows},
        },
    }


def _post(body: dict) -> dict:
    request = urllib.request.Request(SEARCH, data=json.dumps(body).encode(),
                                     headers={"User-Agent": AGENT,
                                              "Content-Type": "application/json"})  # fmt: skip
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read())


def pick() -> None:
    """Every cluster's representative, then :data:`PICK` of them at random (seeded)."""
    found: list[str] = []
    while True:
        page = _post(_query(len(found), 1000))
        found += [hit["identifier"] for hit in page.get("result_set", [])]
        print(f"{len(found)} of {page['group_by_count']}", file=sys.stderr)
        if len(found) >= page["group_by_count"] or not page.get("result_set"):
            break
    rng = np.random.default_rng(SEED)
    chosen = sorted(rng.choice(sorted(set(found)), size=min(PICK, len(found)), replace=False))
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "ids.txt").write_text("\n".join(chosen) + "\n")
    print(f"{len(chosen)} entities picked from {len(set(found))} clusters", file=sys.stderr)


def picked() -> list[tuple[str, str]]:
    """(entry, entity id) of each picked entity."""
    out = []
    for line in (DATA / "ids.txt").read_text().split():
        entry, entity = line.strip().rsplit("_", 1)
        out.append((entry, entity))
    return out


def fetch() -> None:
    """The picked entries' mmCIFs, with curl (Windows' own, whose certificate store knows
    files.rcsb.org's issuer, where Python's may not)."""
    PDB.mkdir(parents=True, exist_ok=True)
    windows = Path(os.environ.get("SYSTEMROOT", "C:/Windows")) / "System32" / "curl.exe"
    curl = str(windows) if windows.exists() else "curl"
    for k, (entry, _entity) in enumerate(picked(), 1):
        path = PDB / f"{entry}.cif"
        if path.exists():
            continue
        packed = path.with_suffix(".cif.gz")
        done = subprocess.run([curl, "-sfL", "-A", AGENT, "-o", str(packed),
                               DOWNLOAD.format(entry)], check=False)  # fmt: skip
        if done.returncode:
            print(f"[{k}] {entry} failed ({done.returncode})", file=sys.stderr)
            continue
        path.write_bytes(gzip.decompress(packed.read_bytes()))
        packed.unlink()
        print(f"[{k}] {entry}", file=sys.stderr)


def training(entry: str) -> bool:
    """Half the entries, fixed by a hash of the ID, train (the Ramachandran grids); the rest
    test (as ``prepare_review.training``)."""
    return zlib.crc32(entry.encode()) % 2 == 0


# -- cut ---------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Cut:
    name: str  # <entry>_<chain>_<first label_seq>_<length>
    entry: str
    chain: str
    first: int  # label_seq of the first residue cut
    length: int
    ss: str  # the cut's secondary structure: coil, helix, strand or mixed
    b: float  # mean B-factor of the cut's backbone atoms
    b_rel: float  # that over the chain's backbone mean
    train: bool


def reference(entry: str) -> gemmi.Structure:
    """The deposit as every method sees it: first model, one location, polymer chains only,
    no hydrogens."""
    st = gemmi.read_structure(str(PDB / f"{entry}.cif"))
    st.setup_entities()
    while len(st) > 1:
        del st[1]
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()
    return st


def usable_chain(st: gemmi.Structure, entry: str, entity: str) -> tuple[str | None, str]:
    """The first chain of *entity* whose residues are all standard amino acids with a whole
    backbone, one per position of the sequence from the first modelled to the last, joined by
    peptide bonds (no gaps), not every residue with its backbone in two locations; or why there
    is none."""
    raw = gemmi.read_structure(str(PDB / f"{entry}.cif"))
    why = "no chain of the entity"
    for chain in st[0]:
        polymer = chain.get_polymer()
        if not len(polymer) or polymer[0].entity_id != entity:
            continue
        seqs = [r.label_seq for r in polymer]
        if len(polymer) < MIN_CHAIN:
            why = "too short"
        elif any(r.name not in AMINO_ACIDS for r in polymer):
            why = "non-standard residue"
        elif any(s is None for s in seqs) or seqs != list(range(seqs[0], seqs[0] + len(seqs))):
            why = "gap in the sequence"
        elif any(any(not r.find_atom(a, "*") for a in BACKBONE) for r in polymer):
            why = "backbone atom missing"
        elif len(_backbone_altlocs(raw, chain.name)) == len(polymer):
            why = "backbone alternate locations"
        elif any(_peptide(polymer[i], polymer[i + 1]) > 1.5 for i in range(len(polymer) - 1)):
            why = "chain break"
        else:
            return chain.name, ""
    return None, why


def _backbone_altlocs(raw: gemmi.Structure, chain: str) -> set[int]:
    """The label_seq of each residue of *chain* with a backbone atom in two locations."""
    return {residue.label_seq for residue in raw[0][chain]
            if any(atom.altloc != "\0" and atom.name in BACKBONE for atom in residue)}  # fmt: skip


def _peptide(a: gemmi.Residue, b: gemmi.Residue) -> float:
    return a.find_atom("C", "*").pos.dist(b.find_atom("N", "*").pos)


def _ss_class(codes: list[int]) -> str:
    from chemsuite.macromolecule.structure.model import COIL, HELIX, STRAND

    names = {COIL: "coil", HELIX: "helix", STRAND: "strand"}
    kinds = {names[int(c)] for c in codes}
    return kinds.pop() if len(kinds) == 1 else "mixed"


def cut() -> None:
    """Find each entry's usable chain, write its reference, and cut :data:`PER_LENGTH` loops
    per length from different entries (seeded), each into its own input file."""
    from chemsuite.macromolecule.structure import io, secondary

    REF.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True, exist_ok=True)
    usable: dict[str, tuple[str, list[int], list[float], list[int], set[int]]] = {}
    reasons: dict[str, str] = {}
    for entry, entity in picked():
        if not (PDB / f"{entry}.cif").exists():
            reasons[entry] = "not fetched"
            continue
        st = reference(entry)
        chain, why = usable_chain(st, entry, entity)
        if chain is None:
            reasons[entry] = why
            continue
        path = REF / f"{entry}.cif"
        if not path.exists():
            st.make_mmcif_document().write_file(str(path))
        s = io.read(path)
        codes = secondary.assign(s)
        ours = [r for r in range(s.n_residues) if s.chain_label(int(s.chain_of[r])) == chain]
        polymer = st[0][chain].get_polymer()
        seqs = [r.label_seq for r in polymer]
        b = [float(np.mean([r.find_atom(a, "*").b_iso for a in BACKBONE])) for r in polymer]
        wobbly = _backbone_altlocs(gemmi.read_structure(str(PDB / f"{entry}.cif")), chain)
        usable[entry] = (chain, seqs, b, [int(codes[r]) for r in ours[: len(seqs)]], wobbly)
    counts: dict[str, int] = {}
    for why in reasons.values():
        counts[why] = counts.get(why, 0) + 1
    print(f"{len(usable)} usable entries; left out: {counts}", file=sys.stderr)
    cuts: list[Cut] = []
    for length in LENGTHS:
        rng = np.random.default_rng([SEED, length])
        order = rng.permutation(sorted(usable))
        taken = 0
        for entry in order:
            chain, seqs, b, codes, wobbly = usable[str(entry)]
            # the cut and MARGIN residues each side: inside the chain, one location each
            starts = [k for k in range(MARGIN, len(seqs) - MARGIN - length + 1)
                      if not wobbly & set(seqs[k - MARGIN : k + length + MARGIN])]  # fmt: skip
            if not starts:
                continue
            at = starts[int(rng.integers(len(starts)))]
            mean = float(np.mean(b))
            segment = b[at : at + length]
            cuts.append(
                Cut(
                    f"{entry}_{chain}_{seqs[at]}_{length}",
                    str(entry),
                    chain,
                    seqs[at],
                    length,
                    _ss_class(codes[at : at + length]),
                    round(float(np.mean(segment)), 2),
                    round(float(np.mean(segment)) / mean, 3),
                    training(str(entry)),
                )
            )
            taken += 1
            if taken == PER_LENGTH:
                break
    for c in cuts:
        path = INPUTS / f"{c.name}.cif"
        if path.exists():
            continue
        st = gemmi.read_structure(str(REF / f"{c.entry}.cif"))
        chain = st[0][c.chain]
        for i in reversed(range(len(chain))):
            seq = chain[i].label_seq
            if seq is not None and c.first <= seq < c.first + c.length:
                del chain[i]
        st.make_mmcif_document().write_file(str(path))
    CUTS.write_text(json.dumps([asdict(c) for c in cuts], indent=1))
    print(f"{len(cuts)} cuts", file=sys.stderr)


def cuts(half: str = "") -> list[Cut]:
    every = [Cut(**c) for c in json.loads(CUTS.read_text())]
    if half == "test":
        return [c for c in every if not c.train]
    if half == "train":
        return [c for c in every if c.train]
    return every


# -- methods -----------------------------------------------------------------------------------
def _half() -> str:
    """``test`` or ``train`` with ``--test`` or ``--train``; "" for every cut."""
    return "test" if "--test" in sys.argv else "train" if "--train" in sys.argv else ""


def _option(name: str) -> str:
    flag = f"--{name}"
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else ""


def run_caps() -> None:
    """chemsuite's Prepare as it is today: every gap longer than one residue capped (repair
    and hydrogens only, no waters, no optimization)."""
    from chemsuite.macromolecule.prepare.pipeline import NONE, Options, prepare
    from chemsuite.macromolecule.structure import io

    _run_ours("caps", lambda s: prepare(s, Options(waters=NONE, optimise=False)), io)


def run_loops() -> None:
    """The loop modeler (``prepare/loops.py``), then Prepare as :func:`run_caps` does, until
    Prepare calls it itself (L2)."""
    from chemsuite.macromolecule.prepare import loops, repair
    from chemsuite.macromolecule.prepare.pipeline import NONE, Options, prepare
    from chemsuite.macromolecule.prepare.templates import Library
    from chemsuite.macromolecule.structure import io

    library = Library()

    def make(s):
        built, _labels = loops.fill_gaps(repair.one_location(s), library)
        return prepare(built, Options(waters=NONE, optimise=False))

    _run_ours("loops", make, io)


def _run_ours(method: str, make, io) -> None:
    folder = OUT / method
    folder.mkdir(parents=True, exist_ok=True)
    todo = [c for c in cuts() if not (folder / f"{c.name}.json").exists()]
    for k, c in enumerate(todo, 1):
        began = time.perf_counter()
        numbers: dict = {"name": c.name}
        try:
            done = make(io.read(INPUTS / f"{c.name}.cif"))
            numbers["seconds"] = time.perf_counter() - began
            io.write(done.structure, folder / f"{c.name}.cif")
            numbers["built"] = sorted(done.built)
        except Exception:  # noqa: BLE001 - a benchmark notes every failure and goes on
            numbers["failed"] = traceback.format_exc()
        (folder / f"{c.name}.json").write_text(json.dumps(numbers))
        print(f"[{k}/{len(todo)}] {c.name}", "failed" if "failed" in numbers else "",
              file=sys.stderr, flush=True)  # fmt: skip


def run_pdbfixer() -> None:
    """PDBFixer's addMissingResidues, in its own environment (``LOOP_BENCH_PYTHON``, default
    the ``cs-loopbench`` conda env) on OpenMM's CPU platform, one process for every cut."""
    env_python = Path.home() / "miniforge3" / "envs" / "cs-loopbench" / "python.exe"
    python = os.environ.get("LOOP_BENCH_PYTHON", str(env_python))
    folder = OUT / "pdbfixer"
    folder.mkdir(parents=True, exist_ok=True)
    todo = [c.name for c in cuts() if not (folder / f"{c.name}.json").exists()]
    if not todo:
        return
    env = {**os.environ, "OPENMM_DEFAULT_PLATFORM": "CPU"}
    subprocess.run([python, "-I", "-u", str(FIXER_SCRIPT), str(INPUTS), str(folder), *todo],
                   check=True, env=env)  # fmt: skip


def run_minimize() -> None:
    """Prepare (as in :func:`run_caps`) and Minimize with the app's default settings, on
    ``--method``'s outputs, the loop's heavy atoms free (``out/<method>-min``)."""
    from chemsuite.macromolecule.prepare.pipeline import NONE, Options, prepare
    from chemsuite.macromolecule.structure import io
    from chemsuite.simulate import build
    from chemsuite.simulate.config import Interpreter
    from chemsuite.simulate.minimise import minimise
    from chemsuite.simulate.runner import Simulator
    from chemsuite.simulate.spec import MinimiseSettings

    method = _option("method")
    if not method:
        sys.exit("--method NAME: whose outputs to minimize")
    python = os.environ.get("CHEMSUITE_OPENMM_PYTHON", "")
    if not python:
        sys.exit("Set CHEMSUITE_OPENMM_PYTHON to the OpenMM env's python")
    worker = Simulator(Interpreter(python))
    settings = MinimiseSettings()
    folder = OUT / f"{method}-min"
    folder.mkdir(parents=True, exist_ok=True)
    todo = [c for c in cuts(_half()) if not (folder / f"{c.name}.json").exists()
            and (OUT / method / f"{c.name}.cif").exists()]  # fmt: skip
    try:
        for k, c in enumerate(todo, 1):
            numbers: dict = {"name": c.name}
            began = time.perf_counter()
            try:
                source = io.read(OUT / method / f"{c.name}.cif")
                heavy = source.subset(source.element != 1)
                done = prepare(heavy, Options(waters=NONE, optimise=False))
                s = done.structure
                free = set(done.built) | _loop_uids(s, c)
                system = build.system_spec(s, done.states)
                relaxed = minimise(worker, s, system, settings.forcefield(), settings.protocol(),
                                   built=free, platform=settings.platform)  # fmt: skip
                io.write(s.with_coords(relaxed.result.coordinates(s)), folder / f"{c.name}.cif")
                numbers["seconds"] = time.perf_counter() - began
            except worker.error as error:
                numbers["failed"] = f"worker: {error}"
                worker.close()
                worker = Simulator(Interpreter(python))
            except Exception:  # noqa: BLE001 - a benchmark notes every failure and goes on
                numbers["failed"] = traceback.format_exc()
            (folder / f"{c.name}.json").write_text(json.dumps(numbers))
            print(f"[{k}/{len(todo)}] {c.name}", "failed" if "failed" in numbers else "",
                  file=sys.stderr, flush=True)  # fmt: skip
    finally:
        worker.close()


def _loop_uids(s, c: Cut) -> set[int]:
    """The uids of the heavy atoms of the residues between the cut's two anchors in *s*
    (found by their author numbers, as PDBFixer renumbers the sequence), which another method
    built."""
    truth = _chain_residues(REF / f"{c.entry}.cif", c.chain)
    at = [r.label_seq for r in truth].index(c.first)
    before, after = truth[at - 1].seqid, truth[at + c.length].seqid
    chain = [r for r in range(s.n_residues) if s.chain_label(int(s.chain_of[r])) == c.chain]
    key = [(int(s.res_seq[r]), str(s.icode[r]) or " ") for r in chain]
    a = key.index((before.num, before.icode or " "))
    b = key.index((after.num, after.icode or " "))
    out: set[int] = set()
    for r in chain[a + 1 : b]:
        span = s.atoms_of(r)
        out |= {int(s.uid[i]) for i in range(span.start, span.stop) if s.element[i] != 1}
    return out


# -- scoring -----------------------------------------------------------------------------------
@dataclass
class Score:
    name: str
    built: bool
    rmsd: float | None = None  # Å, loop backbone, no superposition
    drift: float | None = None  # Å, the farthest any other backbone atom of the chain moved
    peptide: float | None = None  # Å, the largest peptide-bond error at and inside the loop
    angle: float | None = None  # degrees, the largest CA–C–N or C–N–CA error there
    severe: int = 0  # heavy-atom pairs closer than SEVERE
    contacts: int = 0  # closer than CONTACT
    outliers: int = 0  # Ramachandran outliers among the loop's residues and its anchors
    rama: int = 0  # residues whose φ/ψ were judged
    seconds: float | None = None
    failed: str = ""


def _chain_residues(path: Path, chain: str) -> list[gemmi.Residue]:
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    for ch in st[0]:
        if ch.name == chain:
            return [r for r in ch if r.name in AMINO_ACIDS]
    return []


def _xyz(residue: gemmi.Residue, name: str) -> np.ndarray | None:
    atom = residue.find_atom(name, "*")
    return None if atom is None else np.array(atom.pos.tolist())


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u, w = a - b, c - b
    cos = float(np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def score_one(c: Cut, method: str, grids: dict | None) -> Score:
    numbers = json.loads((OUT / method / f"{c.name}.json").read_text())
    if "failed" in numbers:
        return Score(c.name, False, failed=numbers["failed"].strip().splitlines()[-1])
    truth = _chain_residues(REF / f"{c.entry}.cif", c.chain)
    seqs = [r.label_seq for r in truth]
    at = seqs.index(c.first)
    before, after = truth[at - 1], truth[at + c.length]
    cut_truth = truth[at : at + c.length]
    out = _chain_residues(OUT / method / f"{c.name}.cif", c.chain)
    key = [(r.seqid.num, r.seqid.icode) for r in out]
    a = key.index((before.seqid.num, before.seqid.icode))
    b = key.index((after.seqid.num, after.seqid.icode))
    loop = out[a + 1 : b]
    score = Score(c.name, False, seconds=numbers.get("seconds"))
    # how far the rest of the chain moved (matched by number)
    where = {k: r for k, r in zip(key, out, strict=True)}
    moved = []
    for r in truth[:at] + truth[at + c.length :]:
        mine = where.get((r.seqid.num, r.seqid.icode))
        if mine is not None:
            moved += [float(np.linalg.norm(_xyz(mine, n) - _xyz(r, n))) for n in BACKBONE
                      if _xyz(mine, n) is not None]  # fmt: skip
    score.drift = max(moved) if moved else None
    if len(loop) != c.length or [r.name for r in loop] != [r.name for r in cut_truth]:
        return score
    if any(_xyz(r, n) is None for r in loop for n in BACKBONE):
        return score
    score.built = True
    d = [_xyz(r, n) - _xyz(t, n) for r, t in zip(loop, cut_truth, strict=True) for n in BACKBONE]
    score.rmsd = float(np.sqrt(np.mean(np.sum(np.square(d), axis=1))))
    chain = [out[a], *loop, out[b]]
    peptide, angle = 0.0, 0.0
    for x, y in zip(chain[:-1], chain[1:], strict=True):
        cx, nx = _xyz(x, "C"), _xyz(y, "N")
        peptide = max(peptide, abs(float(np.linalg.norm(cx - nx)) - PEPTIDE))
        angle = max(angle, abs(_angle(_xyz(x, "CA"), cx, nx) - CA_C_N),
                    abs(_angle(cx, nx, _xyz(y, "CA")) - C_N_CA))  # fmt: skip
    score.peptide, score.angle = peptide, angle
    score.severe, score.contacts = _clashes(OUT / method / f"{c.name}.cif", c.chain, a, b)
    if grids is not None:
        score.outliers, score.rama = _ramachandran(out, a, b, grids)
    return score


def _clashes(path: Path, chain: str, a: int, b: int) -> tuple[int, int]:
    """Heavy-atom pairs between the loop (amino acids a+1 .. b-1 of *chain*) and any residue
    of the structure that is neither the same nor next to it in the chain (caps included)."""
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    xyz, owner = [], []
    offset = 0
    for ch in st[0]:
        k = -1
        for r in ch:
            if ch.name == chain and r.name in AMINO_ACIDS:
                k += 1
                at = k
            else:  # caps, other chains: never a neighbour of the loop
                at = 10**6 + offset
                offset += 2
            for atom in r:
                xyz.append(atom.pos.tolist())
                owner.append(at)
    xyz, owner = np.array(xyz), np.array(owner)
    in_loop = (owner > a) & (owner < b)
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(CONTACT, output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]
    keep = (in_loop[i] | in_loop[j]) & (np.abs(owner[i] - owner[j]) > 1)
    d = np.linalg.norm(xyz[i[keep]] - xyz[j[keep]], axis=1)
    return int((d < SEVERE).sum()), int(len(d))


def rama_class(residue: str, following: str | None) -> str:
    if residue == "GLY":
        return "gly"
    if residue == "PRO":
        return "pro"
    if following == "PRO":
        return "prepro"
    return "general"


def _ramachandran(out: list[gemmi.Residue], a: int, b: int, grids: dict) -> tuple[int, int]:
    """Outliers (below the class's allowed level) among the anchors and the loop's residues."""
    bad = judged = 0
    for k in range(a, b + 1):
        if k == 0 or k + 1 >= len(out):
            continue
        prev, r, nxt = out[k - 1], out[k], out[k + 1]
        phi = _dihedral(_xyz(prev, "C"), _xyz(r, "N"), _xyz(r, "CA"), _xyz(r, "C"))
        psi = _dihedral(_xyz(r, "N"), _xyz(r, "CA"), _xyz(r, "C"), _xyz(nxt, "N"))
        grid = grids[rama_class(r.name, nxt.name)]
        judged += 1
        bad += int(_grid_value(grid, phi, psi) < grid["allowed"])
    return bad, judged


def _grid_value(grid: dict, phi: float, psi: float) -> float:
    step = grid["step"]
    i = int(((phi + 180.0) % 360.0) // step)
    j = int(((psi + 180.0) % 360.0) // step)
    return grid["logp"][i][j]


def _grids() -> dict | None:
    path = ROOT / "chemsuite" / "macromolecule" / "prepare" / "data" / "ramachandran.json"
    return json.loads(path.read_text())["classes"] if path.exists() else None


def score() -> None:
    sys.stdout.reconfigure(encoding="utf-8")  # the table has Å and ≤ (cp1252 on Windows)
    half = _half()
    chosen = cuts(half)
    grids = _grids()
    methods = [m for m in ("caps", "pdbfixer", "loops", "caps-min", "pdbfixer-min", "loops-min")
               if (OUT / m).exists()]  # fmt: skip
    print(f"# Loop benchmark ({half or 'all'} cuts: {len(chosen)})\n")
    rows: dict[str, list[tuple[Cut, Score]]] = {}
    for method in methods:
        rows[method] = [(c, score_one(c, method, grids)) for c in chosen
                        if (OUT / method / f"{c.name}.json").exists()]  # fmt: skip
    print(
        "| method | length | n | built | median RMSD | mean RMSD | ≤ 2 Å | max drift | "
        "peptide err > 0.1 Å | any severe clash | Rama outliers | median s |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for method, scored in rows.items():
        for length in LENGTHS:
            these = [s for c, s in scored if c.length == length]
            if these:
                print(_row(method, str(length), these))
        print(_row(method, "all", [s for _c, s in scored]))
    print("\n## By the cut's secondary structure (all lengths)\n")
    print(
        "| method | class | n | built | median RMSD | mean RMSD | ≤ 2 Å | max drift | "
        "peptide err > 0.1 Å | any severe clash | Rama outliers | median s |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for method, scored in rows.items():
        for kind in ("coil", "helix", "strand", "mixed"):
            these = [s for c, s in scored if c.ss == kind]
            if these:
                print(_row(method, kind, these))
    failures = [(m, s) for m, scored in rows.items() for _c, s in scored if s.failed]
    if failures:
        print("\n## Failures\n")
        for method, s in failures:
            print(f"- {method} {s.name}: {s.failed}")


def _row(method: str, label: str, scores: list[Score]) -> str:
    built = [s for s in scores if s.built]
    rmsd = np.array([s.rmsd for s in built])
    drift = [s.drift for s in scores if s.drift is not None]
    seconds = [s.seconds for s in scores if s.seconds is not None]
    judged = sum(s.rama for s in built)

    def share(k: int) -> str:
        return f"{100 * k / len(built):.0f}%" if built else "–"

    return " | ".join([
        f"| {method}", label, str(len(scores)), f"{len(built)}",
        f"{np.median(rmsd):.2f}" if len(rmsd) else "–",
        f"{np.mean(rmsd):.2f}" if len(rmsd) else "–",
        share(int((rmsd <= 2.0).sum())) if len(rmsd) else "–",
        f"{max(drift):.2f}" if drift else "–",
        share(sum((s.peptide or 0) > 0.1 for s in built)),
        share(sum(s.severe > 0 for s in built)),
        f"{100 * sum(s.outliers for s in built) / judged:.1f}%" if judged else "–",
        f"{np.median(seconds):.1f} |" if seconds else "– |",
    ])  # fmt: skip


STEPS = {"pick": pick, "fetch": fetch, "cut": cut, "caps": run_caps, "loops": run_loops,
         "pdbfixer": run_pdbfixer, "minimize": run_minimize, "score": score}  # fmt: skip

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in STEPS:
        sys.exit(f"usage: loop_bench.py {{{','.join(STEPS)}}}")
    STEPS[sys.argv[1]]()
