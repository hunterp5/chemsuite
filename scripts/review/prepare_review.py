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
"""Protein preparation P0: run the baselines on the neutron set and score them.

Everything lives under ``samples/prepare/`` (local data, not committed):

- ``neutron/*.cif``: the deposits (``ids.txt`` lists them; docs/prepare/p0-review.md says how
  the list was made);
- ``ccd/``: CCD components for the ligands (their ideal hydrogen counts);
- ``inputs/<ID>.pdb``: heavy atoms only, half the Asn/Gln/His flipped (``prepare_truth``);
- ``out/<method>/<ID>.pdb``: each method's output.

Steps, each skipping work already done::

    python scripts/review/prepare_review.py inputs
    python scripts/review/prepare_review.py reduce    # Reduce built in WSL (~/bench/reduce)
    python scripts/review/prepare_review.py propka    # in samples/prepare/bench-venv
    python scripts/review/prepare_review.py gemmi     # gemmi h with the CCP4 monomer library
    python scripts/review/prepare_review.py protoss   # ProteinsPlus, public PDB IDs only
    python scripts/review/prepare_review.py chemsuite # our own pipeline (--titrate: P2's states)
    python scripts/review/prepare_review.py fit       # fit its weights on half the entries
    python scripts/review/prepare_review.py score > /tmp/p0-scores.md
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
import zlib
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_truth as pt

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "prepare"
NEUTRON = DATA / "neutron"
SEED = 20261112
BENCH_PY = DATA / "bench-venv" / "Scripts" / "python.exe"
REDUCE = "~/bench/reduce/reduce_src/reduce"
REDUCE_DICT = "~/bench/reduce/reduce_wwPDB_het_dict.txt"
FITTED = None  # chemsuite's weights: None for the defaults (energy.Weights, fitted)
METHODS = ("naive", "reduce", "propka", "gemmi", "protoss", "moe", "chemsuite", "titrated")


def entries() -> list[str]:
    return sorted(line.strip() for line in (NEUTRON / "ids.txt").read_text().split() if line)


def deposit(pdb_id: str) -> Path:
    return NEUTRON / f"{pdb_id}.cif"


def crystal_ph(pdb_id: str) -> float:
    """The crystallisation pH, from its own item or the details text; 7.0 when neither says."""
    block = gemmi.cif.read(str(deposit(pdb_id))).sole_block()
    value = pt.crystal_ph(block)
    if value is None:
        details = block.find_value("_exptl_crystal_grow.pdbx_details") or ""
        match = re.search(r"p[HD]\s*=?\s*(\d+(?:\.\d+)?)", details)
        value = float(match.group(1)) if match else None
    return 7.0 if value is None else value


# -- ligands' ideal hydrogens (CCD) ------------------------------------------------------------
def ligand_names(pdb_id: str) -> set[str]:
    structure = gemmi.read_structure(str(deposit(pdb_id)))
    return {
        r.name
        for ch in structure[0]
        for r in ch
        if r.het_flag == "H" and r.name not in pt.WATERS and r.name not in pt.METALS
    }


def ideal(names: set[str]) -> dict[str, dict[str, int]]:
    """Hydrogens per heavy atom of each component, as the CCD draws it (downloaded with curl
    into ``ccd/``)."""
    from chemsuite.protein.structure.components import CcdCache

    folder = DATA / "ccd"
    folder.mkdir(parents=True, exist_ok=True)
    cache = CcdCache(folder, offline=True)
    found: dict[str, dict[str, int]] = {}
    for name in sorted(names):
        path = folder / f"{name}.cif"
        if not path.exists() and not (folder / f"{name}.missing").exists():
            url = f"https://files.rcsb.org/ligands/download/{name}.cif"
            done = subprocess.run(["curl", "-sf", url, "-o", str(path)], check=False)
            if done.returncode:
                (folder / f"{name}.missing").write_bytes(b"")
        component = cache.get(name)
        mol = component.to_mol() if component else None
        if mol is None:
            continue
        names_ = [a[0] for a in component.atoms]
        found[name] = {names_[i]: atom.GetTotalNumHs() for i, atom in enumerate(mol.GetAtoms())}
    return found


_TRUTH: dict[str, pt.Truth] = {}
_IDEAL: dict[str, dict[str, int]] = {}


def truths() -> dict[str, pt.Truth]:
    if not _TRUTH:
        ids = entries()
        _IDEAL.update(ideal(set().union(*(ligand_names(i) for i in ids))))
        for pdb_id in ids:
            _TRUTH[pdb_id] = pt.truth(deposit(pdb_id), _IDEAL)
            _TRUTH[pdb_id].ph = crystal_ph(pdb_id)
    return _TRUTH


# -- inputs ------------------------------------------------------------------------------------
def make_inputs() -> None:
    total = 0
    for pdb_id in entries():
        out = DATA / "inputs" / f"{pdb_id}.pdb"
        flipped = pt.make_input(deposit(pdb_id), out, SEED)
        listed = [[*key, name] for (key, name), flip in flipped.items() if flip]
        out.with_suffix(".flips.json").write_text(json.dumps(listed))
        total += len(listed)
    print(f"{len(entries())} inputs, {total} side chains flipped")


def flipped_sites(pdb_id: str) -> set[tuple[pt.Key, str]]:
    listed = json.loads((DATA / "inputs" / f"{pdb_id}.flips.json").read_text())
    return {((c, n, i), name) for c, n, i, name in listed}


def wsl_path(path: Path) -> str:
    drive, rest = path.drive.rstrip(":").lower(), path.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


# -- baselines ---------------------------------------------------------------------------------
def run_reduce() -> None:
    """One WSL session, four Reduce processes at a time (starting WSL per file is slow)."""
    out = DATA / "out" / "reduce"
    out.mkdir(parents=True, exist_ok=True)
    todo = [i for i in entries() if not (out / f"{i}.pdb").exists()]
    inputs, outputs = wsl_path(DATA / "inputs"), wsl_path(out)
    script = DATA / "run_reduce.sh"  # a file, so no shell on the way expands its $1
    script.write_bytes(
        (f'{REDUCE} -BUILD -DB {REDUCE_DICT} -Quiet "{inputs}/$1.pdb" > "{outputs}/$1.tmp" '
         f'2>/dev/null\nmv "{outputs}/$1.tmp" "{outputs}/$1.pdb"\necho "$1"\n').encode()
    )  # fmt: skip
    (DATA / "reduce_ids.txt").write_bytes(("\n".join(todo) + "\n").encode())
    command = f"xargs -P 4 -n 1 bash {wsl_path(script)} < {wsl_path(DATA / 'reduce_ids.txt')}"
    subprocess.run(["wsl", "-d", "Ubuntu", "--", "bash", "-lc", command], check=False)


def run_gemmi() -> None:
    out = DATA / "out" / "gemmi"
    out.mkdir(parents=True, exist_ok=True)
    for pdb_id in entries():
        target = out / f"{pdb_id}.pdb"
        if target.exists():
            continue
        done = subprocess.run(
            [str(BENCH_PY.with_name("gemmi.exe")), "h", f"--monomers={DATA / 'monomers'}", "--sort",
             str(DATA / "inputs" / f"{pdb_id}.pdb"), str(target)],
            capture_output=True, text=True, check=False,
        )  # fmt: skip
        if done.returncode:
            print(pdb_id, "failed:", done.stderr.strip().splitlines()[-1:], flush=True)


def run_chemsuite() -> None:
    """chemsuite's own hydrogens (``protein/prepare``), with the CCD entries in ``ccd/``; with
    ``--titrate``, P2's protonation states at the crystal pH (in ``out/titrated``)."""
    from chemsuite.protein.prepare.pipeline import Options, prepare
    from chemsuite.protein.prepare.sites import ALL, OFF
    from chemsuite.protein.structure import io
    from chemsuite.protein.structure.components import CcdCache

    titrate = "--titrate" in sys.argv
    out = DATA / "out" / ("titrated" if titrate else "chemsuite")
    out.mkdir(parents=True, exist_ok=True)
    ccd = CcdCache(DATA / "ccd", offline=True)
    redo = "--redo" in sys.argv
    for pdb_id in entries():
        target = out / f"{pdb_id}.pdb"
        if target.exists() and not redo:
            continue
        try:
            options = Options(ph=crystal_ph(pdb_id), weights=FITTED,
                              titrate=ALL if titrate else OFF)  # fmt: skip
            done = prepare(io.read(DATA / "inputs" / f"{pdb_id}.pdb"), options, ccd)
        except (ValueError, RuntimeError, KeyError) as error:
            print(pdb_id, "failed:", error, flush=True)
            continue
        io.write(done.structure, target)


# -- fitting the weights -----------------------------------------------------------------------
def training(pdb_id: str) -> bool:
    """Half the entries, fixed by a hash of the ID, fit the weights; the rest test them."""
    return zlib.crc32(pdb_id.encode()) % 2 == 0


def _site_answers(pdb_id: str, expected: pt.Truth):
    """The tables of one entry, and for each scored site its key and each state's answer."""
    from chemsuite.protein.prepare import energy, sites
    from chemsuite.protein.prepare.hydrogens import add_hydrogens
    from chemsuite.protein.structure import io
    from chemsuite.protein.structure.components import CcdCache

    built = add_hydrogens(io.read(DATA / "inputs" / f"{pdb_id}.pdb"),
                          ccd=CcdCache(DATA / "ccd", offline=True))  # fmt: skip
    found = sites.find_sites(built)
    s = built.structure
    table = energy.tables(s, found, crystal_ph(pdb_id))
    scored = []
    for k, site in enumerate(found):
        if site.kind == sites.ROTOR:
            continue
        r = site.residue
        name = str(s.res_name[r])
        key = (str(s.chain_name[s.chain_of[r]]), int(s.res_seq[r]), str(s.icode[r]).strip())
        where = expected.flips.get((key, name))
        first, second = pt.FLIPS[name][0]
        right = []
        for state in site.states:
            if where is None:
                right.append(None)
                continue
            at = state.xyz[0, state.names.index(first)]
            right.append(bool(np.linalg.norm(at - where[first])
                              < np.linalg.norm(at - where[second])))  # fmt: skip
        scored.append((k, name, right, [st.tautomer for st in site.states],
                       expected.states.get((key, name))))  # fmt: skip
    return table, scored


def _accuracy(cases, weights) -> dict[str, tuple[int, int]]:
    from chemsuite.protein.prepare.solve import solve

    tally: dict[str, list[int]] = {}
    for table, scored in cases:
        picked = solve(table, weights)
        for k, name, right, tautomers, state in scored:
            got = picked[k]
            if right[got] is not None:
                tally.setdefault(f"{name} orientation", [0, 0])[0] += right[got]
                tally[f"{name} orientation"][1] += 1
            if name == "HIS" and state in ("HID", "HIE", "HIP"):
                tally.setdefault("His tautomer", [0, 0])[0] += tautomers[got] == state
                tally["His tautomer"][1] += 1
    return {k: (v[0], v[1]) for k, v in tally.items()}


def _mean(found: dict[str, tuple[int, int]]) -> float:
    return float(np.mean([right / total for right, total in found.values() if total]))


def fit() -> None:
    """A coordinate search over the weights on the training half; prints both halves."""
    from dataclasses import replace

    from chemsuite.protein.prepare.energy import Weights

    found = truths()
    cases = {True: [], False: []}
    for pdb_id in entries():
        try:
            cases[training(pdb_id)].append(_site_answers(pdb_id, found[pdb_id]))
        except (ValueError, RuntimeError, KeyError) as error:
            print(pdb_id, "failed:", error, flush=True)
    weights = Weights()
    best = _mean(_accuracy(cases[True], weights))
    print("start", round(best, 4), weights, flush=True)
    for _round in range(3):
        for name, values in FIT_GRID.items():
            for value in values:
                trial = replace(weights, **{name: value})
                score = _mean(_accuracy(cases[True], trial))
                if score > best + 1e-6:
                    best, weights = score, trial
            print(name, round(best, 4), weights, flush=True)
    for half, label in ((True, "training"), (False, "test")):
        print(label, {k: f"{r / t:.1%} ({r}/{t})" for k, (r, t) in
                      sorted(_accuracy(cases[half], weights).items())})  # fmt: skip


FIT_GRID = {
    "clash": (40.0, 80.0, 160.0, 320.0, 640.0),
    "lone": (4.0, 8.0, 16.0, 32.0, 64.0),
    "water": (0.0, 0.5, 1.0, 2.0),
    "flip": (-0.2, 0.0, 0.1, 0.3),
    "hip": (1.0, 2.0, 3.0, 4.0, 6.0),
    "hid": (0.2, 0.4, 0.8, 1.2, 1.6),
    "ph": (0.0, 0.25, 0.5, 1.0),
    "coulomb": (0.0, 0.5, 1.0, 2.0, 4.0),
    "burial": (0.0, 0.5, 1.0, 2.0, 4.0),
    "metal": (0.0, 1.0, 2.0, 4.0, 8.0),
}


def run_propka() -> None:
    out = DATA / "out" / "propka"
    out.mkdir(parents=True, exist_ok=True)
    for pdb_id in entries():
        if (out / f"{pdb_id}.pka").exists():
            continue
        subprocess.run(
            [str(BENCH_PY), "-m", "propka", "--quiet", str(DATA / "inputs" / f"{pdb_id}.pdb")],
            cwd=out, capture_output=True, check=False,
        )  # fmt: skip
        print(pdb_id, "ok" if (out / f"{pdb_id}.pka").exists() else "failed", flush=True)


def run_protoss() -> None:
    """ProteinsPlus Protoss on the public PDB ID (never an uploaded file). Polls each job;
    stays under the service's 30 jobs a minute."""
    out = DATA / "out" / "protoss"
    out.mkdir(parents=True, exist_ok=True)
    jobs_file = out / "jobs.json"
    jobs: dict[str, str] = json.loads(jobs_file.read_text()) if jobs_file.exists() else {}
    for pdb_id in entries():  # submit every job first, 2.5 s apart (24 a minute)
        if (out / f"{pdb_id}.pdb").exists() or pdb_id in jobs:
            continue
        body = json.dumps({"protoss": {"pdbCode": pdb_id.lower()}})
        posted = _curl(["-X", "POST", "https://proteins.plus/api/protoss_rest",
                        "-H", "Content-Type: application/json", "-H", "Accept: application/json",
                        "-d", body])  # fmt: skip
        if posted.get("location"):
            jobs[pdb_id] = posted["location"]
            jobs_file.write_text(json.dumps(jobs, indent=1))
        else:
            print(pdb_id, "no job:", posted, flush=True)
        time.sleep(2.5)
    for _round in range(120):  # then collect them, 2.5 s between polls
        waiting = [i for i in jobs if not (out / f"{i}.pdb").exists()]
        if not waiting:
            break
        for pdb_id in waiting:
            result = _curl(["-H", "Accept: application/json", jobs[pdb_id]])
            time.sleep(2.5)
            if not result.get("protein"):
                continue
            target = out / f"{pdb_id}.pdb"
            ligands = result.get("ligands") or result.get("ligand")
            if ligands:
                subprocess.run(["curl", "-sf", ligands, "-o", str(out / f"{pdb_id}.sdf")],
                               check=False)  # fmt: skip
            subprocess.run(["curl", "-sf", result["protein"], "-o", str(target)], check=False)
            print(pdb_id, "ok", flush=True)
    for pdb_id in jobs:
        if not (out / f"{pdb_id}.pdb").exists():
            print(pdb_id, "never finished", flush=True)


def _curl(args: list[str]) -> dict:
    done = subprocess.run(["curl", "-s", *args], capture_output=True, text=True, check=False)
    try:
        return json.loads(done.stdout or "{}")
    except json.JSONDecodeError:
        return {"error": done.stdout[:200]}


# -- reading each method's answer --------------------------------------------------------------
def naive(expected: pt.Truth, flipped: set) -> pt.Observed:
    """Standard states at pH 7 (His as HIE), nothing flipped, ligands as the CCD draws them."""
    standard = {"HIS": "HIE", "ASP": "ASP", "GLU": "GLU", "LYS": "LYS", "CYS": "CYS",
                "TYR": "TYR"}  # fmt: skip
    seen = pt.Observed()
    for site in expected.states:
        seen.states[site] = standard[site[1]]
    for site in expected.flips:
        seen.flips[site] = site not in flipped
    for key, name, atom in expected.ligand_h:
        seen.ligand_h[(key, name, atom)] = _IDEAL.get(name, {}).get(atom)
    return seen


def propka_states(pdb_id: str, expected: pt.Truth) -> pt.Observed | None:
    """PROPKA's pKa values turned into states at the crystal pH. PROPKA gives no His tautomer
    or orientation, so those count as not answered."""
    path = DATA / "out" / "propka" / f"{pdb_id}.pka"
    if not path.exists():
        return None
    pkas: dict[tuple[str, int, str], float] = {}
    summary = False
    for line in path.read_text().splitlines():
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            summary = True
            continue
        if summary:
            parts = line.split()
            if len(parts) >= 4 and parts[0] in pt.TITRATABLE:
                try:
                    pkas[(parts[0], int(parts[1]), parts[2])] = float(parts[3])
                except ValueError:
                    continue
            elif line.startswith("-----") and pkas:
                break
    seen = pt.Observed()
    ph = expected.ph or 7.0
    for key, name in expected.states:
        pka = pkas.get((name, key[1], key[0]))
        if pka is None:
            seen.states[(key, name)] = None
            continue
        protonated = pka > ph
        seen.states[(key, name)] = {
            "HIS": "HIP" if protonated else "HIS?",
            "ASP": "ASH" if protonated else "ASP",
            "GLU": "GLH" if protonated else "GLU",
            "LYS": "LYS" if protonated else "LYN",
            "CYS": "CYS" if protonated else "CYM",
            "TYR": "TYR" if protonated else "TYM",
        }[name]
    return seen


def protoss(pdb_id: str, expected: pt.Truth) -> pt.Observed | None:
    """Protoss's protein from its PDB file, its ligands from the SDF (matched to the deposit's
    atoms by position). Its orientations are not scored: it starts from the deposit, whose
    orientation is the answer."""
    from rdkit import Chem

    path = DATA / "out" / "protoss" / f"{pdb_id}.pdb"
    if not path.exists() or path.stat().st_size == 0:
        return None
    seen = pt.observe(path, expected)
    seen.flips = {}
    sdf = path.with_suffix(".sdf")
    placed: list[tuple[tuple[float, float, float], int]] = []
    if sdf.exists() and sdf.stat().st_size > 0:
        for mol in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False):
            if mol is None:
                continue
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() > 1:
                    hs = sum(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
                    p = conf.GetAtomPosition(atom.GetIdx())
                    placed.append(((p.x, p.y, p.z), hs))
    positions = _ligand_positions(pdb_id)
    for site in expected.ligand_h:
        where = positions.get(site)
        best = None
        if where is not None:
            for xyz, hs in placed:
                if sum((a - b) ** 2 for a, b in zip(xyz, where, strict=True)) < 0.25:
                    best = hs
                    break
        seen.ligand_h[site] = best
    return seen


def _ligand_positions(pdb_id: str) -> dict[tuple[pt.Key, str, str], tuple[float, ...]]:
    structure = gemmi.read_structure(str(deposit(pdb_id)))
    found = {}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                site = (pt.key_of(chain, residue), residue.name, atom.name)
                found.setdefault(site, tuple(atom.pos.tolist()))
    return found


def observed(method: str, pdb_id: str, expected: pt.Truth) -> pt.Observed | None:
    if method == "protoss":
        return protoss(pdb_id, expected)
    if method == "naive":
        return naive(expected, flipped_sites(pdb_id))
    if method == "propka":
        return propka_states(pdb_id, expected)
    path = DATA / "out" / method / f"{pdb_id}.pdb"
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pt.observe(path, expected)


def method_scores(method: str) -> tuple[int, dict[str, pt.Tally]]:
    """How many entries *method* answered, and its tally per question."""
    found = truths()
    pairs = []
    for pdb_id in entries():
        seen = observed(method, pdb_id, found[pdb_id])
        if seen is not None:
            pairs.append((found[pdb_id], seen))
    return len(pairs), pt.score(pairs)


def score() -> None:
    print(f"# Scores on {len(entries())} neutron structures\n")
    rows: dict[str, dict[str, pt.Tally]] = {}
    covered: dict[str, int] = {}
    for method in METHODS:
        covered[method], rows[method] = method_scores(method)
    labels = sorted({label for tallies in rows.values() for label in tallies})
    print("| Question | " + " | ".join(f"{m} ({covered[m]})" for m in METHODS) + " |")
    print("|---|" + "---|" * len(METHODS))
    for label in labels:
        cells = [rows[m].get(label, pt.Tally()).text() for m in METHODS]
        print(f"| {label} | " + " | ".join(cells) + " |")


def moe_list() -> None:
    """``moe/entries.csv`` for the user's MOE run: each input with its crystallisation pH.
    Protonate3D's outputs go in ``out/moe/<ID>.pdb`` and are scored like the others."""
    folder = DATA / "moe"
    folder.mkdir(parents=True, exist_ok=True)
    lines = ["id,input,ph,output"]
    for pdb_id in entries():
        lines.append(f"{pdb_id},inputs/{pdb_id}.pdb,{crystal_ph(pdb_id):.1f},out/moe/{pdb_id}.pdb")
    (folder / "entries.csv").write_bytes(("\n".join(lines) + "\n").encode())
    print(f"{len(lines) - 1} entries in {folder / 'entries.csv'}")


STEPS = {"inputs": make_inputs, "moe": moe_list, "reduce": run_reduce, "gemmi": run_gemmi,
         "propka": run_propka, "protoss": run_protoss, "chemsuite": run_chemsuite, "fit": fit,
         "score": score}  # fmt: skip

if __name__ == "__main__":
    STEPS[sys.argv[1]]()
