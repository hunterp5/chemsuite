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
"""Protein preparation P2: experimental pKa values (PKAD-R) against chemsuite and PROPKA.

PKAD-R: Chen et al., *J. Comput. Biophys. Chem.* (2025), doi:10.1142/S2737416525500164,
the curated successor of PKAD-2 (https://compbio.clemson.edu/pkad-r/). It states no licence,
so it is used for benchmarking only and never committed: the CSV, saved by hand, lives at
``samples/prepare/pkad/pkadr.csv``, with the structures (public PDB entries) beside it.

Steps, each skipping work already done::

    python scripts/review/pkad_review.py fetch      # the PDB entries from RCSB
    python scripts/review/pkad_review.py propka     # PROPKA 3.5.1 in samples/prepare/bench-venv
    python scripts/review/pkad_review.py tables     # chemsuite's energy tables, pickled
    python scripts/review/pkad_review.py fit        # the titration weights, on the training half
    python scripts/review/pkad_review.py chemsuite  # pKa values with energy.Weights() (or FITTED)
    python scripts/review/pkad_review.py score      # RMSE per method, both halves

``tables`` and ``chemsuite`` take ``--redo``; ``chemsuite --mc`` uses Monte Carlo.
"""

from __future__ import annotations

import csv
import pickle
import subprocess
import sys
import zlib
from dataclasses import dataclass, replace
from functools import cache
from pathlib import Path

import numpy as np

from chemsuite.protein.prepare.energy import Weights

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "prepare" / "pkad"
TABLE = DATA / "pkadr.csv"
PDBS = DATA / "pdb"
BENCH_PY = ROOT / "samples" / "prepare" / "bench-venv" / "Scripts" / "python.exe"
CCD = ROOT / "samples" / "prepare" / "ccd"
TABLES = DATA / "tables"
#: PROPKA 3's model pKa values: the null model
MODEL = {"ASP": 3.8, "GLU": 4.5, "HIS": 6.5, "CYS": 9.0, "TYR": 10.0, "LYS": 10.5}


@dataclass(frozen=True)
class Entry:
    pdb: str
    chain: str
    number: int
    name: str
    pka: float
    family: str  # proteins over 30% identical share a family

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.number, self.name)


@cache
def entries() -> tuple[Entry, ...]:
    """PKAD-R's wild-type main values with no warning and a single number."""
    with TABLE.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    clean = [r for r in rows if r["pKa Classification"] == "Main" and not r["Warning"]
             and not r["Mut. Pos."].strip() and r["ResName"] in MODEL]  # fmt: skip
    families = _families(clean)
    found = []
    for r in clean:
        try:
            pka = float(r["Expt. pKa"])
            number = int(r["ResID in PDB"])
        except ValueError:
            continue
        pdb = r["PDB"].strip().upper()
        found.append(Entry(pdb, r["Chain"].strip(), number, r["ResName"], pka, families[pdb]))
    return tuple(found)


def _families(rows: list[dict]) -> dict[str, str]:
    """Each PDB ID's family: linked through the "Sequence Identity > 30%" column."""
    parent: dict[str, str] = {}

    def root(x: str) -> str:
        while parent.setdefault(x, x) != x:
            x = parent[x]
        return x

    for r in rows:
        pdb = r["PDB"].strip().upper()
        root(pdb)
        for other in r["Sequence Identity > 30%"].replace(";", ",").split(","):
            other = other.strip().split(".")[0].upper()
            if other:
                parent[root(other)] = root(pdb)
    return {pdb: root(pdb) for pdb in list(parent)}


def training(entry: Entry) -> bool:
    """Half the families, fixed by a hash of the family, fit the weights."""
    return zlib.crc32(entry.family.encode()) % 2 == 0


def structures() -> list[str]:
    return sorted({e.pdb for e in entries()})


# -- inputs and PROPKA -------------------------------------------------------------------------
def fetch() -> None:
    PDBS.mkdir(parents=True, exist_ok=True)
    for pdb in structures():
        target = PDBS / f"{pdb}.pdb"
        if target.exists():
            continue
        url = f"https://files.rcsb.org/download/{pdb}.pdb"
        done = subprocess.run(["curl", "-sf", url, "-o", str(target)], check=False)
        print(pdb, "ok" if done.returncode == 0 else "failed", flush=True)


def run_propka() -> None:
    out = DATA / "propka"
    out.mkdir(parents=True, exist_ok=True)
    for pdb in structures():
        if (out / f"{pdb}.pka").exists() or not (PDBS / f"{pdb}.pdb").exists():
            continue
        subprocess.run([str(BENCH_PY), "-m", "propka", "--quiet", str(PDBS / f"{pdb}.pdb")],
                       cwd=out, capture_output=True, check=False)  # fmt: skip
        print(pdb, "ok" if (out / f"{pdb}.pka").exists() else "failed", flush=True)


def propka_pkas(pdb: str) -> dict[tuple[str, int, str], float]:
    path = DATA / "propka" / f"{pdb}.pka"
    if not path.exists():
        return {}
    found: dict[tuple[str, int, str], float] = {}
    summary = False
    for line in path.read_text().splitlines():
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            summary = True
            continue
        if summary:
            parts = line.split()
            if len(parts) >= 4 and parts[0] in MODEL:
                try:
                    pka = float(parts[3])
                    number = int(parts[1])
                except ValueError:
                    continue
                if pka < 99:  # 99.99: a Cys in a disulfide, which nothing can titrate
                    found[(parts[2], number, parts[0])] = pka
            elif line.startswith("-----") and found:
                break
    return found


# -- chemsuite ------------------------------------------------------------------------------
def build_tables() -> None:
    """Each structure prepared (one location, missing atoms built, hydrogens), its titration
    sites found and their energy tables pickled with each site's key."""
    from chemsuite.protein.prepare import energy, hydrogens, repair, sites, templates
    from chemsuite.protein.structure import io
    from chemsuite.protein.structure.components import CcdCache

    TABLES.mkdir(parents=True, exist_ok=True)
    ccd = CcdCache(CCD)
    library = templates.Library(ccd)
    redo = "--redo" in sys.argv
    for pdb in structures():
        target = TABLES / f"{pdb}.pkl"
        if target.exists() and not redo:
            continue
        try:
            s = repair.build_missing(repair.one_location(io.read(PDBS / f"{pdb}.pdb")), library)
            built = hydrogens.add_hydrogens(s, library, ccd)
        except (ValueError, RuntimeError, KeyError) as error:
            print(pdb, "failed:", error, flush=True)
            continue
        found = sites.find_sites(built, titrate=sites.ALL)
        table = energy.tables(built.structure, found, 7.0)
        keys = [_key(built.structure, site) for site in found]
        with target.open("wb") as handle:
            pickle.dump({"keys": keys, "single": table.single, "pairs": table.pairs,
                         "protons": table.protons}, handle)  # fmt: skip
        print(pdb, len(found), "sites", flush=True)


def _key(s, site) -> tuple[str, int, str]:
    """(chain, number, name) of a site; an N-terminus is NTERM."""
    r = site.residue
    nterm = any(st.form == "NH3+" for st in site.states)
    name = "NTERM" if nterm else str(s.res_name[r])
    return (str(s.chain_name[s.chain_of[r]]), int(s.res_seq[r]), name)


@cache
def loaded() -> dict[str, tuple[list, object]]:
    """Each structure's site keys and tables."""
    from chemsuite.protein.prepare.energy import Tables

    out = {}
    for path in sorted(TABLES.glob("*.pkl")):
        with path.open("rb") as handle:
            d = pickle.load(handle)
        out[path.stem] = (d["keys"], Tables([], d["single"], d["pairs"], 7.0, d["protons"]))
    return out


def chemsuite_run(weights, method: str = "mean field", only: set[str] | None = None) -> dict:
    """pKa values of every PKAD-R entry with tables: {Entry: pKa}."""
    from chemsuite.protein.prepare import populations

    by_pdb: dict[str, list[Entry]] = {}
    for e in entries():
        by_pdb.setdefault(e.pdb, []).append(e)
    found = {}
    for pdb, (keys, table) in loaded().items():
        wanted = [e for e in by_pdb.get(pdb, []) if only is None or e.family in only]
        if not wanted:
            continue
        done = populations.titration(table, weights, method=method)
        at = {k: pka for k, pka in zip(keys, done.pka, strict=True) if pka == pka}
        for e in wanted:
            if e.key in at:
                found[e] = at[e.key]
    return found


FITTED: Weights | None = None  # the weights to try, until the defaults hold them


def chemsuite() -> None:
    method = "monte carlo" if "--mc" in sys.argv else "mean field"
    found = chemsuite_run(FITTED or Weights(), method)
    with (DATA / "chemsuite.csv").open("w", encoding="utf-8", newline="") as handle:
        out = csv.writer(handle)
        out.writerow(["pdb", "chain", "number", "name", "pka"])
        for e, pka in sorted(found.items(), key=lambda kv: (kv[0].pdb, kv[0].key)):
            out.writerow([e.pdb, e.chain, e.number, e.name, f"{pka:.2f}"])
    print(len(found), "values written")


def rmse(found: dict, half: bool | None) -> float:
    error = [p - e.pka for e, p in found.items() if half is None or training(e) == half]
    return float(np.sqrt(np.mean(np.square(error)))) if error else float("inf")


def fit() -> None:
    """The titration weights (P1's stay as they are), fitted by Nelder-Mead on the training
    half's RMSE; prints both halves."""
    from scipy.optimize import minimize

    train = {e.family for e in entries() if training(e)}
    base = FITTED or Weights()
    start = np.array([_weight(base, name) for name in FIT_NAMES])

    def weights_of(x: np.ndarray):
        return replace(base, **{n: float(abs(v)) for n, v in zip(FIT_NAMES, x, strict=True)})

    def loss(x: np.ndarray) -> float:
        score = rmse(chemsuite_run(weights_of(x), only=train), True)
        print(round(score, 4), np.round(np.abs(x), 3).tolist(), flush=True)
        return score

    step = np.maximum(np.abs(start) * 0.5, 0.2)
    simplex = np.vstack([start] + [start + step[k] * np.eye(len(start))[k]
                                   for k in range(len(start))])  # fmt: skip
    done = minimize(loss, start, method="Nelder-Mead",
                    options={"initial_simplex": simplex, "maxfev": 250, "xatol": 0.01,
                             "fatol": 1e-3})  # fmt: skip
    weights = weights_of(done.x)
    print({n: round(_weight(weights, n), 3) for n in FIT_NAMES})
    found = chemsuite_run(weights)
    print("training", round(rmse(found, True), 3), "test", round(rmse(found, False), 3))


def _weight(weights, name: str) -> float:
    """One weight by name."""
    from chemsuite.protein.prepare.energy import TERMS

    return float(weights.vector()[TERMS.index(name)])


FIT_NAMES = ("bump", "repel", "screened", "deep", "ionic", "desolv_acid", "desolv_base", "anti")


# -- scores --------------------------------------------------------------------------------
def predictions(method: str) -> dict[Entry, float]:
    if method == "null":
        return {e: MODEL[e.name] for e in entries()}
    if method == "propka":
        found = {}
        for e in entries():
            pka = propka_pkas(e.pdb).get(e.key)
            if pka is not None:
                found[e] = pka
        return found
    if method == "chemsuite":
        return chemsuite_pkas()
    raise ValueError(method)


def chemsuite_pkas() -> dict[Entry, float]:
    path = DATA / "chemsuite.csv"
    if not path.is_file():
        return {}
    by_key = {(e.pdb, *e.key): e for e in entries()}
    found = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            e = by_key.get((row["pdb"], row["chain"], int(row["number"]), row["name"]))
            if e is not None and row["pka"]:
                found[e] = float(row["pka"])
    return found


def stats(pairs: list[tuple[Entry, float]]) -> str:
    if not pairs:
        return "–"
    error = np.array([p - e.pka for e, p in pairs])
    rmse = float(np.sqrt(np.mean(error**2)))
    within = float(np.mean(np.abs(error) <= 1.0))
    return f"{rmse:.2f} ({within:.0%} ≤1, n={len(pairs)})"


def score() -> None:
    methods = ("null", "propka", "chemsuite")
    found = {m: predictions(m) for m in methods}
    # every method on the same residues: those all of them answered
    shared = [e for e in entries() if all(e in found[m] for m in methods if found[m])]
    print(f"# PKAD-R: {len(entries())} values, {len(shared)} answered by every method\n")
    print("RMSE (share within 1 unit, count)\n")
    print("| Residues | Half | " + " | ".join(methods) + " |")
    print("|---|---|" + "---|" * len(methods))
    for name in ("all", *MODEL):
        for half, label in ((True, "training"), (False, "test")):
            chosen = [e for e in shared if training(e) == half and name in ("all", e.name)]
            cells = [stats([(e, found[m][e]) for e in chosen if e in found[m]]) for m in methods]
            print(f"| {name} | {label} | " + " | ".join(cells) + " |")


STEPS = {
    "fetch": fetch,
    "propka": run_propka,
    "tables": build_tables,
    "fit": fit,
    "chemsuite": chemsuite,
    "score": score,
}

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    STEPS[sys.argv[1]]()
