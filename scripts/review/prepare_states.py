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

"""Protonation states on the neutron set: does titration put each residue in the state the
neutron data show? (docs/prepare/p5-review.md)

Each entry's energy tables are built once with every group titrating (``titrate=ALL``) and
once without (P1), at the crystal's pH, and pickled with each site's observed state. Scoring
a set of weights then needs no structure work: each table is solved as the pipeline does it
(the best states, then each titrating group in its majority form at that pH), in seconds.

    python scripts/review/prepare_states.py build [--jobs N] [--redo]
    python scripts/review/prepare_states.py score [--test] [--fast]  # null, P1, titration
    python scripts/review/prepare_states.py fit     # the titration weights, training half

Per residue type it counts the share of states right, the unusual states found (a protonated
Asp or Glu, a neutral Lys, a Tyr or Cys anion, a charged His) and the wrong switches (a usual
state called unusual). His also has its charge, tautomer and orientation.
"""

from __future__ import annotations

import pickle
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_review as pr

from chemsuite.macromolecule.prepare.pipeline import SWITCH

CASES = pr.DATA / "state_cases.pkl"
TYPES = ("ASP", "GLU", "LYS", "TYR", "CYS")
#: the usual form of each type at pH 7
USUAL = {"ASP": "ASP", "GLU": "GLU", "LYS": "LYS", "TYR": "TYR", "CYS": "CYS"}


@dataclass
class Scored:
    """One site of a table: where it is, what the neutron data show, and each state's form
    and (for His) whether its ring points the deposited way."""

    index: int
    name: str
    truth: str
    forms: list[str]
    right: list[bool | None] = field(default_factory=list)


@dataclass
class Case:
    pdb_id: str
    ph: float
    titrated: object  # energy.Tables with titrate=ALL
    p1: object  # energy.Tables with titrate=OFF
    sites: list[Scored]  # indices into titrated
    p1_sites: list[Scored]  # His only, indices into p1


# pickled under this module's name, so the cases load however the script is started
Scored.__module__ = Case.__module__ = "prepare_states"
sys.modules.setdefault("prepare_states", sys.modules[__name__])


def _scored(s, found, expected, kinds) -> list[Scored]:
    from chemsuite.macromolecule.prepare import sites

    out = []
    for k, site in enumerate(found):
        r = site.residue
        name = str(s.res_name[r])
        if site.kind == sites.ROTOR and not (site.titrating and name in TYPES):
            continue
        if site.kind not in kinds:
            continue
        key = (str(s.chain_name[s.chain_of[r]]), int(s.res_seq[r]), str(s.icode[r]).strip())
        truth = expected.states.get((key, name))
        if truth is None:
            continue
        forms = [st.form or st.tautomer for st in site.states]
        if _nterm(forms):
            continue
        right: list[bool | None] = []
        where = expected.flips.get((key, name)) if name == "HIS" else None
        for st in site.states:
            if where is None:
                right.append(None)
                continue
            at = st.xyz[0, st.names.index("ND1")]
            right.append(bool(np.linalg.norm(at - where["ND1"])
                              < np.linalg.norm(at - where["CD2"])))  # fmt: skip
        out.append(Scored(k, name, truth, forms, right))
    return out


def _nterm(forms: list[str]) -> bool:
    """An N-terminus (its backbone NH3+), which shares its residue's name and number with the
    side chain the neutron data score."""
    return "NH3+" in forms


def _case(pdb_id: str) -> Case | str:
    from chemsuite.macromolecule.prepare import energy, sites
    from chemsuite.macromolecule.prepare.hydrogens import add_hydrogens
    from chemsuite.macromolecule.structure import io
    from chemsuite.macromolecule.structure.components import CcdCache

    try:
        expected = pr.truths()[pdb_id]
        ph = pr.crystal_ph(pdb_id)
        built = add_hydrogens(io.read(pr.DATA / "inputs" / f"{pdb_id}.pdb"),
                              ccd=CcdCache(pr.DATA / "ccd", offline=True))  # fmt: skip
        s = built.structure
        everything = sites.find_sites(built, titrate=sites.ALL)
        p1 = sites.find_sites(built)
        kinds = {sites.HIS, sites.ACID, sites.ROTOR}
        return Case(pdb_id, ph, energy.tables(s, everything, ph), energy.tables(s, p1, ph),
                    _scored(s, everything, expected, kinds),
                    _scored(s, p1, expected, {sites.HIS}))  # fmt: skip
    except (ValueError, RuntimeError, KeyError) as error:
        return f"{pdb_id} failed: {error}"


def build() -> None:
    from concurrent.futures import ProcessPoolExecutor

    from chemsuite.macromolecule.prepare.energy import TERMS

    if CASES.exists() and "--redo" not in sys.argv:
        print(CASES, "exists (--redo to rebuild)")
        return
    ids = pr.entries()
    jobs = int(pr._option("jobs") or 1)
    cases = []
    with ProcessPoolExecutor(jobs) as pool:
        for done in pool.map(_case, ids):
            if isinstance(done, str):
                print(done, flush=True)
            else:
                cases.append(done)
    with CASES.open("wb") as handle:
        pickle.dump({"terms": TERMS, "cases": cases}, handle)
    print(len(cases), "cases written to", CASES)


def cases() -> list[Case]:
    from chemsuite.macromolecule.prepare.energy import TERMS

    with CASES.open("rb") as handle:
        d = pickle.load(handle)
    if tuple(d["terms"]) != tuple(TERMS):
        raise SystemExit("the cases were built with other energy terms: build --redo")
    return d["cases"]


# -- deciding, as the pipeline does --------------------------------------------------------------
def _levels(table, populations_at, ph: float, switch) -> dict[int, np.ndarray]:
    """Each titrating site's allowed states: its proton count chosen as ``pipeline`` chooses it
    (:func:`~chemsuite.macromolecule.prepare.pipeline.protonated`, with the thresholds
    *switch*)."""
    from chemsuite.macromolecule.prepare import sites
    from chemsuite.macromolecule.prepare.pipeline import protonated

    allowed = {}
    for i, (n, pop) in enumerate(zip(table.protons, populations_at, strict=True)):
        if len(n) and n.max() > n.min():
            share = float(pop[n == n.max()].sum())
            name = sites.group(table.sites[i])
            up = protonated(name, share, ph, switch) if name else share >= 0.5
            allowed[i] = n == (n.max() if up else n.min())
    return allowed


def decide(table, weights, ph: float, switch) -> list[int]:
    """The best states, then each titrating group in its form at *ph* (``pipeline.prepare``
    without the report)."""
    from chemsuite.macromolecule.prepare import populations, solve

    picked = solve.solve(table, weights)
    if not any(n.max() > n.min() for n in table.protons if len(n)):
        return picked
    done = populations.titration(table, weights, grid=np.array([ph]))
    allowed = _levels(table, [pop[0] for pop in done.populations], ph, switch)
    if all(mask[picked[i]] for i, mask in allowed.items()):
        return picked
    return solve.solve(table, weights, allowed=allowed)


def decide_fast(table, weights, ph: float, switch) -> list[int]:
    """:func:`decide` read off the mean-field populations at *ph*: each titrating group's form,
    then its most populated state in that form (each other site's most populated state). Close
    to the exact answer and a hundred times faster: for the fit, whose result is then checked
    with :func:`decide`."""
    from chemsuite.macromolecule.prepare import populations

    done = populations.titration(table, weights, grid=np.array([ph]))
    at = [pop[0] for pop in done.populations]
    allowed = _levels(table, at, ph, switch)
    return [int(np.argmax(np.where(allowed[i], p, -1.0) if i in allowed else p))
            for i, p in enumerate(at)]  # fmt: skip


@dataclass
class Counts:
    right: int = 0
    total: int = 0

    def add(self, ok: bool) -> None:
        self.total += 1
        self.right += bool(ok)

    @property
    def share(self) -> float:
        return self.right / self.total if self.total else float("nan")

    def text(self) -> str:
        return f"{100 * self.share:.1f}% ({self.right}/{self.total})" if self.total else "–"


def _tally(out: dict[str, Counts], name: str, truth: str, got: str, right: bool | None) -> None:
    def add(label: str, ok: bool) -> None:
        out.setdefault(label, Counts()).add(ok)

    if name == "HIS":
        add("His charge", (got == "HIP") == (truth == "HIP"))
        add("His tautomer", got == truth)
        if truth == "HIP":
            add("His charged found", got == "HIP")
        else:
            add("His wrongly charged", got == "HIP")
        if right is not None:
            add("His orientation", right)
        return
    add(f"{name} state", got == truth)
    if truth != USUAL[name]:
        add(f"{name} unusual found", got == truth)
    else:
        add(f"{name} wrong switches", got != truth)


def evaluate(chosen: list[Case], weights, method: str, switch=None) -> dict[str, Counts]:
    """Counts per question: ``null`` (every group in its usual state, every His HIE as
    built), ``p1`` (titration off), ``titrated`` or ``fast`` (titrated, read off the
    populations: :func:`decide_fast`)."""
    out: dict[str, Counts] = {}
    for case in chosen:
        scored = [sc for sc in case.sites if not _nterm(sc.forms)]  # cases built before the fix
        if method == "null":
            for sc in scored:
                got = "HIE" if sc.name == "HIS" else USUAL[sc.name]
                right = sc.right[sc.forms.index("HIE")] if sc.name == "HIS" else None
                _tally(out, sc.name, sc.truth, got, right)
            continue
        if method == "p1":
            picked = pr_solve(case.p1, weights)
            for sc in case.p1_sites:
                k = picked[sc.index]
                _tally(out, sc.name, sc.truth, sc.forms[k], sc.right[k])
            for sc in scored:  # P1 keeps every other group in its usual state
                if sc.name != "HIS":
                    _tally(out, sc.name, sc.truth, USUAL[sc.name], None)
            continue
        how = decide_fast if method == "fast" else decide
        picked = how(case.titrated, weights, case.ph, SWITCH if switch is None else switch)
        for sc in scored:
            k = picked[sc.index]
            _tally(out, sc.name, sc.truth, sc.forms[k], sc.right[k])
    return out


def pr_solve(table, weights) -> list[int]:
    from chemsuite.macromolecule.prepare import solve

    return solve.solve(table, weights)


ORDER = ("ASP state", "ASP unusual found", "ASP wrong switches", "GLU state", "GLU unusual found",
         "GLU wrong switches", "LYS state", "LYS unusual found", "LYS wrong switches",
         "TYR state", "TYR unusual found", "TYR wrong switches", "CYS state",
         "CYS unusual found", "CYS wrong switches", "His charge", "His charged found",
         "His wrongly charged", "His tautomer", "His orientation")  # fmt: skip


def table_of(results: dict[str, dict[str, Counts]]) -> str:
    methods = list(results)
    lines = ["| Question | " + " | ".join(methods) + " |", "|---|" + "---|" * len(methods)]
    for q in ORDER:
        cells = [results[m].get(q, Counts()).text() for m in methods]
        lines.append(f"| {q} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def objective(found: dict[str, Counts]) -> float:
    """What the fit maximises: the mean share right over each type's states, His charge, His
    tautomer and His orientation. Every wrong switch costs as much as a find gains."""
    keys = [f"{t} state" for t in TYPES] + ["His charge", "His tautomer", "His orientation"]
    return float(np.mean([found[k].share for k in keys if k in found and found[k].total]))


def score() -> None:
    from chemsuite.macromolecule.prepare.energy import Weights

    weights = FITTED or Weights()
    every = cases()
    half = "--test" in sys.argv
    chosen = [c for c in every if pr.training(c.pdb_id) != half]
    label = "test" if half else "training"
    methods = ("null", "p1", "fast") if "--fast" in sys.argv else ("null", "p1", "titrated", "fast")
    results = {m: evaluate(chosen, weights, m, FITTED_SWITCH) for m in methods}
    p1_default = evaluate(chosen, Weights(), "p1")
    assert {k: v.right for k, v in p1_default.items()} == {
        k: v.right for k, v in results["p1"].items()
    }, "a fitted weight changes titration-off results: it must not be shared with P1"
    print(f"## Neutron set, {label} half ({len(chosen)} entries)\n")
    print(table_of(results))
    print("\nobjective:", {m: round(objective(r), 4) for m, r in results.items()})


FITTED = None  # weights to score, until the defaults hold them
FITTED_SWITCH: dict[str, float] | None = None  # thresholds, until pipeline.SWITCH holds them
#: the weights and thresholds to search ("switch:ASP" is ASP's threshold)
_SHARES = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
FIT_GRID: dict[str, tuple[float, ...]] = {
    **{f"switch:{g}": _SHARES for g in ("HIS", "GLU", "ASP", "LYS", "TYR", "CYS")},
    "screened": (0.0, 0.05, 0.1, 0.2, 0.4),
    "desolv_acid": (0.0, 0.5, 1.0, 2.0),
    "desolv_base": (0.0, 2.0, 4.0, 8.0),
    "bump": (4.0, 8.0, 16.0, 32.0),
    "repel": (0.0, 2.0, 4.0, 8.0),
}  # what p5-review.md searched (ionic stays 0: it also scores P1's His)


def _with(weights, switch: dict[str, float], name: str, value: float):
    if name.startswith("switch:"):
        return weights, {**switch, name.split(":")[1]: value}
    return replace(weights, **{name: value}), switch


def _value(weights, switch: dict[str, float], name: str) -> float:
    if name.startswith("switch:"):
        return switch.get(name.split(":")[1], 0.5)
    return _weight(weights, name)


def fit() -> None:
    """A coordinate search over :data:`FIT_GRID` on the training half's objective."""
    from chemsuite.macromolecule.prepare.energy import Weights

    train = [c for c in cases() if pr.training(c.pdb_id)]
    weights = FITTED or Weights()
    switch = dict(FITTED_SWITCH or SWITCH)
    best = objective(evaluate(train, weights, "fast", switch))
    print("start", round(best, 4), flush=True)
    for _round in range(3):
        moved = False
        for name, values in FIT_GRID.items():
            for value in values:
                w, sw = _with(weights, switch, name, value)
                got = objective(evaluate(train, w, "fast", sw))
                if got > best + 1e-6:
                    best, weights, switch, moved = got, w, sw, True
            print(name, round(best, 4), _value(weights, switch, name), flush=True)
        if not moved:
            break
    print("weights", {n: _value(weights, switch, n) for n in FIT_GRID if ":" not in n})
    print("switch", switch)
    print(table_of({"fast": evaluate(train, weights, "fast", switch)}))


def _weight(weights, name: str) -> float:
    from chemsuite.macromolecule.prepare.energy import TERMS

    return float(weights.vector()[TERMS.index(name)])


STEPS = {"build": build, "score": score, "fit": fit}

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    STEPS[sys.argv[1]]()
