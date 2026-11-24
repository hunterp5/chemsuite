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

"""External pKa predictors on the neutron set's inputs, as per-site CSVs for
``prepare_states.py score --tools`` (docs/prepare/p6-research.md §6).

Runs in WSL, in its own conda env (``pka-bench``: Python 3.11, pKAI 1.2.0 with its torch,
and for KaML-CBtree the KaMLs clone at ``~/bench/KaMLs`` with pycaret, catboost and biopandas;
its bundled 32-bit ``mkdssp`` needs ``libc6:i386``).
It imports nothing from chemsuite, and nothing here reaches the app's environment.

    python scripts/review/pka_tools.py pkai [--jobs N] [--redo]   # pKAI+ (the "+" model)
    python scripts/review/pka_tools.py kaml [--jobs N] [--redo]   # KaML-CBtree, one chain at a time

Every tool gets the input PROPKA got (``samples/prepare/inputs/<ID>.pdb``), cleaned the same
way: ATOM records of the 20 standard amino acids, no hydrogens, the first conformation, and
each chain's residues numbered 1, 2, … (pKAI drops every atom with an insertion code and KaML
reads residue numbers as integers). The values are mapped back to the deposited numbering and
written to ``samples/prepare/out/pka/<tool>/<ID>.csv`` (chain, res_seq, icode, group, pka).
Neither tool sees HETATM records (ligands, metals, waters), by their own design.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "samples" / "prepare"
OUT = DATA / "out" / "pka"
KAML = Path.home() / "bench" / "KaMLs" / "KaML-CBTrees"
STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}  # fmt: skip
#: the residues each tool predicts, which are also our group names
GROUPS = {"ASP", "GLU", "HIS", "CYS", "TYR", "LYS"}
Original = tuple[int, str]  # res_seq, icode


def entries() -> list[str]:
    return sorted(w for w in (DATA / "neutron" / "ids.txt").read_text().split() if w)


def _option(name: str) -> str:
    flag = f"--{name}"
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else ""


# -- the input, cleaned and renumbered -----------------------------------------------------------
def chains(pdb_id: str) -> dict[str, list[tuple[Original, list[str]]]]:
    """Each chain's residues in file order: the deposited number and the ATOM lines."""
    out: dict[str, list[tuple[Original, list[str]]]] = {}
    for line in (DATA / "inputs" / f"{pdb_id}.pdb").read_text().splitlines():
        if not line.startswith("ATOM") or line[17:20] not in STANDARD or line[16] not in " A":
            continue
        element = line[76:78].strip() or line[12:14].strip()
        if element in ("H", "D"):
            continue
        chain, where = line[21], (int(line[22:26]), line[26].strip())
        residues = out.setdefault(chain, [])
        if not residues or residues[-1][0] != where:
            residues.append((where, []))
        residues[-1][1].append(line)
    return out


def write(path: Path, chosen: dict[str, list[tuple[Original, list[str]]]]) -> None:
    """*chosen*'s chains with residues numbered from 1, altloc and insertion code blank."""
    lines, serial = [], 0
    for chain, residues in chosen.items():
        for number, (_, atoms) in enumerate(residues, 1):
            for atom in atoms:
                serial += 1
                lines.append(f"ATOM  {serial % 100000:5d} {atom[12:16]} {atom[17:21]}{chain}"
                             f"{number:4d}    {atom[30:]}")  # fmt: skip
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def save(tool: str, pdb_id: str, rows: list[tuple[str, int, str, str, float]]) -> None:
    out = OUT / tool
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / f"{pdb_id}.tmp"
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("chain", "res_seq", "icode", "group", "pka"))
        writer.writerows(rows)
    tmp.replace(out / f"{pdb_id}.csv")


# -- pKAI+ ---------------------------------------------------------------------------------------
_MODEL = None


def _pkai_start() -> None:
    import pkai
    import torch

    global _MODEL
    sys.path.insert(0, os.path.dirname(pkai.__file__))  # its modules import each other bare
    torch.set_num_threads(1)
    _MODEL = torch.jit.load(f"{os.path.dirname(pkai.__file__)}/models/pKAI+_model.pt")


def pkai_one(pdb_id: str) -> str:
    import contextlib
    import io

    from protein import Protein

    start, found = time.perf_counter(), chains(pdb_id)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / f"{pdb_id}.pdb"
        write(path, found)
        protein = Protein(str(path))
        values = []
        if protein.tit_residues:  # none in a nucleic acid alone
            protein.apply_cutoff()
            with contextlib.redirect_stdout(io.StringIO()):  # it prints every value
                values = protein.predict_pkas(_MODEL, "cpu")
    rows = []
    for chain, number, name, pka in values:
        if name in GROUPS:
            (seq, icode), _ = found[chain][number - 1]
            rows.append((chain, seq, icode, name, float(pka)))
    save("pkai", pdb_id, rows)
    return f"{pdb_id} {len(rows)} sites {time.perf_counter() - start:.1f} s"


# -- KaML-CBtree ---------------------------------------------------------------------------------
#: its feature columns, as KaML-CBtree.py lists them
_NUM = ["com_min_d0_polar", "com_min_d1_polar", "com_min_d0_nonpolar", "com_min_d1_nonpolar",
        "com_min_d0_bb_NH", "com_min_d1_bb_NH", "com_min_d0_neg_O", "com_min_d1_neg_O",
        "com_min_d0_pos_N", "com_min_d1_pos_N", "com_min_d0_hbond_o", "com_min_d1_hbond_o",
        "com_min_d0_hbond_n", "com_min_d1_hbond_n", "com_min_d0_hbond_h", "com_min_d1_hbond_h",
        "com_min_d0_C_SG", "com_min_d1_C_SG", "com_n_hv_6", "com_n_hv_9", "com_n_hv_12",
        "com_n_hv_15", "com_npolar_5", "com_npolar_10", "com_npolar_15", "com_polar_5",
        "com_polar_10", "com_polar_15", "DA1", "DA2", "DD1", "DD2", "flexibility", "n_sc_C",
        "model_pka", "buried_ratio", "metal"]  # fmt: skip
_CAT = ["rss", "rp2ss", "rp4ss", "rm2ss", "rm4ss", "charge", "group_code"]
_KAML: dict = {}


def _kaml_start() -> None:
    """Loads the two finalized models once per worker, and runs DSSP once per chain: KaML runs it
    again for every residue of the same file, and the output does not change."""
    sys.path.insert(0, str(KAML))
    os.chdir(tempfile.gettempdir())  # pycaret writes logs.log to the current folder
    import features
    from pycaret.regression import load_model

    _KAML["acid"] = load_model(str(KAML / "models" / "catboost_acid_finalized"), verbose=False)
    _KAML["base"] = load_model(str(KAML / "models" / "catboost_base_finalized"), verbose=False)
    run, seen = features.subprocess.check_output, {}

    def cached(command, *args, **kwargs):
        if not command.startswith(features.DSSP_PATH):
            return run(command, *args, **kwargs)
        if command not in seen:
            seen.clear()
            seen[command] = run(command, *args, **kwargs)
        return seen[command]

    features.subprocess = SimpleNamespace(check_output=cached, DEVNULL=subprocess.DEVNULL)
    features.bbh_dis = _distance  # its own takes x1 - y1 for x1 - x2
    _KAML["features"] = features


def _distance(x1, y1, z1, x2, y2, z2) -> float:
    """The backbone hydrogen-bond distances (DA1, DA2, DD1, DD2) as KaML's training table has
    them. Its published ``bbh_dis`` subtracts y1 from x1, which puts these four features
    about 14 Å off the training values (checked on its own training proteins: 1KSM, 4LZT, …)."""
    return float(((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5)


def kaml_one(pdb_id: str) -> str:
    """Each chain on its own (KaML takes one chain without missing atoms), in a work folder of
    its own, since KaML writes its scratch files to the current one."""
    start, found = time.perf_counter(), chains(pdb_id)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / "tools").symlink_to(KAML / "tools")
        os.chdir(work)
        try:
            rows = _kaml_chains(pdb_id, work, found)
        finally:
            os.chdir(tempfile.gettempdir())
    save("kaml", pdb_id, rows)
    return f"{pdb_id} {len(rows)} sites {time.perf_counter() - start:.1f} s"


def _kaml_chains(pdb_id: str, work: Path, found) -> list[tuple[str, int, str, str, float]]:

    import pandas as pd

    rows = []
    for chain, residues in found.items():
        sites = [[chain, atoms[0][17:20], str(n)] for n, (_, atoms) in enumerate(residues, 1)
                 if atoms[0][17:20] in GROUPS]  # fmt: skip
        if not sites:
            continue
        path = work / f"{pdb_id}_{chain}.pdb"
        write(path, {chain: residues})
        table, kept = _features(_KAML["features"].generate([path.name, *sites]))
        if not kept:
            continue
        table = table[_NUM + _CAT]
        table[_NUM] = table[_NUM].apply(pd.to_numeric, errors="coerce")  # " " → NaN, imputed
        acid, base = _KAML["acid"].predict(table), _KAML["base"].predict(table)
        for k, row in enumerate(kept):
            _, name, number = sites[row]
            (seq, icode), _ = residues[int(number) - 1]
            pka = base[k] if name in ("HIS", "LYS") else acid[k]
            rows.append((chain, seq, icode, name, round(float(pka), 3)))
    return rows


def _features(text: str):
    """KaML's feature table, read as KaML-CBtree.py reads it, without the rows it could not
    compute (a residue missing a side-chain atom gets an error row with one field too many),
    and the indices of the rows kept. A residue left out gets no value: its form in water."""
    import io

    import pandas as pd

    lines = text.strip().splitlines()
    width = lines[0].count(",")
    kept = [k for k, line in enumerate(lines[1:]) if line.count(",") == width]
    table = pd.read_csv(io.StringIO("\n".join([lines[0]] + [lines[k + 1] for k in kept])))
    return table, kept


STEPS = {"pkai": (pkai_one, _pkai_start), "kaml": (kaml_one, _kaml_start)}


def run() -> None:
    tool = sys.argv[1]
    one, start = STEPS[tool]
    redo = "--redo" in sys.argv
    only = [w for w in _option("ids").split(",") if w]
    todo = [i for i in (only or entries()) if redo or not (OUT / tool / f"{i}.csv").exists()]
    with Pool(int(_option("jobs") or 1), initializer=start) as pool:
        for message in pool.imap_unordered(_guarded, [(one, i) for i in todo]):
            print(message, flush=True)


def _guarded(job) -> str:
    one, pdb_id = job
    try:
        return one(pdb_id)
    except (ValueError, KeyError, IndexError, RuntimeError, OSError,
            subprocess.CalledProcessError) as error:  # fmt: skip
        return f"{pdb_id} failed: {type(error).__name__}: {error}"


if __name__ == "__main__":
    run()
