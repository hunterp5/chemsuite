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

"""Curate the pinned OpenADMET snapshot into one table (Phase C, docs/SOM_PLAN.md).

    .venv/Scripts/python.exe scripts/cyp/curate.py

Reads ``samples/cyp/snapshot`` (see fetch.py) and writes ``samples/cyp/instances.csv``: one row
per distinct compound, with a pIC50 and its standard deviation for each isoform that measured it.
Both are local test data and neither is committed.

What it does, and why:

* **Standardizes with chemsuite's own standardizer**, cleaning up, stripping salts and
  neutralizing, but *not* canonicalizing tautomers and *not* protonating. A feature then means
  the same number the table's own descriptor columns would give for the same compound.
* **Deduplicates by InChIKey**, because the same compound under two names must not land on both
  sides of a split. Where duplicates disagree, the spread is reported rather than hidden: values
  within ``AGREEMENT`` are averaged, and a pair further apart than that is dropped for that
  isoform and counted.
* **Reports train/test overlap.** The 750-compound test set was unblinded on 2026-11-04, so it
  is a sanity check and not a blind one; a compound in both would make even that worthless.
* **Reports the value distribution per isoform.** Checked when this was written: the pIC50s are
  continuous with no pileup at the assay floor (one row at the minimum per isoform), so they are
  a regression target. If a later snapshot reports inactives as a single bound, this is where it
  shows up, and the target would need censoring before anything could be scored honestly.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, rdBase

from chemsuite.analysis.standardize import Options, standardize

REPO = Path(__file__).resolve().parents[2]
SNAPSHOT = REPO / "samples" / "cyp" / "snapshot"
INSTANCES = REPO / "samples" / "cyp" / "instances.csv"

ISOFORMS = ("1A2", "2C9", "2D6", "3A4")
VALUE = "CYP{isoform}_pIC50_direct_inhibition"
SOURCES = {
    "train": "cyp-challenge-TRAIN_inhibition.csv",
    "test": "cyp-challenge-TEST-UNBLINDED_inhibition.csv",
}

#: Two measurements of one compound this far apart (pIC50 units) are averaged; further apart,
#: neither is trusted for that isoform. A tenfold disagreement in potency is not noise.
AGREEMENT = 1.0

#: The preparation: no tautomer step and no protonation, so the structure stays the one a
#: chemist typed and a feature means what the table's own descriptor columns mean.
PREPARATION = Options(cleanup=True, strip_salts=True, neutralise=True)


@dataclass
class Compound:
    """One distinct structure and every measurement that landed on it."""

    smiles: str
    inchikey: str
    name: str
    splits: set[str] = field(default_factory=set)
    values: dict[str, list[float]] = field(default_factory=lambda: defaultdict(list))
    spreads: dict[str, list[float]] = field(default_factory=lambda: defaultdict(list))

    def settled(self, isoform: str) -> tuple[float | None, float | None, bool]:
        """The value to use for *isoform*, its standard deviation, and whether it was dropped.

        One measurement is itself; several within :data:`AGREEMENT` are averaged; a wider
        disagreement is dropped, since there is no way to tell which reading is the real one.
        """
        found = self.values.get(isoform) or []
        if not found:
            return None, None, False
        if max(found) - min(found) > AGREEMENT:
            return None, None, True
        spread = self.spreads.get(isoform) or []
        return (
            statistics.fmean(found),
            statistics.fmean(spread) if spread else None,
            False,
        )


@dataclass
class Report:
    """What curation changed, so the plan can state it rather than imply it."""

    read: int = 0
    unparsed: int = 0
    unstandardized: int = 0
    salts_stripped: int = 0
    merged: int = 0
    disagreed: dict[str, int] = field(default_factory=lambda: defaultdict(int))


def prepared(smiles: str, report: Report) -> tuple[str, str] | None:
    """*smiles* standardized, as (SMILES, InChIKey); None when it cannot be read or prepared."""
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            report.unparsed += 1
            return None
        outcome = standardize(mol, PREPARATION)
        if outcome.mol is None:
            report.unstandardized += 1
            return None
        report.salts_stripped += bool(outcome.salts)
        return Chem.MolToSmiles(outcome.mol), Chem.MolToInchiKey(outcome.mol)


def read_source(path: Path, split: str, compounds: dict[str, Compound], report: Report) -> None:
    """Fold one snapshot file into *compounds*, keyed by InChIKey."""
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            report.read += 1
            made = prepared(row["SMILES"], report)
            if made is None:
                continue
            smiles, key = made
            found = compounds.get(key)
            if found is None:
                found = compounds[key] = Compound(smiles, key, row.get("Molecule_Name", ""))
            else:
                report.merged += 1
            found.splits.add(split)
            for isoform in ISOFORMS:
                column = VALUE.format(isoform=isoform)
                text = (row.get(column) or "").strip()
                if not text:
                    continue
                found.values[isoform].append(float(text))
                deviation = (row.get(f"{column}_std") or "").strip()
                if deviation:
                    found.spreads[isoform].append(float(deviation))


def write(path: Path, compounds: list[Compound], report: Report) -> None:
    """One row per compound, blank where an isoform has no usable value."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["smiles", "inchikey", "name", "split"]
    for isoform in ISOFORMS:
        columns += [f"pIC50_{isoform}", f"std_{isoform}"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for compound in compounds:
            # a compound in both files belongs to train: it can no longer check anything
            split = "train" if "train" in compound.splits else "test"
            row = [compound.smiles, compound.inchikey, compound.name, split]
            for isoform in ISOFORMS:
                value, deviation, dropped = compound.settled(isoform)
                report.disagreed[isoform] += dropped
                row += [
                    "" if value is None else f"{value:.4f}",
                    "" if deviation is None else f"{deviation:.4f}",
                ]
            writer.writerow(row)


def describe(compounds: list[Compound], report: Report) -> None:
    """Print what the plan needs to quote: counts per isoform, overlap, and the spread."""
    both = [c for c in compounds if len(c.splits) > 1]
    print(f"read {report.read} rows -> {len(compounds)} distinct compounds")
    print(f"  unreadable {report.unparsed}, not standardized {report.unstandardized}")
    print(f"  salts stripped {report.salts_stripped}, merged onto an earlier row {report.merged}")
    print(f"  in both the train and the test file: {len(both)}")
    print()
    print("| isoform | train | test | dropped (disagreed) | min | max | mean |")
    print("|---|---|---|---|---|---|---|")
    for isoform in ISOFORMS:
        kept = {"train": [], "test": []}
        for compound in compounds:
            value, _, _ = compound.settled(isoform)
            if value is not None:
                kept["train" if "train" in compound.splits else "test"].append(value)
        every = kept["train"] + kept["test"]
        low, high = (min(every), max(every)) if every else (0.0, 0.0)
        mean = statistics.fmean(every) if every else 0.0
        print(
            f"| {isoform} | {len(kept['train'])} | {len(kept['test'])} | "
            f"{report.disagreed[isoform]} | {low:.2f} | {high:.2f} | {mean:.2f} |"
        )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--out", type=Path, default=INSTANCES)
    arguments = parser.parse_args(argv[1:])

    compounds: dict[str, Compound] = {}
    report = Report()
    for split, name in SOURCES.items():
        path = arguments.snapshot / name
        if not path.exists():
            print(f"{path} is not there; run scripts/cyp/fetch.py first", file=sys.stderr)
            return 1
        read_source(path, split, compounds, report)

    ordered = sorted(compounds.values(), key=lambda c: c.inchikey)
    write(arguments.out, ordered, report)
    describe(ordered, report)
    print(f"\nwrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
