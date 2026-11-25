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

"""Curate the pinned permeability snapshot into one table per endpoint (P1, PERMEABILITY_PLAN).

    .venv/Scripts/python.exe scripts/permeability/curate.py

Reads ``samples/permeability/snapshot`` (see fetch.py) and writes one CSV per endpoint beside
it: ``mdck_er.csv``, ``caco2.csv`` and ``pampa.csv``, each one row per distinct compound with
its standardized SMILES, its InChIKey and its value. All of it is local test data and none of
it is committed.

What it does, and why each part has to be there:

* **Standardizes with chemsuite's own standardizer**, cleaning up, stripping salts and
  neutralizing, but *not* canonicalizing tautomers and *not* protonating. A feature then means
  the same number the table's own descriptor columns would give for the same compound — which
  is the whole point of fitting through the app's own feature pool.
* **Deduplicates by InChIKey**, because the same compound under two names must not land on
  both sides of a scaffold split. Where duplicates disagree the spread is reported rather than
  hidden: values within :data:`Endpoint.agreement` are averaged and a wider pair is dropped for
  that endpoint and counted.
* **Reports the value distribution, and fails on a pile-up at an assay bound.** A censored
  target — everything less permeable than the assay can see reported as one floor value —
  cannot be scored honestly as a regression, and this is where that shows up. Checked when
  this was written: none of the three endpoints is censored.
* **Says what it threw away.** Every count is printed, so docs/PERMEABILITY_PLAN.md can state
  the curation rather than imply it.

**The Esaki exclusions are not applied here**, because none of these three sources needs them:
each is a single named assay with its own column, not a ChEMBL pull. They become necessary the
day the Esaki 2024 ChEMBL set is folded into the Caco-2 endpoint (doi:10.1186/s13321-024-00826-z,
which kept 1,661 entries of a far larger pull by excluding ambiguous cell lines, P-gp/BCRP
overexpression and basolateral-to-apical measurements). That is a follow-up, and the Caco-2
model's size is the reason to do it.
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
SNAPSHOT = REPO / "samples" / "permeability" / "snapshot"
CURATED = REPO / "samples" / "permeability"

#: The preparation: no tautomer step and no protonation, so the structure stays the one a
#: chemist typed and a feature means what the table's own descriptor columns mean.
PREPARATION = Options(cleanup=True, strip_salts=True, neutralise=True)

#: A regression endpoint with more than this share of its rows sitting exactly on its lowest
#: or highest value is censored, not continuous, and curation stops.
BOUND_SHARE = 0.05


@dataclass(frozen=True)
class Endpoint:
    """One endpoint: where its values are, and what counts as agreement between two of them."""

    id: str
    title: str
    source: str
    smiles_column: str
    value_column: str
    task: str  # "regression" or "classification"
    unit: str
    #: two readings of one compound this far apart are averaged; further apart, neither is
    #: trusted. For a class label the only agreement is the same label.
    agreement: float
    citation: str


ENDPOINTS = (
    Endpoint(
        "mdck_er",
        "MDR1-MDCK efflux ratio",
        "biogen_adme_public_3521.csv",
        "SMILES",
        "LOG MDR1-MDCK ER (B-A/A-B)",
        "regression",
        "log10 (B-A/A-B)",
        # half a log of efflux ratio is within what the assay reproduces; a full log is not
        0.5,
        "Fang et al., J. Chem. Inf. Model. 2023, 63, 3263-3274 (MIT)",
    ),
    Endpoint(
        "caco2",
        "Caco-2 permeability",
        "caco2_wang.tab",
        "Drug",
        "Y",
        "regression",
        "log10 cm/s",
        # inter-laboratory Caco-2 scatter is 0.3-0.5 log units, so this is the assay's own noise
        0.5,
        "Wang et al. 2016, via Therapeutics Data Commons (caco2_wang)",
    ),
    Endpoint(
        "pampa",
        "PAMPA permeability",
        "pampa_ncats.tab",
        "Drug",
        "Y",
        "classification",
        "permeable (1) or not (0)",
        0.0,
        "NCATS PAMPA pH 5, via Therapeutics Data Commons (pampa_ncats)",
    ),
)


@dataclass
class Compound:
    """One distinct structure and every reading that landed on it."""

    smiles: str
    inchikey: str
    values: list[float] = field(default_factory=list)

    def settled(self, endpoint: Endpoint) -> tuple[float | None, bool]:
        """The value to use, and whether the readings disagreed too much to use any.

        One reading is itself; several within the endpoint's agreement are averaged (a class
        label is averaged only when every reading is the same one); a wider spread is dropped,
        because there is no way to tell which reading is the real one.
        """
        if not self.values:
            return None, False
        if max(self.values) - min(self.values) > endpoint.agreement:
            return None, True
        if endpoint.task == "classification":
            return self.values[0], False
        return statistics.fmean(self.values), False


@dataclass
class Report:
    """What curation changed, so the plan can state it rather than imply it."""

    read: int = 0
    no_value: int = 0
    unparsed: int = 0
    unstandardized: int = 0
    salts_stripped: int = 0
    merged: int = 0
    disagreed: int = 0


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
        try:
            # a few of these structures have an aromatic ring the InChI layer cannot kekulize;
            # RDKit's MolSanitizeException derives from ValueError
            named = Chem.MolToSmiles(outcome.mol), Chem.MolToInchiKey(outcome.mol)
        except (ValueError, RuntimeError):
            report.unstandardized += 1
            return None
        report.salts_stripped += bool(outcome.salts)
        return named


def read_source(path: Path, endpoint: Endpoint, report: Report) -> dict[str, Compound]:
    """Every usable row of *path*, folded into one entry per InChIKey.

    The file is read with :mod:`csv` and not split on newlines: some of these names and SMILES
    hold embedded newlines, so the line count and the row count are not the same number.
    """
    delimiter = "\t" if path.suffix == ".tab" else ","
    compounds: dict[str, Compound] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter=delimiter):
            report.read += 1
            text = (row.get(endpoint.value_column) or "").strip()
            if not text:
                report.no_value += 1
                continue
            made = prepared(row.get(endpoint.smiles_column) or "", report)
            if made is None:
                continue
            smiles, key = made
            found = compounds.get(key)
            if found is None:
                found = compounds[key] = Compound(smiles, key)
            else:
                report.merged += 1
            found.values.append(float(text))
    return compounds


def distribution(values: list[float], endpoint: Endpoint) -> tuple[str, str]:
    """What the values look like, and the reason to stop if there is one."""
    if endpoint.task == "classification":
        counts: dict[float, int] = defaultdict(int)
        for value in values:
            counts[value] += 1
            share = ", ".join(
                f"{label:g}: {count} ({count / len(values):.0%})"
                for label, count in sorted(counts.items())
            )
        if set(counts) - {0.0, 1.0}:
            return share, f"{endpoint.id}: a class label that is not 0 or 1"
        return share, ""

    low, high = min(values), max(values)
    at_low = sum(1 for value in values if value == low)
    at_high = sum(1 for value in values if value == high)
    text = (
        f"min {low:.2f}  median {statistics.median(values):.2f}  max {high:.2f}  "
        f"mean {statistics.fmean(values):.2f}  at min {at_low}  at max {at_high}"
    )
    worst = max(at_low, at_high) / len(values)
    if worst > BOUND_SHARE:
        return text, (
            f"{endpoint.id}: {worst:.0%} of the values sit on one bound, so this endpoint is "
            f"censored rather than continuous and cannot be fitted as a regression as it is"
        )
    return text, ""


def write(path: Path, endpoint: Endpoint, compounds: list[Compound], report: Report) -> list[float]:
    """One row per compound, in InChIKey order; returns the values written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    values: list[float] = []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["smiles", "inchikey", endpoint.id])
        for compound in compounds:
            value, disagreed = compound.settled(endpoint)
            report.disagreed += disagreed
            if value is None:
                continue
            writer.writerow([compound.smiles, compound.inchikey, f"{value:g}"])
            values.append(value)
    return values


def curate(endpoint: Endpoint, snapshot: Path, into: Path) -> str:
    """Curate one endpoint, printing what happened; returns a reason to stop, or ""."""
    report = Report()
    compounds = read_source(snapshot / endpoint.source, endpoint, report)
    ordered = sorted(compounds.values(), key=lambda compound: compound.inchikey)
    path = into / f"{endpoint.id}.csv"
    values = write(path, endpoint, ordered, report)

    print(f"\n{endpoint.title} ({endpoint.task}, {endpoint.unit})")
    print(f"  source   {endpoint.source}")
    print(f"  cite     {endpoint.citation}")
    print(
        f"  read {report.read}, of which {report.no_value} had no value for this endpoint; "
        f"{report.unparsed} unparsed, {report.unstandardized} unstandardized, "
        f"{report.salts_stripped} salts stripped"
    )
    print(f"  {report.merged} duplicates merged, {report.disagreed} dropped as disagreeing")
    if not values:
        return f"{endpoint.id}: nothing usable was curated"
    text, stop = distribution(values, endpoint)
    print(f"  {len(values)} compounds -> {path.name}")
    print(f"  {text}")
    return stop


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT, help="where fetch.py wrote")
    parser.add_argument("--into", type=Path, default=CURATED, help="where to write the tables")
    parser.add_argument("--only", action="append", default=[], help="an endpoint id; repeatable")
    arguments = parser.parse_args(argv[1:])

    chosen = [e for e in ENDPOINTS if not arguments.only or e.id in arguments.only]
    if not chosen:
        print(f"no such endpoint; try {', '.join(e.id for e in ENDPOINTS)}", file=sys.stderr)
        return 2

    stops = [stop for e in chosen if (stop := curate(e, arguments.snapshot, arguments.into))]
    if stops:
        print("\n" + "\n".join(stops), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
