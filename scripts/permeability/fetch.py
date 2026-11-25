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

"""Download the three permeability sources at pinned digests (docs/PERMEABILITY_PLAN.md P1).

    .venv/Scripts/python.exe scripts/permeability/fetch.py

The files land in ``samples/permeability/snapshot``, which is local test data and is not
committed — the same arrangement as ``samples/cyp/``. What *does* ship is a model fitted from
them, so each source's terms are recorded here and the licence text is downloaded beside the
data rather than summarized from memory.

**Every file is pinned by SHA-256**, as ``scripts/cyp/fetch.py`` pins OpenADMET's snapshot: a
dataset amended upstream would otherwise silently change what a model was measured against.
Rerun with ``--update-digests`` to pin new ones deliberately.

Fetching goes through :func:`chemsuite.core.download.fetch_bytes` rather than ``urllib``
directly, so this script gets the app's certificate handling (``certifi`` when it is installed,
the system store otherwise), its HTTPS-only rule, and its audit-log entry per request. The
files are small enough to read into memory; nothing about the user is sent.

The three sources, and what each one actually is:

* **MDR1-MDCK efflux ratio** — Fang et al., J. Chem. Inf. Model. 2023, 63, 3263–3274
  (doi:10.1021/acs.jcim.3c00160), Biogen's public ADME set at a pinned commit. **MIT**, and the
  licence file is fetched with it. 3,521 compounds, of which 2,642 carry an efflux ratio.
  Continuous, log10, no pile-up at either assay bound. This is the cleanest of the three and
  the one the plan falls back to if the others cannot be cleared.
* **Caco-2 log Papp** — Wang et al. 2016, as Therapeutics Data Commons serves it
  (``caco2_wang``, Harvard Dataverse). 906 compounds, log10 cm/s. TDC redistributes a paper's
  supporting information; measured values are facts rather than authorship, but the citation
  belongs in the NOTICE and does not travel without it.
* **PAMPA** — NCATS' PAMPA pH 5 screen as TDC serves it (``pampa_ncats``). **A caveat that
  changed the plan:** TDC's set is 2,035 compounds with a *binary* label, not the ~6,500-row
  continuous set ADME@NCATS describes. ADME@NCATS serves its data through a JavaScript
  application with no pinnable file, so rather than scrape it this endpoint ships as a
  classification — which chemsuite's QSAR engine supports natively, class sets and all. Work
  by US government employees, so the measurements are not under copyright.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path

from chemsuite.core.download import fetch_bytes

REPO = Path(__file__).resolve().parents[2]
DESTINATION = REPO / "samples" / "permeability" / "snapshot"

#: Biogen's public ADME set, and the commit every digest below was taken at.
BIOGEN_COMMIT = "b00df003de117ce9e5b381afd886095c5f2af2d5"
_BIOGEN = (
    f"https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/{BIOGEN_COMMIT}"
)
#: Therapeutics Data Commons serves its tables as numbered Dataverse files.
_DATAVERSE = "https://dataverse.harvard.edu/api/access/datafile"


@dataclass(frozen=True)
class Source:
    """One file to fetch: where it comes from, and what it must hash to."""

    name: str
    url: str
    sha256: str
    #: rows the file has including its header, as a second check on the digest
    rows: int
    note: str


FILES = (
    Source(
        "biogen_adme_public_3521.csv",
        f"{_BIOGEN}/ADME_public_set_3521.csv",
        "2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b",
        3522,
        "MDR1-MDCK ER, continuous (MIT)",
    ),
    Source(
        "biogen_LICENSE.txt",
        f"{_BIOGEN}/LICENSE",
        "06557ca6eeee99c162df4e5e819167c6103f35faffb51f7b0a4833905280a740",
        21,
        "the MIT licence for the file above; ships beside the model",
    ),
    Source(
        "caco2_wang.tab",
        f"{_DATAVERSE}/4259569",
        "447d9f1af487c06b080145a2361dbe45af5b27afc83e1a975dd29b9ea535d6ab",
        911,  # lines, not rows: some names hold newlines
        "Caco-2 log Papp, continuous (Wang 2016 via TDC)",
    ),
    Source(
        "pampa_ncats.tab",
        f"{_DATAVERSE}/6695858",
        "18b0ae755310cad7e4b6c3cdbe4a3c9190609995103c92536d41a4536511d854",
        2036,
        "PAMPA, BINARY (NCATS via TDC)",
    ),
)


def count_rows(data: bytes) -> int:
    """Lines in *data*, header included. Not a CSV parse: a cheap second opinion on the digest."""
    return data.count(b"\n") + (1 if data and not data.endswith(b"\n") else 0)


def download(source: Source, into: Path) -> bytes:
    """Fetch *source* and write it under *into*, returning what was fetched."""
    data = fetch_bytes(source.url, purpose="permeability training data")
    into.mkdir(parents=True, exist_ok=True)
    (into / source.name).write_bytes(data)
    return data


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--into", type=Path, default=DESTINATION, help="where to put the files")
    parser.add_argument(
        "--update-digests",
        action="store_true",
        help="print the digests of what was downloaded instead of checking them, to pin them",
    )
    arguments = parser.parse_args(argv[1:])

    print(f"three permeability sources -> {arguments.into}")
    wrong = []
    for source in FILES:
        data = download(source, arguments.into)
        found, rows = hashlib.sha256(data).hexdigest(), count_rows(data)
        if arguments.update_digests:
            print(f'        "{found}",\n        {rows},  # {source.name}')
            continue
        bad = found != source.sha256 or rows != source.rows
        if bad:
            wrong.append(
                f"  {source.name}\n    expected {source.sha256} ({source.rows} lines)"
                f"\n    found    {found} ({rows} lines)"
            )
        size = len(data) / 1e6
        print(f"  {source.name:34s} {size:5.1f} MB  {'MISMATCH' if bad else 'ok':8s} {source.note}")

    if wrong:
        print(
            "\nWhat is served now is not what these digests were taken from. That is the pin\n"
            "doing its job: find out what changed upstream before measuring anything against\n"
            "it, then rerun with --update-digests to pin the new files deliberately.\n"
            + "\n".join(wrong),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
