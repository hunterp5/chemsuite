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

"""Download OpenADMET's CYP direct-inhibition data at a pinned commit (Phase C, docs/SOM_PLAN.md).

    .venv/Scripts/python.exe scripts/cyp/fetch.py

The measurements land in ``samples/cyp/``, which is local test data and is not committed — the
same arrangement as ``samples/som/``. Nothing downloaded here ever ships.

**Pinned, because the set was amended during the challenge.** The organizers refreshed the
training and test data with the finalised assay results on 2026-08-06 and unblinded the test
set on 2026-11-04, so "the current version" names different numbers depending on the day. This
script pins one commit and the SHA-256 of each file it reads, and refuses a file whose digest
does not match rather than quietly measuring something else.

**The licence is not settled.** The dataset declares ``apache-2.0`` in its Hugging Face metadata
and nowhere else: the card carries no licence text, no copyright holder, and no statement about
the measurements themselves. Training on it locally is fine. Shipping anything fitted from it is
a question to answer first — the standard CLAUDE.md already sets for pKa sources.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DESTINATION = REPO / "samples" / "cyp" / "snapshot"

#: The dataset, and the commit every digest below was taken at (2026-11-05, after the
#: 2026-08-06 "finalised assay results" refresh).
DATASET = "openadmet/cyp-challenge-train-test"
COMMIT = "5ab9dd9daee4dd98840ad585455a2b5aa386f5d0"

#: file -> (SHA-256, rows the dataset card claims). Only the direct-inhibition endpoint: the
#: TDI, Emax and single-concentration sets are different endpoints and are out of scope.
FILES = {
    "cyp-challenge-TRAIN_inhibition.csv": (
        "b8f79addd266fb6f9f4c222c5e4e73d926362328b6a8d2841871a54e46bd2278",
        4905,
    ),
    "cyp-challenge-TEST-UNBLINDED_inhibition.csv": (
        "6551eba506c90cf51434bc1764424e928cfa19b45854148d692247b8eb49f74c",
        750,
    ),
}

#: The four isoforms measured by dose-response, and the column holding each one's value.
ISOFORMS = ("1A2", "2C9", "2D6", "3A4")
VALUE = "CYP{isoform}_pIC50_direct_inhibition"
#: Beside each value: the 95% confidence bounds and the standard deviation.
SPREAD = ("_conf_low", "_conf_high", "_std")


def url_for(name: str) -> str:
    """Where *name* is served at the pinned commit."""
    return f"https://huggingface.co/datasets/{DATASET}/resolve/{COMMIT}/{name}"


def digest(path: Path) -> str:
    """*path*'s SHA-256, read in blocks so a large file needs no more memory than a small one."""
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            sha.update(block)
    return sha.hexdigest()


def download(name: str, into: Path) -> Path:
    """Fetch *name* at the pinned commit; the file itself is not read here."""
    into.mkdir(parents=True, exist_ok=True)
    path = into / name
    request = urllib.request.Request(url_for(name), headers={"User-Agent": "chemsuite-fetch"})
    with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as handle:
        while block := response.read(1 << 20):
            handle.write(block)
    return path


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--into", type=Path, default=DESTINATION, help="where to put the files")
    parser.add_argument(
        "--update-digests",
        action="store_true",
        help="print the digests of what was downloaded instead of checking them, to pin a commit",
    )
    arguments = parser.parse_args(argv[1:])

    print(f"{DATASET} at {COMMIT[:7]} -> {arguments.into}")
    wrong = []
    for name, (expected, rows) in FILES.items():
        path = download(name, arguments.into)
        found = digest(path)
        size = path.stat().st_size / 1e6
        if arguments.update_digests:
            print(f'  "{name}": (\n        "{found}",\n        {rows},\n    ),')
            continue
        if found != expected:
            wrong.append(f"  {name}\n    expected {expected}\n    found    {found}")
        print(f"  {name}  {size:.1f} MB  {'ok' if found == expected else 'DIGEST MISMATCH'}")

    if wrong:
        print(
            "\nThe files served at this commit are not the ones these digests were taken from.\n"
            "That is the pin doing its job: check what changed before measuring anything, then\n"
            "rerun with --update-digests to pin the new ones deliberately.\n" + "\n".join(wrong),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
