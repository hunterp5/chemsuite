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
"""Write ``chemsuite/macromolecule/align/data/matrices.json`` from NCBI's substitution matrices.

BLOSUM45/62/80 (Henikoff & Henikoff, Proc. Natl. Acad. Sci. USA 1992, 89, 10915) and
PAM30/70/250 (Dayhoff, Schwartz & Orcutt, Atlas of Protein Sequence and Structure 1978, 5,
suppl. 3, 345), as NCBI BLAST distributes them. The files carry no licence text; shipping them
as data, cited and pinned, was decided with the user (CLAUDE.md "Decided", 2026-11-26). Only the
numbers are kept, with each source file's SHA-256.

    python scripts/build_matrices.py path/to/folder   # holding BLOSUM45 ... PAM250

Download them from https://ftp.ncbi.nih.gov/blast/matrices/ into a folder of their own first.
Afterwards update ``SHA256`` in ``chemsuite/macromolecule/align/matrices.py``.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "chemsuite" / "macromolecule" / "align" / "data"
SOURCE = "https://ftp.ncbi.nih.gov/blast/matrices/{name}"
NAMES = ("BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250")
#: the letters every matrix is stored in (NCBI's order; B, Z, X and * as NCBI gives them)
LETTERS = "ARNDCQEGHILKMFPSTWYVBZX*"
CITATIONS = {
    "BLOSUM": "Henikoff S, Henikoff JG. Proc. Natl. Acad. Sci. USA 1992, 89, 10915-10919.",
    "PAM": "Dayhoff MO, Schwartz RM, Orcutt BC. Atlas of Protein Sequence and Structure 1978, "
    "5 (suppl. 3), 345-352.",
}


def read(path: Path) -> list[list[int]]:
    """The matrix in ``LETTERS`` order from an NCBI matrix file."""
    lines = [ln for ln in path.read_text(encoding="ascii").splitlines() if ln and ln[0] != "#"]
    header = lines[0].split()
    rows = {ln.split()[0]: [int(v) for v in ln.split()[1:]] for ln in lines[1:]}
    index = {letter: i for i, letter in enumerate(header)}
    table = [[rows[a][index[b]] for b in LETTERS] for a in LETTERS]
    assert all(table[i][j] == table[j][i] for i in range(24) for j in range(24)), path
    return table


def main(folder: str) -> None:
    matrices = {}
    for name in NAMES:
        path = Path(folder) / name
        matrices[name] = {
            "source": SOURCE.format(name=name),
            "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "citation": CITATIONS["BLOSUM" if name.startswith("BLOSUM") else "PAM"],
            "scores": read(path),
        }
    document = {
        "about": "Amino-acid substitution matrices as NCBI BLAST distributes them, in 1/2- or "
        "1/3-bit units as the source states; rows and columns in the order of 'letters'. "
        "Written by scripts/build_matrices.py.",
        "letters": LETTERS,
        "matrices": matrices,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    text = _dumps(document)
    target = OUT / "matrices.json"
    target.write_bytes(text.encode("ascii"))
    print(target, hashlib.sha256(text.encode("ascii")).hexdigest())


def _dumps(document: dict) -> str:
    """JSON with each matrix row on one line, so the file reads as a table."""
    text = json.dumps(document, indent=1)
    for entry in document["matrices"].values():
        body = ",\n".join("    " + json.dumps(row) for row in entry["scores"])
        text = text.replace(_indented(entry["scores"]), "[\n" + body + "\n   ]", 1)
    return text + "\n"


def _indented(table: list[list[int]]) -> str:
    """How ``json.dumps(document, indent=1)`` writes a matrix four levels deep."""
    inner = ["[\n" + ",\n".join("     " + json.dumps(v) for v in row) + "\n    ]" for row in table]
    return "[\n" + ",\n".join("    " + r for r in inner) + "\n   ]"


if __name__ == "__main__":
    main(*sys.argv[1:])
