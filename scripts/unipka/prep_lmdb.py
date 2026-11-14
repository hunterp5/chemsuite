#!/usr/bin/env python3
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

"""Uni-pKa's preprocess_pka.py with a time limit per row (docs/UNIPKA_PLAN.md, U2).

    prep_lmdb.py <rows.tsv> <out.lmdb> [--workers 6] [--timeout 180]

Runs in WSL, in the training env (~/bench/unipka/.venv, with pebble). Each row's conformers
come from upstream's own ``smi2coords``, so the recipe is unchanged; a row that takes longer than
*timeout* seconds (a large macrocycle such as rifampicin retries each failed embed 5,000 times)
is dropped and listed in ``<out>.dropped.tsv``. Rows keep their input order, which the
cross-validation split depends on, and the lmdb is committed every 500 rows.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from concurrent.futures import TimeoutError as FutureTimeout
from pathlib import Path

import lmdb
from pebble import ProcessPool

sys.path.insert(0, str(Path.home() / "bench/unipka/src/scripts"))
from preprocess_pka import smi2coords  # upstream's, Apache-2.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("rows", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=180.0)
    args = parser.parse_args()

    with args.rows.open(encoding="utf-8") as f:
        rows = [(r["SMILES"], float(r["TARGET"])) for r in csv.DictReader(f, delimiter="\t")]
    args.out.unlink(missing_ok=True)
    env = lmdb.open(str(args.out), subdir=False, lock=False, readahead=False, meminit=False,
                    max_readers=1, map_size=int(100e9))  # fmt: skip
    dropped: list[tuple[int, str, str]] = []
    written = 0
    txn = env.begin(write=True)
    with ProcessPool(max_workers=args.workers) as pool:
        futures = [pool.schedule(smi2coords, (row,), timeout=args.timeout) for row in rows]
        for i, (future, row) in enumerate(zip(futures, rows, strict=True)):
            try:
                blob = future.result()
            except FutureTimeout:
                blob, why = None, "timeout"
            except Exception as error:  # noqa: BLE001  (a worker that died: drop the row)
                blob, why = None, type(error).__name__
            else:
                why = "failed"
            if blob is None:
                dropped.append((i, why, row[0]))
            else:
                txn.put(f"{written}".encode("ascii"), blob)
                written += 1
            if (i + 1) % 500 == 0 or i + 1 == len(rows):
                txn.commit()
                txn = env.begin(write=True)
                print(f"{i + 1}/{len(rows)} rows, {written} written, {len(dropped)} dropped",
                      flush=True)  # fmt: skip
    txn.commit()
    env.close()
    with Path(f"{args.out}.dropped.tsv").open("w", encoding="utf-8", newline="") as f:
        csv.writer(f, delimiter="\t").writerows([("row", "why", "SMILES"), *dropped])
    print(f"done: {written} written, {len(dropped)} dropped", flush=True)
    os._exit(0)  # skip joining workers stuck in a C embed


if __name__ == "__main__":
    main()
