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
"""Bundle the AutoDock Vina binaries chemsuite ships (``chemsuite/protein/dock/bin``).

Vina is Apache-2.0, which a GPL-3-or-later project may ship as long as its licence and notice
travel with it, so ``LICENSE-AutoDock-Vina.txt`` and ``NOTICE`` sit beside the binaries.

The binaries are checked in; this script only puts them there, so bundling is reproducible and
the version is easy to bump. Every file is verified against the SHA-256s pinned in
``bin/checksums.json``::

    python scripts/fetch_vina.py              # download what is missing and verify it
    python scripts/fetch_vina.py --check      # verify what is there, download nothing
    python scripts/fetch_vina.py --refresh    # pin the hashes of what the release serves now
    python scripts/fetch_vina.py --version 1.2.8 --refresh   # bump the version

``--refresh`` trusts the download, so read the release page before using it: it is how a new
version (or a first bundling) is pinned, not a way past a failed check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "chemsuite" / "protein" / "dock" / "bin"
PINS = BIN / "checksums.json"
RELEASE = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v{version}/{asset}"
LICENCE = "https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/v{version}/LICENSE"
DEFAULT_VERSION = "1.2.7"
#: where each platform's binary goes, and the release asset it comes from
PLATFORMS = {
    "win_amd64/vina.exe": "vina_{version}_win.exe",
    "macos_arm64/vina": "vina_{version}_mac_aarch64",
    "linux_x86_64/vina": "vina_{version}_linux_x86_64",
}
NOTICE = """\
AutoDock Vina
Copyright 2021 Diogo Santos-Martins, Jerome Eberhardt, Stefano Forli and others
Copyright 2006-2010 Oleg Trott, The Scripps Research Institute

Licensed under the Apache License, Version 2.0 (see LICENSE-AutoDock-Vina.txt).

The binaries in this folder are the unmodified releases of AutoDock Vina {version}, from
https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v{version}. chemsuite itself is
GPL-3.0-or-later and does not link against Vina: it runs these programs and reads their output.

If you use docking results from chemsuite, cite:

  J. Eberhardt, D. Santos-Martins, A. F. Tillack, S. Forli. AutoDock Vina 1.2.0: New Docking
  Methods, Expanded Force Field, and Python Bindings. J. Chem. Inf. Model. 2021, 61, 3891-3898.

  O. Trott, A. J. Olson. AutoDock Vina: improving the speed and accuracy of docking with a new
  scoring function, efficient optimization and multithreading. J. Comput. Chem. 2010, 31, 455-461.
"""


def pins() -> dict:
    if not PINS.exists():
        return {"version": DEFAULT_VERSION, "files": {}}
    return json.loads(PINS.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download(url: str, into: Path) -> None:
    """Fetch *url* to *into* with curl (its TLS store is the system's, which Python's is not
    on Windows)."""
    into.parent.mkdir(parents=True, exist_ok=True)
    part = into.with_suffix(into.suffix + ".part")
    command = ["curl", "-sSfL", "--retry", "3", "--max-time", "300", "-o", str(part), url]
    done = subprocess.run(command, capture_output=True, text=True)
    if done.returncode != 0:
        part.unlink(missing_ok=True)
        raise SystemExit(f"could not download {url}: {done.stderr.strip() or done.returncode}")
    part.replace(into)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="")
    parser.add_argument("--check", action="store_true", help="verify only, download nothing")
    parser.add_argument("--refresh", action="store_true", help="pin what the release serves now")
    args = parser.parse_args(argv[1:])
    pinned = pins()
    version = args.version or pinned["version"]
    changed = version != pinned["version"]
    if changed and not args.refresh:
        raise SystemExit(f"{version} is not the pinned version: use --refresh to bump it")
    files: dict[str, str] = {} if changed else dict(pinned["files"])
    problems = []
    for name, asset in PLATFORMS.items():
        path = BIN / name
        want = files.get(name, "")
        if not path.exists():
            if args.check:
                problems.append(f"{name}: missing")
                continue
            download(RELEASE.format(version=version, asset=asset.format(version=version)), path)
        found = digest(path)
        if args.refresh or not want:
            files[name] = found
            print(f"pinned {name} {found[:16]}… ({path.stat().st_size:,} bytes)")
        elif found != want:
            problems.append(f"{name}: sha256 {found} is not the pinned {want}")
        else:
            print(f"ok {name} ({path.stat().st_size:,} bytes)")
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    if args.check:
        return 0
    licence = BIN / "LICENSE-AutoDock-Vina.txt"
    if not licence.exists() or args.refresh:
        download(LICENCE.format(version=version), licence)
    (BIN / "NOTICE").write_text(NOTICE.format(version=version), encoding="utf-8")
    PINS.write_text(json.dumps({"version": version, "files": files}, indent=2) + "\n",
                    encoding="utf-8")  # fmt: skip
    print(f"AutoDock Vina {version} bundled in {BIN.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
