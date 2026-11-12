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

"""Rebuild chemsuite/protein/prepare/data/templates.json from the wwPDB CCD:

    python scripts/build_prepare_templates.py

The file is checked in, so users never need this. It keeps, for each residue in
`templates.BUNDLED`, every atom (hydrogens too) with its ideal coordinates, formal charge and
bonds. The CCD is in the public domain (wwPDB, CC0 1.0).
"""

from __future__ import annotations

import json
import sys
import tempfile

from chemsuite.protein.prepare import templates
from chemsuite.protein.structure.components import CcdCache


def main() -> int:
    found = {}
    with tempfile.TemporaryDirectory() as folder:
        ccd = CcdCache(folder)
        for code in templates.BUNDLED:
            block = ccd.block(code)
            template = templates.from_block(block) if block is not None else None
            if template is None:
                print(f"{code}: not found", file=sys.stderr)
                return 1
            found[code] = template.as_json()
            print(code, len(template.atoms), "atoms")
    templates.DATA.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(found, separators=(",", ":"), sort_keys=True) + "\n"
    templates.DATA.write_text(text, encoding="utf-8", newline="\n")
    print(f"wrote {len(found)} templates to {templates.DATA}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
