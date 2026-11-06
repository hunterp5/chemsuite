# Sketcher P0: licence review of Ketcher material (2026-11-05)

**Reference commit:** github.com/epam/ketcher `2b031cc214364fb0be7fe42c902ea772e981fc34`
(master, 2026-11-05). Every rule and data file the sketcher takes from Ketcher is cited at this
commit.

## What we take, and on what terms

| Material | Path | Licence evidence | Verdict |
|---|---|---|---|
| Chemistry and editing rules, re-implemented in Python (no code copied) | `packages/ketcher-core/src/**`, `packages/ketcher-react/src/script/**` | Each file carries EPAM's Apache-2.0 header ("Copyright 2021 EPAM Systems, Licensed under the Apache License, Version 2.0"). The repository's LICENSE is Apache-2.0. | **OK.** Cite the source file beside each rule and keep Ketcher's NOTICE in `chemsuite/sketch/data/`. |
| CTfile default valence model, **ported line by line** (P1) | `ketcher-core/src/domain/entities/atom.ts` (`calcValence`), `domain/constants/elements.ts` | Same Apache-2.0 header | **OK.** `chemsuite/sketch/ctfile.py` keeps EPAM's copyright and licence notice, states that it was translated and cut down (Apache-2.0 §4(b)), and points to `chemsuite/sketch/data/LICENSE-Ketcher.txt` and `NOTICE-Ketcher.txt`. |
| Functional groups, 62 records | `packages/ketcher-react/src/templates/fg.sdf` | In the Apache-2.0 repository with no separate licence. File history starts 2022-12-22, when the files moved into `ketcher-react` (#1983); last change 2023-07-06. | Licence OK, but **not shipped**: P2 found 9 unusable records and a wrong Mes (`abbreviation-review.md`). chemsuite uses its own verified table. |
| Salts and solvents, 135 records | `…/templates/salts-and-solvents.sdf` | Same; last change 2025-03-31 | Same |
| Structure templates, 276 records | `…/templates/library.sdf` | Same; last change 2026-10-14 | Same |

Ketcher's NOTICE lists other licences (MIT, BSD, LGPL-2.1, MPL and so on) for its *dependencies*.
We take no dependency code, so they do not apply.

Apache-2.0 is compatible with chemsuite's GPL-3.0 when chemsuite is distributed as a whole.
Apache-2.0 section 4 requires us to:
- keep the licence and NOTICE with the data;
- state that files were changed if we change them.

## Data quality: Ketcher's data must be checked, not trusted

The MolManager review (`molmanager-review.md` §2) found that Ketcher's **"Mes"** record is not
mesityl. It has 8 carbons and attaches through a CH₂, which makes it 3-methylbenzyl.
Mesityl is 2,4,6-trimethylphenyl, C₉H₁₁.

So, before any record ships (phase P2):
- each abbreviation is checked against an independent definition: its formula, and its
  attachment atom given by the name;
- corrections are recorded in a changes file beside the data, as Apache-2.0 §4(b) requires.

**Format note:** Ketcher's SDF files put a blank line after each `$$$$`. That is not standard
SDF, and a strict reader takes the blank line as the record name. Our loader must allow for it.
