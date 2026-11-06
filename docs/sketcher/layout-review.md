# Sketcher P2: layout review for Clean Structure (2026-11-06)

**Question (SKETCHER_PLAN §4.6):** should Clean Structure use RDKit's own 2D layout
(`rdDepictor.Compute2DCoords`) or CoordGen (`rdCoordGen.AddCoords`, Schrödinger's open-source
layout engine shipped with RDKit)?

**Method.** `scripts/review/layout_review.py fda_approved.sdf 500` (RDKit
2026.09.1):
- the first 500 structures of the team's FDA-approved set, largest fragment only, coordinates
  removed (25.4 heavy atoms on average);
- each is laid out by both methods and by chemsuite's Clean, which runs both and keeps the
  better layout;
- every layout is scored in units of its median bond length:
  - **overlaps:** non-bonded atoms closer than 0.3;
  - **crossings:** bonds with no shared atom that cross;
  - **stretched bonds:** bonds more than 10% off the median length.

## Results

| Method | With overlaps | With crossings | With a bond off by >10% | Overlaps | Crossings | Bond spread, median / worst | ms / structure |
|---|---|---|---|---|---|---|---|
| RDKit `Compute2DCoords` (ring templates on) | 2 | 10 | 57 | 2 | 15 | 0.000 / 1.271 | 0.3 |
| CoordGen | 1 | 9 | 65 | 1 | 13 | 0.000 / 4.279 | 1.4 |
| **Best of both (chemsuite Clean)** | **1** | **6** | **48** | **1** | **8** | 0.000 / 1.271 | 2.4 |

- Neither method wins everywhere. CoordGen has fewer overlaps and crossings overall, but on
  some macrocycles it stretches a bond to four times its length. RDKit's ring templates draw
  some bridged systems better.
- Ranking each structure's two layouts (overlaps, then crossings, then stretched bonds,
  then the worst bond) and keeping the better one improves every column at ~2 ms a structure,
  which is unnoticeable for one drawing.
- The stretched bonds that remain are almost all in macrocycles and bridged rings, which no
  2D layout draws with even bonds.

**Stereo.** The script also draws each of the 500 structures, moves every atom by up to 0.3
bond lengths at random, runs Clean, and compares isomeric SMILES before and after:
- **500 of 500 kept the same molecule and stereo; none was refused.**

**Decision.** Clean Structure (`chemsuite/sketch/layout.py`) works as follows:
1. Lay the drawing out both ways and rank the two layouts as above.
2. Turn and move (never mirror) each cleaned fragment onto where it was drawn, by a
   least-squares rigid fit (Kabsch).
3. Keep the drawn wedges if the molecule reads the same.
4. Otherwise let RDKit place wedges (`WedgeMolBonds`) for the stereo read before the clean.
5. Try the next layout if neither keeps the molecule. If none does, leave the drawing
   unchanged and say why.

The tests are `tests/sketch/test_layout.py`: drug cases with stereocentres, E/Z bonds, a
bridged bicycle and a salt; a drawing with no sensible coordinates; place and turn kept; no
mirroring; re-wedging.
