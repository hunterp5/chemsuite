# Reactions R1: round-trip review on USPTO-50k (2026-11-06)

**Question (REACTIONS_PLAN R1):** does a reaction keep everything when it goes from reaction
SMILES into the sketcher, then out to an RXN file, and back again? That means:
- the same components on each side, in order;
- the same atom maps;
- the same stereo.

**Method.** `scripts/review/reaction_review.py data_processed.csv`, with RDKit 2026.09.1.
- **Data:** USPTO-50k, from github.com/connorcoley/retrosim, MIT. These are the 50,016
  atom-mapped patent reactions of Schneider et al. (J. Chem. Inf. Model. 2016, 56, 2336),
  drawn from Lowe's CC0 extraction.
- **Each `rxn_smiles` goes through four steps:**
  1. `chemsuite.sketch.reaction.from_text` draws it;
  2. the drawing is written back as reaction SMILES;
  3. the drawing is written as an RXN file;
  4. the RXN file is drawn and written again.
- **Pass:** both results equal RDKit's own reaction SMILES of the input.

## Results

| Result | Reactions |
|---|---|
| Identical reaction SMILES, both ways | **49,986** (99.94%) |
| Same molecules by InChIKey; RDKit writes the aromatic form differently | 27 |
| Different | 3 |

- **The 27.** All are large fused heteroaromatics, such as thiazino-fused and
  pyrrolo-fused pyridones.
  - The drawing holds Kekulé bonds. Read back, RDKit perceives a different but equivalent
    aromatic form.
  - The molecules are identical by InChIKey.
- **The 3.** These are explicit hydrogen counts on phosphorus that no valence rule gives
  (`[PH4]`, `[PH2]=O`).
  - The sketcher draws hydrogens from its valence model (docs/sketcher/valence-review.md),
    so it cannot hold them.
  - Two of the 3 are the same reaction recorded twice.

## Bugs the review found, now fixed

1. **RXN files lost aromatic NH hydrogens.**
   - RDKit wrote aromatic bond type 4, which drops the H of an [nH] in fused systems (an
     azaindole), and the file could not be read back.
   - RXN files are now written with Kekulé bonds, as molfiles are (`reaction.to_rxn_block`).
2. **E/Z was invented from the layout.** This also affected plain molecules pasted as
   SMILES: `CC=CC` came back as `C/C=C/C`.
   - A double bond that SMILES leaves unspecified is now drawn "either" (CTfile stereo 3),
     as RDKit's molfile writer marks it (`convert._undefined_double_bonds`).
   - Atom maps count when deciding whether a bond can be E/Z.
   - Cumulated bonds (azides, allenes) never can.
3. **Reaction templates lost E/Z and stereocentres.** RDKit reads reaction templates without
   stereo perception, so each template is now read again as a molecule
   (`reaction._template_sketch`).
4. **A lone ion such as `[H+]` got NaN coordinates.** CoordGen, which the app prefers
   globally, places a single atom at NaN. RDKit's own layout is used then.
5. **Bridged stereocentres were lost when laying out text.** For a quinuclidinium,
   RDKit's wedges on its preferred layout read back as contradictory. When the drawing does
   not read back with the molecule's stereo, `convert.from_mol` now tries the other layout.

Tests: `tests/sketch/test_reaction.py` and `tests/sketch/test_layout.py`
(`test_text_without_a_drawing_keeps_bridged_stereo`).
