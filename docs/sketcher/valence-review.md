# Sketcher P0: valence and implicit hydrogens (2026-11-05)

**Question.** Which rule should the sketcher use to work out implicit hydrogens and flag bad
valences?

**Candidates:**
1. Ketcher's model: the BIOVIA CTfile valence model, as implemented in `Atom.calcValence`
   (`ketcher-core/src/domain/entities/atom.ts` at `2b031cc2`) in its default mode.
2. RDKit's model, the one chemsuite already uses for the table and search.

**Reproduce:** `python scripts/review/valence_review.py`. The script contains a line-by-line
Python port of Ketcher's rules, for review only.

## Method

1. **Reference species.** 61 real species whose hydrogen count on one atom is not in doubt.
   Each is held the way a sketch holds it: heavy atoms, charges and bonds, with no hydrogen
   information.
   - Both models must supply the hydrogens.
   - RDKit is asked the way chemsuite will ask it: the sketch is written as a V2000 molfile
     without hydrogen counts and read back with sanitisation. That step normalises hypervalent
     drawings such as nitro and perchlorate.
   - The set covers neutral organics, onium ions, anions, radicals and carbenes, boron and
     aluminium ate complexes, hypervalent P, S and halogens, common counter-ions, and
     organometallics drawn with covalent bonds.
2. **Exhaustive sweep.** 35 elements × charges −2 to +2 × 0–2 radical electrons × bond-order
   sums 0–6, which is 3,675 cases, to map where the models differ.

## Results

**Reference species: RDKit 59/61, Ketcher 57/61.**

| Species | Drawn as | Expected H | Ketcher | RDKit |
|---|---|---|---|---|
| hexafluorosilicate | `F[Si-2](F)(F)(F)(F)F` | 0 | error ✗ | 0 |
| hexafluorophosphate | `F[P-](F)(F)(F)(F)F` | 0 | error ✗ | 0 |
| hexafluoroarsenate | `F[As-](F)(F)(F)(F)F` | 0 | error ✗ | 0 |
| hexafluoroantimonate | `F[Sb-](F)(F)(F)(F)F` | 0 | error ✗ | 0 |
| chlorine trifluoride | `FCl(F)F` | 0 | 0 | error ✗ |
| iodine heptafluoride | `FI(F)(F)(F)(F)(F)F` | 0 | 0 | error ✗ |

Every other species agreed with the expected count in both models.

**Sweep:** 1,752 of 3,675 cases differ.
- Nearly all are metals, multiply charged atoms or radicals on heavy main-group elements.
  Ketcher returns "bond-order sum, no hydrogens" for most metals; RDKit has element-specific
  valence lists.
- For C, O, S, B and Se at charges −1 to +1 without radicals, the models agree in every case.
- The organic cases that differ are hypervalent N, P and halogens. For example, Ketcher gives
  a neutral nitrogen with four bonds one hydrogen (valence 5), where RDKit reports an error.

## Decision

**The sketcher uses RDKit's valence model**, run through the same molfile and sanitisation path
the rest of chemsuite uses.

1. **It gets more real species right.** Ketcher's failures include PF₆⁻, the counter-ion in
   HATU and HBTU, two coupling reagents medicinal chemists draw daily. RDKit's two failures
   (ClF₃, IF₇) are interhalogens a medicinal chemist will rarely draw.
2. **It is consistent.** A drawing is valid in the sketcher exactly when the table, search and
   descriptors can use it. With Ketcher's model, some drawings would look fine in the sketcher
   and then fail in the table.

**How it shows while drawing.** A bad valence is shown, not refused: an atom that RDKit rejects
is marked, and the reason appears on hover. Ketcher does the same with its `badConn` flag.

**Exporting a molfile** writes hydrogen counts explicitly for any atom where the two models
differ, so other programs (ChemDraw, Marvin) cannot read it differently.

**Known gaps** (recorded in tests as expected failures):
- ClF₃ and IF₇ are rejected by RDKit.
- Ketcher's "pentavalent nitro" drawing is accepted, and normalised to `[N+](=O)[O-]` on
  sanitisation. That is correct chemistry, but the user sees the drawing change. The sketcher
  says so.

**What we keep from Ketcher here:** its interaction, marking a bad atom rather than refusing it.
We do not use its valence table.
