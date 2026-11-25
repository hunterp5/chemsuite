# L0 review: the loop benchmark

2026-11-25. Plan: [LOOPS_PLAN.md](../LOOPS_PLAN.md), milestone L0. Script:
`scripts/review/loop_bench.py`; data under `samples/loops/` (not committed). The full tables are
`samples/loops/l0-test.md` (`loop_bench.py score --test`).

L1 (`prepare/loops.py`) was built alongside the benchmark rather than after it, so this review
scores it too. Its numbers also stand in for L3's table, which the plan defines as this same
benchmark on the test half. L3 still needs its run by hand in the app.

## Set

- 560 cuts: 80 at each length (2, 3, 4, 6, 8, 10 and 12 residues), taken from high-resolution
  (≤ 1.8 Å), nonredundant chains (one per 30% sequence cluster), away from chain ends.
- Train/test halves by CRC32 of the ID. **All numbers below are from the test half, 285 cuts.**
  The Ramachandran grids were built from the train half only.
- Methods:
  - **caps:** Prepare as it is today, ACE/NME on both ends; the floor;
  - **PDBFixer:** `addMissingResidues`, run in the OpenMM env, benchmark only;
  - **ours:** `loops.py`, then Prepare;
  - each of the three again after our Minimize (`-min`), with the built atoms and the cut's
    residues free.

## Results (test half)

Median backbone RMSD (N, CA, C, O) to the deposit, in Å, with no superposition: the anchors are
fixed, so none is needed. The ≤ 2 Å column is the share of loops within 2 Å.

| length | n | ours | ours-min | PDBFixer | PDBFixer-min | ours-min ≤ 2 Å | PDBFixer-min ≤ 2 Å |
|---|---|---|---|---|---|---|---|
| 2 | 42 | 0.29 | 0.22 | 2.03 | 1.65 | 100% | 58% |
| 3 | 40 | 0.41 | 0.25 | 2.44 | 2.41 | 100% | 32% |
| 4 | 34 | 0.56 | 0.45 | 3.47 | 3.80 | 91% | 15% |
| 6 | 43 | 1.67 | 1.57 | 4.05 | 4.08 | 62% | 2% |
| 8 | 37 | 2.83 | 2.41 | 4.62 | 4.80 | 38% | 3% |
| 10 | 46 | 5.14 | 5.24 | 5.64 | 5.69 | 15% | 2% |
| 12 | 43 | 5.23 | 5.34 | 6.77 | 6.74 | 2% | 0% |
| all | 285 | 1.52 | 1.38 | 3.99 | 4.05 | 57% | 16% |

| measure (all lengths) | ours | ours-min | PDBFixer | PDBFixer-min |
|---|---|---|---|---|
| built | 285/285 | 283/285 | 272/285 | 269/285 |
| any severe clash (< 2.2 Å) | 13% | 0.4% (1 loop) | 95% | 3% |
| peptide-bond error > 0.1 Å | 18% | 0% | 32% | 2% |
| Ramachandran outliers (per residue) | 0.0% | 1.1% | 35.8% | 25.0% |
| median seconds | 5.3 | 5.4 | 3.1 | 8.8 |

By secondary structure, ours-min does best on strands (median 0.34 Å, all within 2 Å), then
helices (median 0.23 Å, 74% within 2 Å) and coil (0.52 Å, 69%). It is weakest where a cut spans several
classes ("mixed", 158 cuts: 2.66 Å, 41%), which is where most of the long cuts fall.

## What it shows

- **Ours beats PDBFixer at every length**, before and after Minimize. The gap is widest where it
  matters most, at 2–6 residues: 0.2–1.6 Å against 1.7–4.1 Å.
- **PDBFixer's loops are unusable without Minimize** (95% clash) and still poor after it (a
  quarter of its loop residues are Ramachandran outliers). It is not a target worth matching.
- **Every loop closes.** Ours built all 285. Our two `-min` failures are Minimize refusing the
  structure itself, which caps hit in the same way (see Failures).
- **Minimize finishes the job.** It removes the junction strain (peptide errors 18% → 0%) and
  all but one clash (13% → 1 loop).
- **Long loops close cleanly but are guesses.** At 10–12 residues the loops are clash-free and
  Ramachandran-clean, but the median RMSD is about 5 Å and a few are 8–11 Å from the deposit
  (11 of the 89 loops of length 10–12). They fill the hole and fix the electrostatics. They do
  not say where the loop really is.
- **Drift is by design.** Before Minimize, the only backbone atom outside the loop that ours
  moves is anchor *a*'s O (median 0.22 Å, max 2.17 Å). The plan frees ψ of *a* and rebuilds
  that O in the new peptide plane; every other atom stays exactly where it was. After Minimize,
  the rest of the chain moves up to 2 Å, as it does with caps alone (1.86 Å): that comes from
  Minimize, not from the loops.

## Gates

The plan's L3 gates, applied to ours-min at each length:

| length | all closed | no severe clash | Rama outliers < 5% | median RMSD ≤ PDBFixer-min's | passes |
|---|---|---|---|---|---|
| 2 | yes | yes | 0.0% | 0.22 ≤ 1.65 | yes |
| 3 | yes | yes | 0.0% | 0.25 ≤ 2.41 | yes |
| 4 | yes | yes | 0.5% | 0.45 ≤ 3.80 | yes |
| 6 | yes | yes | 0.6% | 1.57 ≤ 4.08 | yes |
| 8 | yes | yes | 0.8% | 2.41 ≤ 4.80 | yes |
| 10 | yes | **no** (1 of 46: `1CS6_A_166_10`, 3 severe pairs) | 2.4% | 5.24 ≤ 5.69 | **no** |
| 12 | yes | yes | 1.2% | 5.34 ≤ 6.74 | yes |

"All closed" counts the loop builder's own output. The 4PUQ and 1I1J cuts that Minimize refused
fail the same way with caps.

**Recommended `loop_max` default for L2: 8**, the longest length up to which every length
passes. Lengths 10 and 12 should stay available in the dialog (the plan's spin box goes to 12),
with the RMSD caveat above in its tooltip.

What L1's unit tests should hold to, set from these numbers so that a regression shows:
- closure under 0.1 Å on the three anchor atoms;
- no severe clash after side chains are placed, for the fixture's loop;
- no Ramachandran outliers in the built residues;
- RMSD under 2 Å for a cut of length ≤ 4 on the fixture.

## Failures (7 of 855 minimizations)

None come from the loop builder.

- **4PUQ_A_40_12 (caps, PDBFixer, ours):** the entry has DNA, and Minimize cannot simulate
  nucleic acids yet. Expected.
- **1I1J_A_10_6 (caps, ours):** "B/MET 3: its atoms do not match the NMET template (15 atoms,
  the template 19)". A different chain's N-terminal Met, not the loop. This looks like a
  Prepare/Minimize bug with a partial N-terminal residue, worth its own look.
- **8U79_A_129_3, 1XG4_A_114_10 (PDBFixer only):** OpenMM returned NaN coordinates, from
  PDBFixer's clashing loops.

## Open

- **Lone residues (m = 1).** The benchmark has no length-1 cuts, so whether `_fill` gives way to
  `loops.py` for single gaps is still undecided. A length-1 set costs one `cut` run and one
  `loops` run.
- **The 1I1J Met template mismatch** above.
- **L2** (Prepare, Inspect and the dialog), and **L3**'s by-hand run of Prepare + Minimize on an
  8-residue gap in the app.
