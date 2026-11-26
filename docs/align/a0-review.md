# Alignment A0 review: licences, kernel speed, benchmark

Part of [ALIGN_PLAN.md](../ALIGN_PLAN.md). In progress.

## Licences (decided with the user, 2026-11-26)

Recorded in CLAUDE.md "Decided": the BLOSUM/PAM matrices ship as data with citations and
SHA-256 pins; BAliBASE, PREFAB, SABmark, HOMSTRAD and qscore are used for local measurement
only (`samples/align/`, never committed); MAFFT, FAMSA and MUSCLE 5 run only in the WSL conda
env `cs-align-bench`.

## The benchmark (in progress)

Sets: drive5's collection (`scripts/review/align_bench.py fetch`): BAliBASE 3 (386 families,
ten categories), BALIS (its 180 families with structures), SABRE (423, from SABmark). Core
columns are the reference's upper case. Families are split in halves by CRC32 of the ID;
everything is tuned on the train half and the gates are judged on the test half.

**The scorer.** Our Q and TC agree to three decimals with Edgar's `qscore` (built locally in
WSL from its source, never committed) on the five families checked.

**Comparators** (WSL env `cs-align-bench`, one thread each, run side by side): MAFFT 7.526
(L-INS-i, FFT-NS-2, `--auto`), FAMSA 2.4.1, MUSCLE 5.3.

**BAliBASE 3, train half, mean Q / TC** (n families; 2026-11-26):

| Method | RV11 (BB11) | RV12 (BB12) | all |
|---|---|---|---|
| MAFFT L-INS-i | 0.687 / 0.460 (19) | 0.943 / 0.856 (23) | 0.874 / 0.628 (170 of 193 done) |
| MAFFT `--auto` | 0.687 / 0.460 | 0.943 / 0.856 | 0.871 / 0.617 |
| MUSCLE 5 | 0.714 / 0.480 | 0.946 / 0.869 | (62 done) |
| FAMSA | 0.658 / 0.419 | 0.928 / 0.841 | 0.840 / 0.532 |
| MAFFT FFT-NS-2 | 0.580 / 0.340 | 0.914 / 0.796 | 0.830 / 0.517 |
| chemsuite, progressive + refinement (no library) | 0.623 / 0.384 | 0.926 / 0.841 | 0.831 / 0.526 |
| chemsuite, with the consistency library (A3) | 0.618 / 0.374 | 0.943 / 0.868 | 0.849 / 0.559 |

With the library, chemsuite matches L-INS-i on RV12 and passes FAMSA overall, in under a third
of L-INS-i's total time (480 s against 1,680 s for the same 193 families, ours using several
threads). **RV11 (equidistant, under 20% identity) does not move: 7 Q and 8.6 TC short of
L-INS-i**, so the proposed RV11 gate fails as things stand. The library is built from
identity-weighted pairwise alignments, which are least reliable exactly there.

What it took to get there, measured on a 58-family sample of the train half:
- The first progressive version scored Q 0.664: residues were set against gappy profile
  columns almost for free while gaps paid full price. Weighting each gap step by the occupancy
  of the column it holds fixed it (one family went from Q 0.44 to 0.96).
- BLOSUM45 with open 15 / extend 1 reached ClustalW (0.780 / 0.443 against its 0.769 / 0.450).
- Refinement on the weighted sum-of-pairs score: 0.804 / 0.491 (two passes; more gain nothing).
- The consistency library (A3), weight 4 with BLOSUM62 11 / 1: 0.832 / 0.567, against
  L-INS-i's 0.862 / 0.591 on the same sample.

## Kernel micro-benchmark (2026-11-26, Windows PC, Python 3.12, numpy 2.5.4, numba 0.68.1)

`python scripts/review/align_kernel_bench.py`: a Gotoh global alignment with traceback
pointers on random 400-residue sequences, 50 sequences (1,225 pairs), best of three.

| Kernel | One pair, 400 × 400 | All 1,225 pairs |
|---|---|---|
| numpy, row by row | 8.7 ms | 10.6 s (estimated) |
| numpy, batched over all pairs | | 7.9 s |
| numba (`njit`, plain loops) | 0.53 ms | 0.68 s (plus 1.6 s compile, once) |

The numpy and numba scores agree exactly.

**What it means for the speed gate (50 × 400 aa in 10 s or less).** In numpy, the all-pairs
alignments that the guide tree, the consistency library (A3) and the structure library need
already take 8–11 s on their own, before any profile alignment or refinement; batching over
pairs gains only 25%, because the row loop is bound by memory traffic, not Python overhead.
A consistency pass needs forward, backward and posterior passes over every pair (about three
times this), so A3 is out of reach in numpy. In numba the same all-pairs work takes 0.7 s, so
the whole aligner with consistency fits the gate with room to spare.

numba is already installed with the base dependencies (umap-learn and pacmap require it; BSD
licence), but chemsuite does not depend on it directly. Using it in the aligner would make it a
direct dependency. **Decided with the user (2026-11-26): numba is a direct dependency**
(`numba>=0.60` in `pyproject.toml`) and the DP kernels are written with `njit(nogil=True)`, so
consistency (A3) and the structure library stay in the plan.
