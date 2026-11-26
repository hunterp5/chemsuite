# Alignment A0 review: licences, kernel speed, benchmark

Part of [ALIGN_PLAN.md](../ALIGN_PLAN.md). In progress.

## Licences (decided with the user, 2026-11-26)

Recorded in CLAUDE.md "Decided": the BLOSUM/PAM matrices ship as data with citations and
SHA-256 pins; BAliBASE, PREFAB, SABmark, HOMSTRAD and qscore are used for local measurement
only (`samples/align/`, never committed); MAFFT, FAMSA and MUSCLE 5 run only in the WSL conda
env `cs-align-bench`.

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
