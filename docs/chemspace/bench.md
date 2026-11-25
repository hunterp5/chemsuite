# TMAP: above 49k, MHFP6, and beside tmap-silicon (ROADMAP §7.0d)

Measured 2026-11-19 on the Windows PC (Windows 10, Python 3.12, RDKit 2026.09.2, 8 workers).
Reproduce with `scripts/review/chemspace_bench.py` and `scripts/review/tmap_compare.py` (their
docstrings give the steps).

**The sets.** Up to 49,189 compounds: BindingDB's distinct compounds (the local sample holds no
more, its full 502 MB SDF included). Above that: the first *n* distinct compounds of ChEMBL 37's
`chembl_37_chemreps.txt.gz` (CC BY-SA 3.0, local test data under `samples/chembl/`, never
committed).

## TMAP above 49k

k = 10, every structure valid. The TMAP column excludes the fingerprints (ECFP4 takes seconds;
MHFP6 is timed in its own column); peak memory is the whole process's, the fingerprints
included (100,000 measured before the change below, 49,189 and 250,000 after).

| compounds | fingerprint | fingerprints | neighbours | layout | TMAP | peak memory |
|---|---|---|---|---|---|---|
| 49,189 (BindingDB) | ECFP4 | | 6 s | 2.5 s | 9 s | 382 MB |
| 49,189 (BindingDB) | MHFP6 | 41 s | 13 s | 3 s | 16 s | |
| 100,000 | ECFP4 | | 14 s | 6 s | 20 s | 709 MB |
| 100,000 | MHFP6 | 61 s | 26 s | 6 s | 32 s | 819 MB |
| 250,000 | ECFP4 | | 36 s | 15 s | 52 s | 1,175 MB |
| 250,000 | MHFP6 | 152 s | 63–66 s | 15 s | 82 s | 1,516 MB |

These times are from before the layout change of 2026-11-24, which takes the layout column
to about twice as long (27 s at 250,000 on ECFP4; see "The layout, 2026-11-24" below).

The forest's candidate table (rows × 16 trees × 20 places) was int64; it is int32 now, which
took the 250k peaks down from 1,496 and 1,836 MB with the same trees (same weights) and times.
There is no hard limit: 250,000 compounds take under a minute on ECFP4 with about 1.2 GB, and
about four minutes on MHFP6, most of it encoding.

## MHFP6: how many values

MHFP6 is RDKit's encoder (`rdMHFPFingerprint`, radius 3), its MinHash values put straight into
the forest. Recall is against every pair's agreement for 300 sampled rows of BindingDB's 49,189.

| values | encoding | neighbours | recall of the ten nearest | standard error at 0.37 |
|---|---|---|---|---|
| 128 | 28 s | 9.7 s | 0.931 | 0.043 |
| 256 | 41 s | 18 s | 0.925 | 0.030 |
| 512 | 39 s | 31 s | 0.921 | 0.021 |

Recall against each count's own estimate barely moves; what more values buy is a less noisy
estimate of the true overlap. 256 is kept: half the forest time of 512, and noise a third below
128's. Encoding is MHFP6's cost: about 6 ms a molecule in RDKit's encoder (it writes each
atom's substructures as SMILES), against well under a millisecond for ECFP4, so MHFP6 TMAPs of
large sets take minutes, not seconds. The user guide says so.

## Beside tmap-silicon

`tmap-silicon` 1.0.19 (the maintained build of Probst and Reymond's `tmap`: LSH Forest, its MST
and OGDF's FM³ layout) in a virtual environment of its own, given the same 256 MHFP6 values of
BindingDB's 49,189. It counts each compound among its own nearest, so it was asked for 11 to
give 10 others (asked for 10 it had 9, which first made it look worse than it is).

| | chemsuite | tmap-silicon, 16 trees | tmap-silicon, 64 trees |
|---|---|---|---|
| recall of the ten nearest | 0.925 | 0.618 | 0.721 |
| tree weight (Σ 1 − agreement; lower: closer neighbours joined) | 21,568 | 21,867 | 21,752 |
| separate trees | 1 | 1 | 1 |
| neighbours | 12.6 s | 26.1 s | 18.6 s |
| whole TMAP | 15.5 s (18.6 s after the 2026-11-24 layout change) | 53.5 s | 36.4 s |

chemsuite finds more of the true neighbours (it looks ten places either side in every band,
then among neighbours' neighbours), builds a tree that joins closer compounds, and is about
2–3× faster.

**The layout was where tmap-silicon was better** (fixed 2026-11-24, below). Its pictures
(`samples/chemspace/compare_tmap_silicon_l64_49189.png` against
`compare_chemsuite_49189_before.png`) spread the tree over the page with every branch apart.
chemsuite's multilevel force layout packed branches into dense clumps joined by long links,
crossed branches over each other, and left part of the page empty. The trees are nearly the
same; the drawing was not.

## The layout, 2026-11-24

**Cause.** Repulsion reached only 3 edge lengths, at every level. Coarsening was healthy (each
level halves the tree, down to 32 nodes on the 49k tree), so the coarse layouts had few nodes
but could not push branches apart beyond that reach.

**Fix** (`chemsuite/analysis/graphs/layout.py`):
- repulsion reaches 16 edge lengths on the coarse levels, where the branches are arranged and
  the nodes are few, and 8 on the finest;
- each finer level starts every node on its group's place moved halfway toward the groups its
  edges lead to, instead of on the group's place with noise.

**Measured** on the same trees, old and new back to back with nothing else running. Crossings
are counted on 8,000 sampled edges (each crossing pair once). Coverage is the share of a
100 × 100 grid over the picture's square that holds a compound.

| tree | layout, old → new | crossings per edge, old → new | coverage, old → new |
|---|---|---|---|
| 49,189 MHFP6 | 3.1 → 6.5 s | 0.059 → 0.026 | 0.14 → 0.31 |
| 49,189 ECFP4 | 3.0 → 5.6 s | 0.060 → 0.030 | 0.11 → 0.23 |
| 249,997 ECFP4 | 15.9 → 26.9 s | 0.049 → 0.024 | 0.13 → 0.27 |

tmap-silicon (64 trees) on the 49,189 MHFP6 values: 0.020 crossings per edge, coverage 0.31.
Peak memory is unchanged (1,177 MB at 250,000 on ECFP4, 1,517 MB on MHFP6).

Chemical Space Networks share the layout. At 0.65 (5,000 compounds, largest part 84) the
pictures are unchanged. At 0.45 (largest part 1,230) the clusters branch off the spine instead
of overlapping it; the layout takes 2.0 s instead of 1.2 s.

**Tried and dropped:**
- unlimited repulsion (the far field on an FFT grid, the near field exact): everything was
  driven to a hollow ring with stretched edges, crossings rose to 0.088 per edge, and it was
  seven times slower;
- the long reach on the finest level as well (16): no better than 8 there, about 30 % slower;
- 32 on the coarse levels: coverage fell and edges stretched;
- more steps per level: within the noise between seeds, for up to twice the time.

Pictures: `compare_chemsuite_49189.png` (now) beside `compare_tmap_silicon_l64_49189.png`;
`tmap_ecfp4_249997.png` and `tmap_mhfp6_249997.png` with their `_before` copies.

**Still open:** at 250,000 the 50 or so small separate trees are packed in a thin strip beside
the giant one (`pack` lines them up after it), so they are hard to see.
