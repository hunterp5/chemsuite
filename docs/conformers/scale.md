# Hundreds of conformers per molecule: what they cost (ROADMAP §7.0f)

Measured 2026-11-17 on the Windows PC (Windows 10, Python 3.12, RDKit 2026.09.2, 8+ cores),
before and after the fixes below. Reproduce with:

    python scripts/review/conformer_scale.py data      # once: the three documents below
    python scripts/review/conformer_scale.py run       # the Qt-free steps
    pytest -m bench -s tests/app/table/test_bench_conformers.py   # the app steps

**The set.** 1,000 distinct drug-like molecules from BindingDB (20–45 heavy atoms, 3–10
rotatable bonds, neutral; seed 2026), with CONFORGE conformers asked for at 10, 100 and 300
per molecule (no energy window). CONFORGE's default pruning (0.5 Å) gave 10, 79 and 190 on
average (at most 10, 100 and 300). The 3D column holds 7, 55 and 134 MB of blobs; the
documents are 4, 23 and 53 MB on disk.

## Against the targets

"1k × 100" and "1k × 300" are the 100 and 300 sets. RSS is the process's memory over the
blobs themselves, with every molecule of the 3D column loaded once.

| what | target | 1k × 100 before → after | 1k × 300 before → after |
|---|---|---|---|
| first screen of the 3D column (cold) | < 300 ms | 128 → 147 ms ✓ | 139 → 135 ms ✓ |
| a page down, pictures not cached | < 150 / 200 ms | 100 → 91 ms ✓ | 108 → 100 ms ✓ |
| sort by another column and repaint | as with 10 conformers | 340 → 128 ms (10: 123) ✓ | 317 → 133 ms ✓ |
| memory over the blobs | < 300 / 400 MB | — → 200 MB ✓ | 434 → 238 MB ✓ |
| open / save the document | < 5 s / < 10 s | 0.24 / 1.1 s ✓ | 0.70 / 2.5 s ✓ |
| longest UI stall while saving | < 100 ms | 49 → 27 ms ✓ | 69 → 74 ms ✓ |
| Browser: step a conformer | < 50 ms | 2.9 → 1.6 ms ✓ | 5.8 → 2.3 ms ✓ |
| Browser: move to the next row | < 150 / 200 ms | 3.7 → 1.9 ms ✓ | 4.2 → 2.7 ms ✓ |
| Browser: Overlay (UI thread) | < 100 ms | 220 → 26 ms ✓ | 559 → 51 ms ✓ |
| apply a new 1k-row 3D column / undo / redo | < 2 s / < 300 ms | 37 / 13 / 14 ms ✓ | 36 / 11 / 12 ms ✓ |
| substructure search on the 3D column | ≤ 2× the 2D column | 152 vs 12 ms ✗ | 416 vs 14 ms ✗ |
| export every conformer (SDF, a job) | none set | 48 → 15 s | 114 → 37 s |

Generating them is not the problem: CONFORGE made the 1k × 300 set in 153 s on 8 workers
(ETKDG would take about 40 s a molecule at 300).

## What changed

- **The molecule cache is limited by size** (`LoadedCache` in `core/table/columns.py`): at
  most 1,024 molecules and 64 MB of blobs. It was 1,024 molecules whatever their size; a
  190-conformer molecule takes about 3× its blob in memory, so a full cache of the 1k × 300
  set held about 430 MB.
- **Overlays work out the bond graph once per molecule** (`molview.geometry.topology`): the
  bond list, the drawn orders and the plane of each double bond were worked out again for every
  conformer drawn, and every bond was a Python loop. Now single bonds are built at once and
  double and triple bonds together; the drawing is identical (checked part for part against
  the previous code, 240 cases).
- **Exporting every conformer** writes a V2000 molfile once and puts each conformer's
  coordinates in it: `MolToMolBlock` copied the whole ensemble for each of the 190,000
  records. Molecules with wedge flags (which RDKit derives from each conformer's geometry) are
  written from a one-conformer copy. The file is byte for byte what the previous code wrote.
  The export now checks for cancel and reports progress every 16 rows instead of every 2,048.

## Still above target: 2D work on a 3D column

Anything that needs only the 2D structure (substructure search, fingerprints, 2D
descriptors, standardize, pKa) loads every conformer of every row it reads, in the worker
processes:

| on the 1k × 300 set | 2D column | 3D column |
|---|---|---|
| substructure screen fingerprints | 0.48 s | 3.0 s |
| substructure search (amide) | 14 ms | 416 ms |
| MolWt + TPSA | 0.07 s | 1.7 s |

Shipping the blobs is not the cost (pickling 134 MB takes 0.2 s); loading them is: 1.2 ms a
molecule with 190 conformers (1.7 ms with their energies) against 0.1 ms for one conformer.
These run as background jobs, so the UI stays responsive, but they take 6–30 times longer.
Two ways to close it, for a decision:

1. **Give 2D-only work the cell's SMILES.** A 3D cell's display text is its canonical SMILES
   (with stereo), which parses in about 0.3 ms. A `ToolContext` helper would hand those tools
   the text instead of the blob when a cell has many conformers. Contained, but every 2D tool
   that reads a 3D column has to be moved to it, and charges, isotopes and stereo must read
   the same from the text as from the 3D structure (tests for each).
2. **Keep the conformers outside the cell** (the storage change ROADMAP §7.0f names): cells
   hold the first conformer, and a per-column store keeps the others as float32 arrays, read
   lazily and saved as their own document part. It fixes this at the root and halves the
   memory again, but touches every tool that reads conformers (Superpose, Similarity 3D, the
   Browser, export, the conformer energies) and the document format.

The numbers so far do not call for (2): at 1k × 300 everything the user waits on is within
target. (1) is the smaller step if the 2D jobs on 3D columns matter in practice.
