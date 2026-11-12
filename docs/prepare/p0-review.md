# Protein preparation P0: the benchmark and the baselines (2026-11-12)

Plan: [`docs/PREPARE_PLAN.md`](../PREPARE_PLAN.md). Code: `scripts/review/prepare_truth.py`
(ground truth and scoring, tested in `tests/protein/prepare/test_prepare_truth.py`) and
`scripts/review/prepare_review.py` (inputs, baselines, scores). Data and tools sit under
`samples/prepare/`, which is not committed.

## 1. Licences

| Material | Use | Licence evidence | Verdict |
|---|---|---|---|
| PDB entries (neutron set) | Ground truth, local only | wwPDB data are CC0 1.0 | **OK** |
| CCD components | Ideal hydrogens of each ligand | wwPDB, CC0 1.0 | **OK** (already used by the viewer) |
| PROPKA 3.5.1 | Baseline, in a separate venv | LGPL-2.1 (repository LICENSE; Debian, Nixpkgs agree) | **OK** as a dependency too (already decided) |
| Reduce 4.17, `rlabduke/reduce` commit `6749e64` | Baseline, built in WSL | Apache-2.0 (GitHub licence API: LICENSE, NOTICE) | **OK**; its rules could be reused with citation and NOTICE |
| CCP4 monomer library, `MonomerLibrary/monomers` commit `0af022b` | `gemmi h` baseline | LGPL-3.0 | **OK** for the benchmark; shipping it is a separate question (P1) |
| gemmi 0.7.5 (`gemmi h`) | Baseline | MPL-2.0 (already a dependency) | **OK** |
| Protoss via ProteinsPlus REST | Baseline, public PDB IDs only (decided 2026-11-11) | The site states no licence or terms beyond a 30-jobs-a-minute rate limit | **OK** for benchmarking public entries; outputs kept locally, not redistributed |
| MOE Protonate3D | Optional baseline, not run | Commercial licence held by the user | Outputs kept locally if ever run |
| PKAD-2 (Clemson) | pKa benchmark (P2) | No licence stated; the page asks for citation (Ancona, Bastola & Alexov, *J. Comput. Biophys. Chem.* 2023) | Benchmark use only, not redistributed; cite it |

## 2. The neutron set

- RCSB search, `exptl.method` = `NEUTRON DIFFRACTION`, on 2026-11-11: **262 entries** (143
  joint X-ray/neutron, 28 more neutron + X-ray, 90 neutron only, 1 with NMR). 124 hold a
  ligand larger than 100 Da other than common additives. Among them are the plan's named cases:
  4PDJ (DHFR–DHF–NADP⁺), 7D6G (DHFR–folate–NADP⁺), 4G0C (carbonic anhydrase II–acetazolamide),
  2ZYE (HIV-1 protease–KNI-272), 5MO0/5MOQ (trypsin–benzamidine), 4JEC/5T8H (HIV protease).
  The list is `samples/prepare/neutron/ids.txt`; downloading it again from the PDB recreates
  the set.
- **What counts as truth** (`prepare_truth.truth`):
  - Exchangeable sites are usually modelled as an H/D pair whose occupancies add up. A polar
    atom's hydrogen count is the summed occupancy of the H and D atoms within 1.25 Å; a proton
    counts at ≥ 0.5.
  - A hydrogen nobody modelled is not the same as one refined to be absent. Lys, Tyr and Cys
    count only when every site was modelled (three on NZ, one on OH or SG); a Cys in a
    disulfide or on a metal is left out. His needs a proton on at least one ring N.
  - Ligand atoms: an atom the CCD draws with hydrogens counts only if the depositors modelled
    some; ligands not in the CCD are left out.
  - Orientations: an Asn/Gln counts only where its amide N carries two hydrogens (so the
    neutron data showed which atom is the N); His orientations always count.
- **Inputs** (`prepare_truth.make_input`): heavy atoms only, first conformation, waters as
  HOH, and **half the Asn/Gln/His side chains flipped at random** (a sequence per entry), so a
  method that never flips scores 50% on orientation. The first run shared one random sequence
  across all entries, which flipped 62% of the scored sites; fixed before these numbers.
- Every structure is scored at its crystallisation pH (`_exptl_crystal_grow.pH`, else the pH
  in the details text, else 7.0). Neutron crystals are often in D₂O; the deposited value is
  used as it is.

## 3. Scores

`python scripts/review/prepare_review.py score`, 2026-11-12. Each cell is right/scored; the
number after a method is how many of the 262 entries it answered. "–" means the method gives
no answer to that question:
- PROPKA places no hydrogens;
- Protoss started from the deposited entry (public IDs only), whose orientations are the
  answer;
- naive and PROPKA leave ligands as the CCD draws them.

MOE Protonate3D is not scored: its outputs are hard to get, so the comparison goes ahead
without it (decided 2026-11-12).

| Question | naive (262) | reduce (262) | propka (262) | gemmi (256) | protoss (259) |
|---|---|---|---|---|---|
| ASN orientation | 50% (1496/3011) | **85%** (2571/3011) | – | 50% (1472/2961) | – |
| GLN orientation | 48% (1020/2109) | **79%** (1664/2109) | – | 48% (996/2067) | – |
| HIS orientation | 51% (903/1774) | **88%** (1567/1774) | – | 51% (881/1734) | – |
| His charge (HIP or neutral) | 71% (1176/1653) | **74%** (1231/1653) | 71% (1102/1549) | 39% (638/1617) | 73% (1191/1642) |
| His tautomer (HID/HIE/HIP) | 45% (742/1653) | **62%** (1033/1653) | 50% (132/266) | 39% (638/1617) | 59% (969/1642) |
| ASP state | **99%** (3916/3956) | **99%** (3916/3956) | 95% (3626/3798) | **99%** (3856/3896) | 98% (3830/3896) |
| ASP uncommon state found | 0% (0/40) | 0% (0/40) | 35% (14/40) | 0% (0/40) | **40%** (16/40) |
| GLU state | **98%** (3737/3810) | **98%** (3737/3810) | 92% (3259/3537) | **98%** (3674/3747) | **98%** (3674/3766) |
| GLU uncommon state found | 0% (0/73) | 0% (0/73) | **42%** (26/62) | 0% (0/73) | 12% (9/73) |
| LYS state | 98% (3114/3164) | 98% (3114/3164) | 98% (2929/2993) | 98% (3055/3105) | 98% (3053/3120) |
| LYS uncommon state found | 0% (0/50) | 0% (0/50) | 0% (0/50) | 0% (0/50) | 0% (0/50) |
| TYR state | 99% (2285/2312) | 99% (2285/2312) | 99% (2167/2196) | 99% (2246/2279) | 99% (2252/2279) |
| TYR uncommon state found | 0% (0/27) | 0% (0/27) | 0% (0/27) | 0% (0/27) | 0% (0/27) |
| CYS state | 96% (340/356) | 96% (340/356) | 95% (308/325) | 95% (323/340) | 95% (331/347) |
| CYS uncommon state found | 0% (0/16) | 0% (0/16) | 0% (0/16) | 0% (0/16) | 0% (0/16) |
| Ligand polar-atom H count | **96%** (1471/1528) | 91% (1383/1519) | – | 92% (1356/1468) | 94% (1374/1462) |
| Ligand atoms protonated unlike the CCD | 0% (0/57) | **37%** (21/57) | – | 6% (3/52) | 31% (15/49) |

The Cys rows were re-scored on 2026-11-12 (during P1): a hydrogen counted as bonded within
1.25 Å, which missed S–H bonds at their nuclear length (1.34 Å), in the deposits and the
outputs alike. S and Se now allow 1.5 and 1.6 Å (`prepare_truth.bond_limit`), and 356 Cys are
scored instead of 218. No other row changed.

The "uncommon" questions are the neutron sites that are not in the textbook state at the
crystal's pH: a protonated Asp or Glu, a neutral Lys, a Tyr or Cys anion. The row above each
one scores every site of that residue type.

`tests/protein/test_bench_prepare.py` (`pytest -m bench`) checks the harness: the method that
flips nothing must land at 45–55% on orientation, and Reduce must keep these scores.

## 4. What the numbers say

1. **Flips are solved well only by Reduce.** It gets 85% of Asn, 79% of Gln and 88% of His
   orientations right. gemmi never flips, so it scores the 50% that the random flips give.
   P1 must at least match Reduce here; its rules (Apache-2.0) can be cited and reused.
2. **His is the hardest residue and nobody is good at it.**
   - On tautomer, the best is 62% (Reduce), then Protoss at 59%. Calling every His HIE
     scores 45%.
   - On charge, every method sits at 71–74%, barely above calling every His neutral (71%).
     So about one His in four is charged or neutral against expectation. A desolvation term
     and the H-bond network together should do better, and this is where the coupled model
     has to earn its place.
   - PROPKA gives a tautomer for only 266 of 1,653 His and gets half of those right.
   - gemmi's 39% on both questions comes from one fixed His form.
3. **Uncommon states: PROPKA finds some but pays for them.**
   - The neutron data hold 40 protonated Asp, 73 protonated Glu, 50 neutral Lys, 27 Tyr
     anions and 16 Cys anions.
   - PROPKA finds 35% of the Asp and 42% of the Glu, but its wrong calls cost more than that.
     Asp drops from 99% to 95% overall (172 wrong of 3,798, against 40 for the textbook
     state) and Glu from 98% to 92%.
   - Protoss finds 40% of the Asp at nearly no cost (98% overall), but only 12% of the Glu.
   - No method finds a neutral Lys or an anionic Tyr or Cys.

   The target for P2 is PROPKA's or Protoss's hit rate without the overall state falling
   below the textbook 99%/98%.
4. **Ligands: the CCD's own protonation is the score to beat.** Ligands as the CCD draws them
   get 96% of polar-atom hydrogen counts right. Reduce (91%), gemmi (92%) and Protoss (94%)
   are all worse overall, though Reduce and Protoss catch 37% and 31% of the 57 atoms whose
   protonation differs from the CCD. Two caveats:
   - the set is small (57 exceptions in 124 ligand-bearing entries);
   - the CCD form is often the neutral one, and PDB depositors tend to model the form they
     expect.

   P3's coupled model has to keep the 96% while catching the exceptions; MolGpKa in the
   environment is the tool for that.
5. **Coverage matters as much as accuracy.** `gemmi h` fails on 6 entries (missing monomers)
   and Protoss on 3. Reduce and the trivial methods answer everything. Our pipeline must never
   stop on an unknown ligand: a ligand without a template gets RDKit hydrogens and a ⚠
   instead.

**Acceptance thresholds** carried into P1–P3, from the best baseline on each question:

| Question | Threshold |
|---|---|
| Asn / Gln / His orientation | ≥ 85% / 79% / 88% |
| His tautomer | ≥ 62% |
| His charge | ≥ 74% |
| Asp and Glu state overall | ≥ 99% and 98%, finding ≥ 40% of the uncommon ones |
| Ligand polar-atom H count | ≥ 96%, catching ≥ 37% of the atoms protonated unlike the CCD |

MOE Protonate3D is optional. If its outputs (`samples/prepare/moe/entries.csv`) are ever run,
it joins the table, and where it beats these thresholds its number becomes the threshold.

## 5. Found on the way

- `Component.to_mol` raised on the CCD's element "X" (an atom of unknown element) instead of
  returning None. Fixed in `chemsuite/protein/structure/components.py`, with a test.
- `gemmi h` stops on any ligand missing from the monomer library (6 of 262 entries). This is
  what the plan feared: gemmi's hydrogens depend on a restraint dictionary for every ligand.
- Protoss returns an empty ligand SDF for some entries; the reader skips those ligands.
