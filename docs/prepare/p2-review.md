# Protein preparation P2: titration (2026-11-12)

Plan: [`docs/PREPARE_PLAN.md`](../PREPARE_PLAN.md); P1: [`p1-review.md`](p1-review.md). Code:
`chemsuite/protein/prepare/` (`sites`, `energy`, `populations`, `solve`, `pipeline`). The
benchmark is `scripts/review/pkad_review.py`.

**Outcome:** P2 works end to end and reaches PROPKA's accuracy on experimental pKa values. It
does not beat the null model, and with titration on the neutron marks drop. So titration is
**off by default** (`pipeline.Options(titrate=True)` turns it on). P1 is unchanged: with
titration off, all 262 neutron outputs are byte-identical to P1's.

## 1. What runs

- **Forms** (`sites.find_sites(p, titrate=True)`). Each form and its rotor or flip variants are
  states of one site:

  | Residue | Forms |
  |---|---|
  | Asp, Glu | charged, or neutral with H on either O, syn or anti |
  | His | HID, HIE, HIP, each as built or flipped |
  | Lys, N-terminus | NH₃⁺ (8 rotor states), or NH₂ (lone pair in 15° steps) |
  | Cys (free) | SH (24 rotor states), or S⁻ |
  | Tyr | OH (2 in plane), or O⁻ |

  Not titrated: Arg; a Cys in a disulfide or on a metal; a His bonded to another residue; the
  C-terminus (PKAD-R has no termini to test it on); ligands (P3).
- **Free energies.** Each state carries its free energy in water from PROPKA 3's model pKa
  values: Asp 3.8, Glu 4.5, His 6.5, Cys 9.0, Tyr 10.0, Lys 10.5, N-terminus 8.0. It also
  carries its proton count. The tables are built once; any pH adds 1.364 · protons · pH.
  - Microstate corrections make an isolated group titrate at its model pKa: the two oxygens
    of an acid, HIE:HID 4:1, and kT·ln(variants) for rotor states. The tests check this for
    every residue type.
- **Environment terms**, linear in fitted weights as in P1:
  - Coulomb on per-state charges (ε = 4r), plus a burial-scaled copy;
  - PROPKA's buried fraction as acid and base desolvation;
  - charged H-bonds;
  - a titrating site's clash and lone-pair contacts, moved to their own `bump` and `repel`
    weights. P1's hard clash (160 × overlap²) turned crystal contacts into 20–370 kcal/mol
    against HIP or NH₃⁺.

  For a titrating group, the water term counts every N, O and S, since its forms all have
  their water in the model pKa.
- **Populations** (`populations.titration`). Mean field over a pH grid (0–14, step 0.25), each
  site updated only while it or its neighbours move. Metropolis Monte Carlo is the check.
  - On T4 lysozyme the two agree within 0.1 pH unit; mean field takes 0.1 s, Monte Carlo
    6.6 s. All 104 PKAD-R structures take 12 s.
  - The tests also check both against an exact Boltzmann sum.
- **Pipeline.** With `titrate=True`, `prepare` writes each group in its majority form (the
  best state within it) and returns `Prepared.pka` and `Prepared.population`. The report:
  - names each uncommon form (ASH, GLH, LYN, TYM, CYM, NH2, HIP) with its protonated share;
  - marks ⚠ a group whose forms are mixed (between 10% and 90%).

## 2. PKAD-R

- **Data.** PKAD-R (Chen et al., *J. Comput. Biophys. Chem.* 2025,
  doi:10.1142/S2737416525500164), saved by hand. Its main wild-type values with no warning:
  581 values in 104 structures.
- **Split.** Training and test halves by family, proteins over 30% identical kept together
  (80 families).
- **PROPKA.** 3.5.1 on the same PDB files. A Cys in a disulfide (PROPKA's 99.99) is no answer,
  so 553 residues are answered by every method.

`pkad_review.py fit`, a coordinate search and then Nelder–Mead, chose the defaults now in
`energy.Weights`:

| Weight | Value |
|---|---|
| bump | 1 |
| repel | 0 |
| screened | 0.02 |
| desolv_acid | 1 |
| anti | 4 |
| deep, ionic, desolv_base | 0 |

RMSE in pH units (share within 1 unit):

| Residues | Half | Null model | PROPKA | chemsuite |
|---|---|---|---|---|
| all | training | 1.06 (80%) | 1.06 (81%) | 1.09 (74%) |
| all | **test** | **1.15** (72%) | 1.33 (68%) | 1.32 (66%) |
| Asp | test | 1.21 | 1.09 | 1.28 |
| Glu | test | 0.79 | 1.03 | 0.94 |
| His | test | 1.27 | 1.38 | 1.59 |
| Cys | test | 3.18 | 4.49 | 3.67 |
| Tyr | test | 0.82 | 0.87 | 0.39 |
| Lys | test | 0.66 | 0.63 | 0.77 |

- **Matching PROPKA** (test half): met, 1.32 against 1.33. `tests/protein/test_bench_prepare.py`
  holds it.
- **Beating the null model**: not met, by chemsuite or PROPKA. The fit drives every
  environment weight towards zero: as built, the terms add as much noise as signal.
- **The worst errors:**
  - catalytic Cys thiolates (measured 4.7–5.6, predicted about 9.5): nothing stabilises the
    anion;
  - His on a metal in the crystal, measured without it;
  - Asp and Glu with strongly shifted pKa values.

## 3. Neutron set with titration on

`prepare_review.py chemsuite --titrate`, then `score` (the `titrated` column):

| Question | Titration off (P1) | Titration on | P2 mark |
|---|---|---|---|
| Asp state overall | 99.0% | 98.1% | ≥ 99% ✗ |
| Glu state overall | 98.1% | 94.7% | ≥ 98% ✗ |
| Uncommon Asp found | 0% (0/40) | 12% (5/40) | ≥ 40% ✗ |
| Uncommon Glu found | 0% (0/73) | 49% (36/73) | ≥ 40% ✓ |
| His orientation | 89.6% | 87.1% | keep P1 ✗ |
| His tautomer | 63.6% | 57.2% | keep P1 ✗ |
| His charge | 73% | 66% | > 74% ✗ |
| Asn / Gln orientation | 86.2 / 80.4% | 86.0 / 80.3% | keep P1 |

His suffers most: with titration on it loses P1's fitted priors and hard clash. Hence the
default.

## 4. Open (the next round, ROADMAP §7a 3)

- PROPKA-style determinants per residue type:
  - side-chain and backbone H-bonds to a carboxylate, a thiolate or an imidazolium;
  - the helix dipole for Cys.

  Then refit, aiming to beat the null model and to turn titration on by default.
- His with titration on: keep P1's orientation terms for the flip, titrate only the charge.
- Titrate the C-terminus.
- Bench and runs: `pytest -m bench tests/protein/test_bench_prepare.py -s`;
  `pkad_review.py tables`, `fit`, `chemsuite`, `score`.
