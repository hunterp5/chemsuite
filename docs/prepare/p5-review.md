# Protein preparation: can titration decide states? (2026-11-19)

Plan: ROADMAP §7a 3 and [`p2-review.md`](p2-review.md) §4. Code: `chemsuite/macromolecule/prepare/`
(`energy`, `sites`, `pipeline`). Benchmarks: `scripts/review/prepare_states.py` (the neutron
set) and `scripts/review/pkad_features.py` (PKAD-R).

**Outcome:** no. Retargeted at the states preparation writes, titration still does not beat
titration off (P1) on the neutron set's test half (0.884 against 0.893), so it **stays off by
default**. What came out of the round:
- a benchmark of states, not pKa values;
- a His that keeps P1's hard clash when titration is on;
- a decision rule that leaves the reported pKa values alone.

Three options are kept for the next attempt (§5).

## 1. pKa error was the wrong measure

Preparation for docking and MD needs each residue in the right form at the pH, not its pKa.
On PKAD-R's training half (286 residues, metal-bound ones left out), at pH 7:

| | States right | Unusual states caught (of 21) | Wrong switches |
|---|---|---|---|
| Null (every group in its form in water) | 265 | 0 | 0 |
| PROPKA 3.5.1 | 267 | 2 | 0 |
| chemsuite, P2 titration | 260 | 2 | 7 |

- **93% of residues are in the usual form.** Of the 21 that are not, 14 are His charged at
  pH 7, 4 are Cys thiolates and 3 are Asp/Glu.
- **A margin does not help.** Requiring the prediction to clear pH 7 by 0.5 units removes
  chemsuite's 7 wrong switches, and its 2 catches with them.
- **PKAD-R is too small to tune on**, so the round moved to the neutron set (§3).

## 2. What carries signal (PKAD-R, cross-validated)

`pkad_features.py` gives each residue heavy-atom features, PROPKA's kinds of determinant:
- backbone and side-chain H-bond partners;
- charged partners;
- the buried fraction;
- Coulomb from the charged groups within 12 Å;
- a helix's first turn near an anion.

It then fits the shift from the model pKa by ridge regression per type, leaving one family
out at a time, on the training half only:

| Model of the shift | RMSE, all types | Without Cys |
|---|---|---|
| None (the null) | 1.145 | 0.996 |
| An offset per type | 1.085 | 1.006 |
| The features | 1.055 | 0.952 |
| The features and an offset | 1.039 | 0.949 |

- **The signal is real but small.** Against the fair null (offsets fitted), the features
  gain 0.046.
- **Which features matter:**
  - Asp: backbone and side-chain H-bonds (r ≈ −0.33 each);
  - Glu and Lys: burial (r = +0.41 and −0.47);
  - His: Coulomb (r = −0.37);
  - Cys: backbone N–H (r = −0.65, on 11 residues).
- **Charged partners and the helix dipole carry almost none.** That is why P2's `ionic` term
  fit to 0.

## 3. The neutron set, scored on states

`prepare_states.py build` makes each entry's energy tables once (titration on, and off) at
the crystal pH, with each site's observed state, and pickles them.

- `score` decides states as `pipeline.prepare` does: the best states, then each titrating
  group's form from its populations at that pH.
- `fit` searches weights and thresholds.
- `--fast` reads states off the mean-field populations. It agrees with the exact decision on
  every residue type, and on His within about a point, and takes 4 s a pass instead of 9 min.

Questions per type:
- the share of states right;
- the unusual states found (ASH, GLH, LYN, TYM, CYM, HIP);
- the wrong switches;
- His charge, tautomer and orientation.

The objective is the mean of each type's share right plus the three His questions, so a wrong
switch costs as much as a find gains. The training and test halves are those of
`prepare_review.py` (CRC32 of the PDB ID).

**Found on the way:** an N-terminal Lys or Cys has two sites with one name and number, the
backbone NH₃⁺ and the side chain. The harness scored both against the side chain's state,
which made 16 Lys and 2 Cys "wrong switches" that no weight could fix. N-termini are now left
out.

### What was tried, on the training half (fast)

| | Objective |
|---|---|
| Null | 0.826 |
| P1 (titration off) | 0.901 |
| P2 titration, as built | 0.865 |
| P2 titration, its own weights refit | 0.894 |
| Plus H-bond determinants and a per-type threshold | 0.901 (a tie with P1) |

- **P2 as built** switches 77 Glu, 22 Asp and 16 Lys that are not unusual, and charges 187 His
  that are neutral.
- **The refit** sets `bump` to 32 (a near-hard clash on HIP), which brings His orientation
  and tautomer back to P1's.
- **The per-type threshold** is the share the other form needs before a group leaves its form
  in water. The fit took Asp, Glu and His to 0.95 and Lys and Tyr to 0.99, so the model has to
  be nearly certain. Even then, Glu found 11 of 39 at the cost of 19 wrong switches, and Asp
  found none.
- **First tried as an energy term** (a bar on the unusual form): it shifts every reported pKa
  (each Asp would read about 0.9) and holds Asp charged even at pH 3, so it was moved to the
  decision instead.
- **The H-bond determinants** (fitted weights 1–2) changed the objective by less than 0.001,
  so they were not kept.

### The test half, scored once

| Question | Null | P1 | Titration, fitted |
|---|---|---|---|
| Asp state | 99.3% | 99.3% | 98.9% (0/13 found, 8 wrong) |
| Glu state | 98.2% | 98.2% | 97.2% (5/34 found, 24 wrong) |
| Lys state | 98.4% | 98.4% | 98.4% |
| Tyr state | 99.3% | 99.3% | 99.3% |
| Cys state | 95.9% | 95.9% | 95.9% |
| His charge | 68.7% | **72.6%** | 69.5% |
| His tautomer | 43.2% | **63.3%** | 60.9% |
| His orientation | 50.3% | **87.6%** | 87.4% |
| Objective | 0.817 | **0.893** | 0.884 |

None of the round's gates is met:
- Asp and Glu finding 40% of their unusual states;
- His charge above 74%;
- every type at the null's share or better.

P1 stays the default. Titration stays an option (Titrate: Pocket or All) for pKa values and
for a look at the groups the model calls mixed.

## 4. What changed in the code

- **A titrating His keeps P1's hard clash** on its neutral states, which decide its flip and
  tautomer. Only HIP's contacts move to `bump` and `repel` (`energy.soft_states`). Before,
  titration on softened the whole His, and His orientation and tautomer dropped.
- **`pipeline.protonated` decides a protein group's form** at the pH: its form in water there,
  unless the other form's share reaches the group's threshold in `pipeline.SWITCH`.
  - The populations and pKa values are untouched, and the rule is the same at any pH.
  - `SWITCH` stays a plain majority (0.5): the fitted thresholds did not help.
  - `sites.group(site)` names the group.
- **Titration off is unchanged**: P1 scores the same as before on every question. `score`
  asserts it, since a weight shared with P1 (as `ionic` is) would change it.

## 5. Kept for the next attempt

Nothing tried here finds the unusual states without paying more in wrong switches. The next
step is to **run existing implementations through `prepare_states.py`, as benchmark-only tools
(the way PROPKA and Protoss are), before building any**:

1. **Poisson–Boltzmann titration** (the H++, PypKa and MCCE family): continuum electrostatics,
   in place of this model's Coulomb and desolvation terms.
2. **An ML protein-pKa model**, trained in WSL and run in numpy as chemsuite pKa is. Measured
   protein pKa values are few (about 1–2k), so these mostly learn from computed ones.
3. **Constant-pH MD for the sites the report marks ⚠**, in the OpenMM worker. It decides
   states by simulation, not prediction, for the few that matter in an MD run. OpenMM's
   support is still young.

Licences are checked before any of them is used, even as a benchmark.

Also open:
- a pocket-only tally, since for docking only the pocket's states matter (the cases would
  need ligand positions);
- titrating the C-terminus (always charged at pH 7; not measured here).
