# Protein preparation: existing pKa tools, licenses first (2026-11-19)

Follows [`p5-review.md`](p5-review.md) §5: before building Poisson–Boltzmann or ML titration,
run existing implementations through `scripts/review/prepare_states.py` as benchmark-only tools
(the way PROPKA and Protoss are). This page is the survey: what exists, its license, what it
outputs, and the order to try them in. No code yet.

## 1. The bar, stated first

Every candidate below predicts **per-site pKa values** (PypKa and MCCE also give
populations). None decides His tautomer or flip. Those stay with P1's hydrogen-bond optimizer,
which is where P1 already wins (tautomer 63.3%, orientation 87.6%). So a tool can only move the
states benchmark through `pipeline.protonated`'s charge decision: Asp/Glu protonated, Lys
neutral, Tyr/Cys anions, His charged.

- **That is a narrow target.** On PKAD-R at pH 7, 93% of residues are in the usual form. PROPKA
  3.5.1 caught 2 of 21 unusual states and made no wrong switches (p5 §1).
- **PROPKA has never been scored on the neutron set's states.** If PROPKA, at about 1.1 RMSE
  on PKAD, does not beat P1 there, a tool at 0.8–1.0 RMSE needs to be much better on the
  *unusual* states specifically, not on average.
- **What makes a candidate worth its cost:** recall of the unusual states without wrong
  switches. KaML is the only one that reports this directly (§3).

## 2. The candidates

| Tool | Method | License | Use here | Runs on | Outputs |
|---|---|---|---|---|---|
| **KaML-CBtree** (Shen lab, [KaMLs](https://github.com/JanaShenLab/KaMLs)) | CatBoost trees on 37 structure features | MIT, models included | **Benchmark; could ship** (trees in numpy) | Linux (WSL); CPU; DSSP and RIDA binaries | pKa per Asp, Glu, His, Cys, Tyr, Lys |
| **pKAI+** (Bayer, [pKAI](https://github.com/bayer-science-for-a-better-life/pKAI)) | MLP on N/O/S atoms within 15 Å | MIT, TorchScript weights in the repo | **Benchmark; could ship** (MLP in numpy) | Any; torch, CPU | pKa per Asp, Glu, His, Lys, Cys, Tyr, termini |
| **MCCE4** with **NextGenPB** ([MCCE4-Alpha](https://github.com/GunnerLab/MCCE4-Alpha), [NextGenPB](https://github.com/concept-lab/NextGenPB)) | PB + Monte Carlo over conformers | MIT + GPL-3.0 | Benchmark only (separate interpreter, like OpenMM) | Linux/WSL via Apptainer; Windows not supported natively | pKa, occupancy of each conformer (charge and tautomer) at a pH |
| **PypKa** ([PypKa](https://github.com/mms-fcul/PypKa)) | PB + Monte Carlo, with tautomers | LGPL-3.0, **but needs DelPhi** (proprietary; free to registered non-profit labs) | Benchmark only, **if** DelPhi's terms allow it | Linux | pKa, tautomer populations, protonated PDB at a pH |
| **openmm-cph** ([openmm-cph](https://github.com/openmm/openmm-cph)) | Constant-pH MD (explicit solvent, implicit-solvent MC moves) | MIT | Later, in the OpenMM worker (p5 §5.3) | Wherever OpenMM runs | States sampled by simulation |
| KaML-ESM ([kaml-esm](https://github.com/JanaShenLab/kaml-esm)) | ESM embeddings + MLP | **CC BY-NC 4.0**; ESM-C through a paid API | Excluded | 24 GB GPU, 128 GB RAM | pKa |
| pKa-ANI ([pKa-ANI](https://github.com/isayevlab/pKa-ANI)) | ANI representations | **Non-commercial, no sharing; the licensor owns modifications** | Excluded, even as a benchmark | — | — |
| PDB2PKA | PB (APBS) | BSD | Excluded: removed from PDB2PQR 3.x, which is itself dropped | — | — |
| ME-pKa ([ME-pKa](https://github.com/yzjyg215/ME-pKa)) | ESM + structure, multimodal | MIT | Not now: no published weights; needs CHARMM and Chimera to featurize | CUDA | — |

Notes:
- **openmm-cph** is still prerelease ("not yet well validated", no releases, last commit
  2026-08-14) and not in OpenMM itself. It supports `ForceField`-built systems only, which
  matches the worker. It stays the MD-side option and is not a preparation step.
- **H++, Karlsberg+, DeepKa and the DelPhiPKa web server** are web services. Sending structures
  to them is publishing them and they cannot be scripted at benchmark scale, so they are left
  out.

## 3. What each claims

- **KaML-CBtree** (Shen & Shen, 2024–2025, [PMC11601431](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601431/)),
  trained on **PKAD-3** (1,167 experimental values, 247 proteins):
  - RMSE: acids 0.76 against PROPKA3's 1.28, bases 0.79 against 0.95. Per type: Asp 0.75,
    Glu 0.60, His 0.85, Cys 1.50, Tyr 1.24, Lys 0.70.
  - **States:** protonated Asp recall 0.77 (PROPKA3 0.28); protonated Glu precision 1.0
    (0.67); Cys thiolate precision/recall 0.73/0.76 and neutral Lys 1.0/0.8 (PROPKA3 finds
    neither).
  - **Caveat:** the split is random by residue, not by protein, so these numbers are
    optimistic. PKAD-3 contains PKAD-R's residues, so PKAD-R cannot test it. The neutron set
    is independent of its training data and is the fair test.
  - Input must be a single chain without missing atoms. Our prepared models are complete, so
    the bench can feed it one chain at a time (losing inter-chain contacts).
- **pKAI+** (Reis et al., JCTC 2022, [PMC9369009](https://pmc.ncbi.nlm.nih.gov/articles/PMC9369009)):
  - Trained on about 3 million PypKa-computed pKa values, not on experiment; pKAI+ adds a
    penalty on large shifts. RMSE 0.98 on 736 PKAD residues (pKAI 1.15).
  - Heavy atoms only (N, O, S within 15 Å), no hydrogens, 4,008 inputs, three hidden layers.
  - Its shrinkage toward the model pKa is what made it accurate on average, and it is also
    what makes it unlikely to call unusual states. The bench will show which wins.
- **MCCE4 and PypKa** are the "Poisson–Boltzmann titration" family that p5 §5.1 names. They
  replace exactly the Coulomb and desolvation terms P2/P5 found weak. Neither publishes run times on
  the pages read; PB with Monte Carlo is expected to take minutes per protein (MCCE with
  side-chain conformers longer), to be measured. That is fine for a benchmark of about 100
  structures run once, but not interactive.

## 4. Order of work (proposed)

1. **One hook in the harness.** Add a method to `prepare_states.py` that takes per-site pKa
   values from a file (`<ID>.csv`: chain, residue, name, pKa) and scores them through the same
   decision as P1's titration (`pipeline.protonated` at the crystal's pH, tautomer and flip
   from P1). One hook serves every tool below.
2. **PROPKA through it first.** It is already a bench dependency, and it sets the bar.
3. **KaML-CBtree and pKAI+** in one WSL bench venv (catboost, pycaret, torch CPU, DSSP, RIDA;
   nothing reaches the app's environment).
4. **MCCE4 + NextGenPB** in WSL through its Apptainer setup, if step 3 shows that pKa quality
   moves the states score at all.
5. **PypKa** only if you decide DelPhi's registration covers this use (your call: it is free
   to non-profit labs that register with Columbia).

If a tool beats P1 on the test half, shipping it is a separate plan:
- **KaML-CBtree:** CatBoost exports oblivious trees to JSON, which evaluate in numpy in a few
  lines. The work is reimplementing its 37 features Qt-free (DSSP-like secondary structure and
  the disorder flag included). MIT models can ship with a NOTICE, as chemsuite pKa's do.
- **pKAI+:** an MLP, so a numpy port is direct (the chemsuite pKa route; never torch in the
  app). Its featurization is simple.

Data to check before any training of our own (not proposed now): PKAD-3's download terms. Its
site did not load on 2026-11-19, and the paper names no data license.

## 5. Steps 1 and 2 done: the hook, and PROPKA as the bar (2026-11-19)

**The hook.** `prepare_states.py score --tools a,b` adds a column per tool, read from
`samples/prepare/out/pka/<tool>/<ID>.csv` (chain, res_seq, icode, group, pka; group a key of
`sites.MODEL_PKA`, so `NTERM` for an N-terminus).
- Each value gives its group's protonated share at the crystal pH (Henderson–Hasselbalch).
  `pipeline.protonated` turns the share into a form, as titration does.
- A group the file leaves out keeps its form in water at that pH.
- The titrated tables are then solved with those forms fixed, so His tautomer and flip are
  chosen as titration on chooses them.
- Only entries every listed tool covers are scored.
- The cases now carry each site's key (`Case.keys`), so `build --redo` was needed once.
- `prepare_states.py propka` turns `prepare_review.py propka`'s `.pka` files into these CSVs
  (262 entries).

A tool with no values does not reproduce the null column, by design. The null holds every
group in its pH 7 form, while the hook (like the pipeline) uses the form at the crystal pH: at
pH 5.5 every His is charged. Checked on all 15,243 scored sites: each lands in its form in
water at the crystal pH, except the 150 His in entries at pH 6.5 (His's model pKa), where the
tie goes to the charged form, as it does in the pipeline.

**PROPKA 3.5.1 on the test half** (124 entries, default weights, `score --test --tools propka`):

| Question | Null | P1 | Titration | PROPKA |
|---|---|---|---|---|
| Asp state | 99.3% | 99.3% | 98.4% | 95.7% |
| Asp unusual found | 0/13 | 0/13 | 1/13 | 4/13 |
| Asp wrong switches | 0 | 0 | 19 | 73 |
| Glu state | 98.2% | 98.2% | 94.3% | 91.5% |
| Glu unusual found | 0/34 | 0/34 | 17/34 | 14/34 |
| Glu wrong switches | 0 | 0 | 90 | 140 |
| Lys, Tyr, Cys unusual found | 0/40 | 0/40 | 0/40 | 0/40 |
| His charge | 68.7% | **72.6%** | 63.5% | 66.7% |
| His charged found | 0/261 | 72/261 | 169/261 | 62/261 |
| His wrongly charged | 0 | 39 | 212 | 78 |
| His tautomer | 43.2% | **63.3%** | 58.0% | 57.7% |
| His orientation | 50.3% | **87.6%** | 84.8% | 86.4% |
| Objective | 0.817 | **0.893** | 0.866 | 0.864 |

The training half agrees (P1 0.901, PROPKA 0.876).

- **PROPKA loses to P1, and slightly to our own titration.** It finds about a third of the
  unusual Asp and Glu, at about 3 wrong switches for each find.
- **Its His calls are worse than P1's on every His question.** P1 charges a His only where
  its neighbors make that the best state. PROPKA charges fewer correct ones and more wrong
  ones, and each wrong HIP also costs a tautomer.
- **The bar for the next tools is therefore P1's 0.893, not PROPKA's.** For a pKa tool to pass
  it, wrong switches must stay below about one per find. KaML-CBtree's reported precision
  (§3) is the only claim in reach of that. Lys, Tyr and Cys anions were not found by anything
  so far (40 on the test half).
- A margin (a switch only when the pKa clears the pH by some units) is the obvious knob for any
  tool. `SWITCH` already takes a threshold per group, so it can be tried without new code by
  passing thresholds to `decide_external`.
