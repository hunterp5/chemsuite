# Sites of metabolism (plan approved 2026-11-17)

chemsuite trains its own site-of-metabolism model on ATTNSOM's design and runs it in numpy, with
no torch. It scores every atom of a molecule for each of nine human P450 isoforms. CYP inhibition
follows on the same package (phase C below) once this ships.

## Why

- **Nothing in chemsuite says anything about metabolic fate.** It predicts pKa, descriptors and
  properties; "where does this get oxidized" is the next question a med chemist asks, and the one
  that decides whether a series is worth making.
- **This is not the dropped metabolite prediction.** `docs/ROADMAP.md` records "Metabolite
  prediction (BioTransformer) is dropped from the plan" and `docs/PLAN.md` lists it under "Left
  out … can come back later as tools without changing the core". That item is rule-based
  metabolite *structure* enumeration through an external Java tool. This is an ML model that
  scores atoms on a molecular graph: no external tool, no new runtime dependency, no generated
  structures. Metabolite structure generation stays out of scope.
- **ATTNSOM** (Kim et al., *Bioinformatics* 42(Supplement_1), 2026,
  doi:10.1093/bioinformatics/btag314, ISMB; preprint arXiv:2601.20891; dmis-lab/ATTNSOM, MIT) is
  the best available basis, and the only strong one that is both openly licensed and cheap to port:
  - graph only — RDKit 2D atom and bond features, no 3D conformers and no protein language model,
    so the whole network is matmuls and scatter-adds, about 8.6M parameters;
  - it beats the alternatives on the one benchmark with exact atom-resolved labels
    (AZ-ExactSOM, MCC: ATTNSOM 0.52, Chemprop 0.48, SMARTCyp 0.25, FAME 3 0.20) and matches
    XenoSite on Zaretzki top-2 accuracy (0.871 vs 0.870) while reporting a far better MCC than
    any graph baseline (0.710 vs Chemprop's 0.560);
  - it ships both datasets and **no weights**, so we train our own either way.
- **Rejected alternatives:**
  - **FAME 3** / FAME.AL (Šícho et al., JCIM 2019): "free for non-commercial and academic
    research", Java. chemsuite is used inside a company, so non-commercial terms rule it out —
    the same bar that ruled out Uni-pKa's weights (`docs/UNIPKA_PLAN.md`).
  - **GLMCyp** (JCIM 2025) and **MetaCYP** (2026): both need a Uni-Mol 3D encoder *and* ESM-2
    protein embeddings, so each is two more ports on top of this one, for a model that cannot run
    without 3D embedding per molecule. MetaCYP has no license file. ATTNSOM also reports that
    GLMCyp's published numbers came from an evaluation set that overlapped its training set, so
    they are not comparable.
  - **DeepCYP** (NAR 2026) is a web server; **DeepMetab** has no released code.
  - **SMARTCyp** and **XenoSite**: beaten, and XenoSite is a web service.

## Name and credit

It is not ATTNSOM, and the app does not call it that, since users would expect the paper's
accuracy from weights we trained ourselves on our own split. In the app it is **Sites of
Metabolism** (Chemistry ▸ Calculate), and in the docs **chemsuite SoM**. The user guide and the
docstrings cite:

- **ATTNSOM** (Kim, Lee, Chung, Park, Baek, Kang; *Bioinformatics* 2026,
  doi:10.1093/bioinformatics/btag314; MIT): the architecture, the atom and bond featurization, and
  the training recipe. The featurizer and forward pass are ports of its code, so its notice ships
  (`chemsuite/metabolism/data/LICENSE-ATTNSOM.txt`).
- **Zaretzki et al.** (JCIM 2012/2013, RS-Predictor and XenoSite): the training labels.
- **Chen et al.** (2025): the AZ-ExactSOM test set.

The weights are chemsuite's own, trained from scratch — there is no pretrained backbone, so unlike
the pKa weights the chain carries no upstream model license.

## License chain

| Part | Source | License | Use |
|---|---|---|---|
| Architecture and training code | dmis-lab/ATTNSOM (*Bioinformatics* 2026) | MIT (© 2026 DMIS Laboratory, Korea University) | Trained with, in WSL; featurizer and forward pass ported to numpy |
| Training labels | Zaretzki set: 679 substrates, 2,003 molecule–isoform instances, 9 isoforms (`cyp_dataset/*.sdf`, `PRIMARY_SOM` property) | ACS supporting information, redistributed in an MIT repo — see **Open question** below | Trained on; **never committed and never shipped** |
| Held-out test set | AZ-ExactSOM: 120 public compounds, exact atom-resolved sites from human hepatocytes (Chen et al. 2025; `cyp_dataset/az/az_120_compounds.csv`) | as published, in the MIT repo | Evaluation, and a second trained head only if the licensing question above forces it |
| Shipped weights | chemsuite's own training run | GPL-3.0-or-later (release `som-weights-1`, with a NOTICE) | Fetched on first use, pinned by SHA-256 |
| **Not used:** FAME 3, FAME.AL | nerdd.univie.ac.at | free for non-commercial and academic research only | – |
| **Not used:** MetaQSAR | Pedretti et al. | academic only | – |
| **Not used:** GLMCyp, MetaCYP, DeepCYP weights | GitHub, Hugging Face, web server | MetaCYP states none; all need Uni-Mol + ESM-2 | – |

**Open question (M0, to settle before M1 trains on it):** the Zaretzki labels began as ACS
supporting information and reach us through an MIT-licensed repository. We train on them and never
redistribute them, which is how pKaHub's excluded sources are handled in `docs/UNIPKA_PLAN.md`. If
the original terms turn out to forbid even that, the fallback is AZ-ExactSOM plus any
openly-licensed set, at a measured cost in accuracy. Record the decision here either way.

## Steps

- **M0 (in progress 2026-11-17):** this document, the license chain, the name, and the roadmap
  entry saying ML SoM is not the dropped BioTransformer item.
- **M1 step 2 reproduced, 2026-11-17.** 10-fold CV on Zaretzki, upstream's defaults, one RTX 4070,
  **16 minutes** for all ten folds (`~/bench/attnsom/runs/repro-zaretzki`):

  | iso | MCC | top-1 | top-2 | top-3 | precision | recall |
  |---|---|---|---|---|---|---|
  | 1A2 | 0.673 ± 0.081 | 0.767 | 0.875 | 0.908 | 0.638 | 0.775 |
  | 2A6 | 0.735 ± 0.089 | 0.856 | 0.905 | 0.924 | 0.711 | 0.828 |
  | 2B6 | 0.784 ± 0.066 | 0.828 | 0.914 | 0.941 | 0.753 | 0.857 |
  | 2C8 | 0.792 ± 0.057 | 0.851 | 0.916 | 0.958 | 0.723 | 0.903 |
  | 2C9 | 0.718 ± 0.108 | 0.783 | 0.908 | 0.965 | 0.665 | 0.826 |
  | 2C19 | 0.762 ± 0.045 | 0.821 | 0.913 | 0.963 | 0.708 | 0.860 |
  | 2D6 | 0.704 ± 0.099 | 0.752 | 0.852 | 0.915 | 0.671 | 0.789 |
  | 2E1 | 0.736 ± 0.091 | 0.785 | 0.875 | 0.944 | 0.732 | 0.805 |
  | 3A4 | 0.624 ± 0.061 | 0.648 | 0.781 | 0.848 | 0.597 | 0.702 |
  | **all** | **0.701 ± 0.028** | 0.761 | **0.865** | 0.917 | 0.664 | 0.789 |

  The paper reports MCC 0.710 and top-2 0.871; we get 0.701 and 0.865, both inside one fold
  standard deviation (MCC ± 0.028).

  **AZ-ExactSOM** (10-fold, `~/bench/attnsom/runs/repro-az`), against the paper's Table 4:

  | | MCC | precision | recall | top-1 | top-2 | top-3 |
  |---|---|---|---|---|---|---|
  | ours | 0.502 ± 0.096 | 0.462 | 0.645 | 0.627 | 0.777 | 0.841 |
  | paper | 0.517 | 0.521 | 0.595 | 0.639 | 0.797 | 0.872 |

  Every metric is within one fold standard deviation; we trade a little precision for recall.
  **Both reproductions pass**, so the architecture and the training recipe are sound and M1 can
  proceed to the production run. For scale: AZ-ExactSOM is 120 compounds, so its ± is wide and it
  stays a sanity check rather than the figure to quote.

- **The published numbers are not new-molecule numbers (2026-11-17).** Upstream cross-validates
  over molecule-**isoform instances**:

  ```python
  skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
  for trainval_idx, test_idx in skf.split(np.arange(len(graphs)), y_strat):
  ```

  A molecule appears once per isoform it was measured against (2,003 instances, 679 molecules,
  up to 9 each), so the same structure can be trained on for 1A2 and tested on for 3A4.
  `apply_no_leakage_to_dataloaders` masks only the auxiliary attention targets, and says so:
  *"The atom-level labels `y` are untouched."*

  Measured (`scripts/som/leakage_check.py`, upstream's exact split): **86.2% of test instances
  (1,727 of 2,003) share their molecule with the training split under a different isoform**,
  81–89% per fold. Only 12.8% of instances come from single-isoform molecules and so cannot leak.
  Isoforms per molecule: 257 molecules appear once, 121 twice, and 18 all nine times.

  Since the paper's own premise is that SoM annotations are highly similar across isoform pairs,
  this is exactly the kind of leakage that flatters the metric. **MCC 0.70 and top-2 0.87
  therefore describe "another isoform of a molecule we have seen", not "a molecule we have never
  seen" — and the latter is what a chemist typing a new structure gets.**

  **Measured both ways, 2026-11-17** (10-fold, identical settings, only the split differs;
  `~/bench/attnsom/runs/{repro,grouped}-zaretzki`):

  | split | MCC | top-1 | top-2 | top-3 | precision | recall | molecule exact |
  |---|---|---|---|---|---|---|---|
  | instance (upstream, and the paper) | 0.701 ± 0.028 | 0.761 | 0.865 | 0.917 | 0.664 | 0.789 | 0.550 |
  | **molecule-grouped (honest)** | **0.524 ± 0.061** | **0.618** | **0.743** | **0.811** | 0.536 | 0.586 | 0.377 |

  **The leakage is worth about 0.18 MCC and 0.12 top-2 accuracy.** Per isoform, grouped MCC:

  | iso | 1A2 | 2A6 | 2B6 | 2C8 | 2C9 | 2C19 | 2D6 | 2E1 | 3A4 |
  |---|---|---|---|---|---|---|---|---|---|
  | MCC | 0.499 | 0.517 | 0.548 | 0.560 | 0.517 | 0.600 | 0.516 | 0.445 | 0.518 |
  | top-2 | 0.727 | 0.741 | 0.776 | 0.733 | 0.764 | 0.772 | 0.760 | 0.723 | 0.716 |

  Grouping flattens the spread: every isoform lands near 0.5, and 3A4 is no longer the outlier
  (0.624 → 0.518 while 1A2 falls 0.673 → 0.499), which says the instance-split ranking between
  isoforms was mostly an artifact of how much each one shared molecules with the others.

  **What we can honestly claim:** for a molecule the model has never seen, the true site of
  metabolism is the top-ranked atom about **62%** of the time, in the top two about **74%**, and
  in the top three about **81%**. It is still a clearly useful tool — FAME 3 managed MCC 0.20 on
  exact atom-resolved annotations — but it is not "0.87 top-2".

  **Decided with the user, 2026-11-18: the user guide quotes the grouped numbers only**, per
  isoform, and does not print the instance-split figures at all. We do not publish a comparison
  against SMARTCyp, XenoSite or FAME either, because every published figure in this field is on
  the inflated basis and ranking against them would need re-running each one grouped. The guide
  says what the number means — the site is in the top two about three times in four for a
  molecule the model has not seen — and leaves it there.

  So: `scripts/som/apply_group_split.py` patches a working copy of upstream's `main.py` with a
  `--group_split` flag that keeps every instance of a molecule in one fold
  (`StratifiedGroupKFold` on canonical SMILES; the clone stays pristine and the flag off
  reproduces upstream exactly). The grouped run is the number the user guide quotes, and the
  instance-split number is reported beside it only to show why it is higher. The independent
  `som-benchmark` repository (WillDeSnoo) makes the same point with scaffold splits and is worth
  reading before the production run.

  Note **3A4 is the weakest isoform** (MCC 0.624, top-2 0.781) despite having the most data (475
  molecules) — it is the promiscuous one with the largest substrates. It is also the isoform
  chemists will look at first, so the user guide should not oversell it, and the per-isoform
  numbers belong in the guide rather than one headline figure.

- **M1:** train in WSL on the local GPU. `scripts/som/` holds the environment, the training
  driver, `convert_weights.py` (torch state dicts → one `.npz` keyed by PyTorch parameter name,
  folds keyed `f0/…` as the pKa weights are) and `golden.py` (upstream torch on ~40 diverse
  molecules → `tests/metabolism/som_golden.json`, plus per-layer tensors for the first two so the
  port can be debugged layer by layer). Also build the applicability-domain index: Morgan
  atom-environment (radius 2) fingerprints of every training atom, packed into the same `.npz`.

  **The run design, settled by the findings below:**
  1. Reproduce upstream's instance-split CV, to prove the recipe. **Done; it passes.**
  2. Re-run it molecule-grouped (`--group_split`), which is the number we quote.
  3. Decide the fold count. This **cannot** be answered from a k-fold run, because for any fold's
     test molecules every *other* fold's model has trained on them, so an ensemble cannot be
     scored without leakage. It needs a true external holdout: set aside ~10% of **molecules**,
     run grouped k-fold inside the remaining 90%, then compare one fold against the k-fold
     ensemble on the holdout. That answers both "how many folds" and "how big is the download"
     (M3) with one run.
  4. Train the production weights with the fold count and epoch budget those runs justify.
     **Epoch budget, measured 2026-11-18:** across the ten grouped folds the best validation
     epoch was 17, 17, 18, 19, 27, 32, 35, 37, 40, 42 — mean 28, median 32, and **none reached
     the 50-epoch cap**, so 50 is ample headroom and the inner-validation early stop is doing the
     work. For the port's development the working model is **grouped fold 1** (trained on 90% of
     molecules, properly early-stopped); it is chosen first-and-arbitrarily rather than by score,
     so the quoted 10-fold CV mean remains an honest estimate of it. Before release, retrain on
     ~all the data at a fixed 32 epochs and regenerate the golden file.
  5. Only then try the `SECONDARY_SOM` idea, and keep it only if step 2's metric improves.
  - **The `SAGPooling` risk is smaller than it looked, but there is a real decision in it**
    (read from `model.py` 2026-11-17). The FiLM graph vector is pooled two different ways
    depending on the batch:

    ```python
    if batch.max() == 0 and batch.size(0) > 1:
        g = gap(h, batch)                                  # plain mean over all atoms
    else:
        x_p, _, _, b_p, _, _ = self.sagpool(h, edge_index, edge_attr, batch)
        g = gap(x_p, b_p)                                  # mean over the kept 80%
    ```

    The mean branch fires only when **exactly one graph** is in the batch. Training and
    evaluation always batch, and upstream's own `inference.py: build_batch` batches **nine
    isoform copies of the one molecule**, so in practice everything takes the SAGPool branch.

    **Settled by measurement, 2026-11-17** (`scripts/som/batch_check.py`, fold 1 of the
    reproduction, 12 molecules × 9 isoforms):

    | comparison | largest absolute logit difference | top-1 site changed |
    |---|---|---|
    | nine-isoform batch vs each graph alone | **7.5566** | **32 of 108 (30%)** |
    | nine-isoform batch vs one 108-graph batch | 0.000008 | – |

    So the two branches are **not interchangeable**: taking the "simpler" mean path would have
    changed the top-ranked site in nearly a third of cases, away from the configuration whose
    accuracy was actually measured. The SAGPool branch, in turn, is independent of how much else
    is in the batch (8e-6 is float noise).

    **Therefore:** the port **must implement `SAGPooling`** (GraphConv score, then per-graph
    top-k at ratio 0.8 — about 30 lines of numpy) and **must always run the isoforms as one
    batch**, exactly as `predict_sites(mol, isoforms)` is designed to. The one-graph mean branch
    is dead code we deliberately never reach, and the golden file must be generated through the
    batched path so the port is checked against the right numbers. No retraining is needed.
- **M2:** the Qt-free `chemsuite/metabolism/` package — `features.py` (33 atom features, 13 bond
  features, every upstream quirk documented as `ionisation/molgpka.py` documents its own; the
  quirks found while reading the source on 2026-11-17 are listed under **Upstream quirks** below),
  `graphnet.py` (the numpy forward pass), `sites.py` (`ISOFORMS`, `MetaboliteSite`,
  `predict_sites`, `predict_on` remapping atom indices through `_smilesAtomOutputOrder`),
  and `predict.py` (the chunked, cancellable batch driver, copied from `ionisation/predict.py`).
  `data/som_weights.npz` is loaded by `sites.py` through `importlib.resources` and
  `functools.cache`, as `ionisation/molgpka.py` loads its own archives — **there is no
  `weights.py`**, since nothing is fetched. Added to `QT_FREE` in
  `tests/test_import_boundaries.py`.
- **M3:** the weights. A trained fold is **34.15 MB** (8.5M float32 parameters, measured
  2026-11-17), so:
  - the 10-fold ensemble is about 341 MB float32 or 171 MB float16 — a download, as the pKa
    weights are (88 MB), with a GitHub release `som-weights-1`, a SHA-256 pin, and the
    once-per-session consent of `tools/macromolecule/pka_model.py`;
  - **one fold at float16 is about 17 MB, which could be vendored in-tree** beside
    `ionisation/data/molgpka_*.npz` (5.5 MB each). That removes the download, the consent
    prompt, the cache folder and the network-policy path altogether — the model would just
    work, offline, first time.

  **Measured 2026-11-17** (`scripts/som/holdout_eval.py`). A k-fold run cannot answer this —
  every other fold trained on any given fold's test molecules — so the external set is
  AZ-ExactSOM: 120 hepatocyte compounds, **zero overlap** with Zaretzki, 95 of them with exact
  sites. The ten grouped folds, averaged over random subsets of each size:

  | folds | float32 | float16 | top-1 | top-2 | MCC |
  |---|---|---|---|---|---|
  | 1 | 34 MB | **17 MB** | 0.427 ± 0.032 | 0.503 ± 0.035 | 0.225 |
  | 2 | 68 MB | 34 MB | 0.462 ± 0.035 | 0.517 ± 0.028 | 0.233 |
  | **3** | 102 MB | **51 MB** | **0.480 ± 0.027** | **0.525 ± 0.031** | 0.250 |
  | 5 | 171 MB | 85 MB | 0.474 ± 0.021 | 0.520 ± 0.018 | 0.236 |
  | 7 | 239 MB | 120 MB | 0.475 ± 0.015 | 0.539 ± 0.019 | 0.242 |
  | 10 | 342 MB | 171 MB | 0.495 | 0.558 | 0.234 |

  So **the ensemble is worth something, and three folds capture most of it**: 1 → 3 buys +0.053
  top-1 (about 1.6 single-fold standard deviations), while 3 → 10 buys a further +0.015 top-1 and
  +0.033 top-2 for 240 MB more. MCC is flat throughout, so the benefit is in *ranking* the sites,
  which is exactly what the tool displays.

  These absolute numbers are cross-dataset transfer with no isoform labels, so they are a
  pessimistic floor, not the figure to quote — only the shape of the curve is being used here.

  **float16 verified harmless, 2026-11-18** (`scripts/som/convert_weights.py --check`, 25
  molecules × 9 isoforms). Storing every tensor of ≥65,536 elements as float16 and the rest as
  float32 gives **17.37 MB** (98.3% of parameters halved; float32 throughout is 34.12 MB):

  | | |
  |---|---|
  | largest logit difference | 0.001471 |
  | largest probability difference | 0.000357 |
  | top-1 site changed | 1 of 225 — and its float32 top-1/top-2 gap was **exactly 0.0** |

  The single reordering is between two atoms the float32 model scored *identically*, so it is a
  tie-break, not a loss. **float16 ships.**

  That tie is itself worth noting for M4: **symmetry-equivalent atoms get identical scores**
  (equivalent methyls, ring positions), so ties are common and the displayed "Top Site" must not
  wobble between runs. The port should break ties by lowest atom index, and the ranked-sites
  column should ideally group equivalent atoms rather than pick one arbitrarily.

  **Decided with the user, 2026-11-18: one fold, vendored in-tree.** About 17 MB at float16,
  committed to `chemsuite/metabolism/data/` beside the `LICENSE-ATTNSOM.txt` notice, loaded with
  `importlib.resources` and `functools.cache` exactly as `ionisation/molgpka.py` loads its own
  5.5 MB archives. The tool then works offline, first time, for everyone.

  **So most of M3 disappears**: no GitHub release, no `SHA256`/`SIZE` pin, no `cached`/`fetch`
  pair, no `ctx.cache_folder("som")`, no `ctx.network_allowed()` check, no once-per-session
  consent prompt, no `weights.py` at all, and no download path to test. `metabolism/weights.py`
  is dropped from the M2 file list; the loader lives in `sites.py`. That is a large simplification
  bought for about five points of top-1, and a later `som-weights-2` release can still add the
  three- or ten-fold ensemble as an optional download if anyone wants it.
- **M4:** the tool, `chemsuite/tools/metabolism/`, on the two-file Predict pKa pattern: `Tool(id=
  "predict_som", title="Sites of Metabolism…", menu=("Chemistry", "Calculate"))`, a kit dialog with
  isoform choice (3A4, 2D6, 2C9, 1A2 on by default; the other five in a collapsed section) and an
  `OUTPUTS` dict writing, per isoform, the ranked sites, the top site, its score and an
  applicability-domain value.
- **M5:** showing the sites on the structure. `core/depict.py: render_png` already takes
  `highlight_atoms`, but `app/table/table_view.py`'s `Highlight` is a column-wide SMARTS query and
  `StructureImages`' cache key is the molecule blob, so per-row highlights would collide. Add an
  index-provider `Highlight` variant with a per-row cache-key component, and `highlightAtomColors`
  for graded scores. Shared rendering, so the full suite runs at the end of this step.
- **M6:** `docs/USER_GUIDE.md` in the same commit as the tool, and the results above.

## Phase C — CYP inhibition

After M, reusing `metabolism/features.py`, `graphnet.py` and `predict.py` with a
molecule-level head instead of a per-atom one. No open-weights model is worth shipping, but the
open data now is: OpenADMET's CYP challenge set (`openadmet/cyp-challenge-train-test`, Apache-2.0;
measured pIC50 for 1A2, 2C9, 2D6 and 3A4, about 6.4k compounds, plus `tdi` and 62.9k
single-concentration rows) **pinned to a dated snapshot**, since the organizers amended it
mid-challenge; and the Veith panel (Veith et al., Nat. Biotechnol. 2009) pulled straight from
PubChem BioAssay rather than through TDC, so redistribution terms never arise — this is what adds
2C19, and gives about 12.1–12.6k compounds per isoform.

The winning open entry in that challenge (thinkMolecularTech) is a multitask Chemprop D-MPNN plus
LightGBM ensemble on ECFP4, 200 RDKit 2D descriptors **and predicted pKa** — all three of which
chemsuite already computes (`analysis/fingerprints.py` with its per-table `FingerprintCache`,
`descriptors/catalog.py`, and our own pKa model). So: measure a scikit-learn baseline on those
features first (sklearn is already a base dependency, so it needs no port at all), and only train
the D-MPNN if it beats the baseline by enough to be worth a second weights release. ATTNSOM's code
already offers a `chemprop` encoder, so the numpy runtime from M2 serves both.

## The data, checked (2026-11-17)

Read straight out of the clone's SDFs with chemsuite's own RDKit (2026.09.2, not upstream's
2025.3.6), which parsed every record:

| | 1A2 | 2A6 | 2B6 | 2C8 | 2C9 | 2C19 | 2D6 | 2E1 | 3A4 |
|---|---|---|---|---|---|---|---|---|---|
| molecules | 271 | 105 | 151 | 142 | 226 | 218 | 270 | 145 | 475 |
| primary sites | 410 | 157 | 210 | 200 | 325 | 299 | 368 | 218 | 666 |

2,003 molecule–isoform instances over 679 unique canonical SMILES, which is exactly what the paper
claims. Every instance has a usable `PRIMARY_SOM`, no site index is out of range, and **no record
carries an explicit hydrogen** — confirming the heavy-atom-only finding below.

**The SDFs also carry `SECONDARY_SOM` and `TERTIARY_SOM`, and upstream reads neither.** Only
`PRIMARY_SOM` becomes a positive label, so known-but-secondary metabolism sites are trained and
scored as negatives. That matters twice over: it caps the achievable precision (some "false
positives" are real sites), and it is the most obvious cheap improvement available to us — soft
labels, or a secondary-site auxiliary head. M1 measures the plain reproduction first, then tries
it, and the plan records whether it helped.

## How we compare — and an open challenge to the architecture choice (2026-11-18)

**There is one directly relevant independent benchmark**: `WillDeSnoo/som-benchmark` (MIT), which
evaluates eight approaches on **the same Zaretzki set** with **Bemis–Murcko scaffold splits** (10
seeds) *and* a second split policy, all heavy atoms as candidates, one multi-isoform model per
approach, ties broken by seeded permutation. It does not redistribute the data (it pulls it from
ATTNSOM's repo at a pinned commit) and reaches the same conclusion we did about the usual splits.

Its top-1 figures, beside our grouped ATTNSOM:

| approach | random/grouped top-1 | scaffold top-1 | scaffold top-2 |
|---|---|---|---|
| atom-type frequency (baseline) | 0.419 | 0.425 | 0.607 |
| **SMARTCyp (no training at all)** | 0.629 | **0.635** | **0.770** |
| fingerprint + GBM | 0.622 | 0.601 | 0.741 |
| atom/bond features + GBM | – | 0.660 | 0.795 |
| **Hückel physics + GBM** | 0.684 | **0.660** | **0.811** |
| **Hückel + BDE + GBM** | 0.706 | **0.674** | **0.812** |
| message-passing GNN | 0.682 | 0.643 | 0.788 |
| **our ATTNSOM, molecule-grouped** | **0.618** | – | **0.743** |

**This is uncomfortable, and it has to be resolved before M2.** On its face, a plain gradient
boosting machine over Hückel features beats our 8.5M-parameter graph network on the harder split
(0.660/0.811 vs 0.618/0.743), and so does **SMARTCyp, which has no training data whatsoever**.
Their trained models lose only 2–4 points going from their easier to their harder split, where
ATTNSOM lost 14 going from instance to molecule-grouped — which is itself a hint that ATTNSOM
leans more on memorizing molecules than on learning chemistry.

**Caveats, which are why this is not yet a verdict:**
- Their top-k denominator is not fully specified on the page ("per-molecule top-k", pooled, but
  the record-vs-molecule aggregation is not written down), while ours is per molecule-isoform
  record from upstream's `metrics.py`. If theirs is per *molecule*, the numbers are not comparable.
- Their tie handling (seeded permutation, plus a pessimistic variant) differs from upstream's
  `argsort`, and we have just shown ties are common on symmetric molecules.
- Their "random" split may already be molecule-grouped; the tiny gap for untrained SMARTCyp
  (0.629 → 0.635) suggests it is, but it is not stated.

**So M1 gains a step 6: run `som-benchmark` ourselves and evaluate our model inside its harness.**
It is MIT, it pulls the data from the clone we already have, and it would give one table with
ATTNSOM, SMARTCyp and the GBMs measured identically. Only that settles whether to keep ATTNSOM.

**If the GBM result holds up, the plan should change, and change for the better:** Hückel theory
is an eigendecomposition of a connectivity matrix — a page of numpy — and gradient boosting is
scikit-learn, which is already a base dependency (`docs/PLAN.md`). That means **no neural port, no
`graphnet.py`, no `SAGPooling`, no 17 MB of vendored weights, and a model measured in kilobytes**,
for better accuracy. The BDE variant is the one to avoid: its bond dissociation energies come from
ALFABET, which is a neural net of its own, and Hückel alone gives up only 0.014 top-1.

## Commercial tools, for context

No current, directly comparable commercial numbers exist, and the honest answer is that we cannot
rank ourselves against them:

- **MetaSite vs StarDrop** (Shin et al., *Drug Metab. Lett.* 2011): on CYP3A4/2C9/2D6 substrates,
  StarDrop 89% (2C9) and 93% (2D6) against MetaSite 63–70%, scored "first choice" and rank-weighted
  rather than top-k. Two versions of the same abstract report different MetaSite figures.
- **Meteor / MetaSite / StarDrop** (2011): roughly 50% sensitivity and precision for MetaSite and
  StarDrop on their own set.
- **SMARTCyp's own paper** (2012) put CYP2C9 top-1/2/3 at 66/86/87% against StarDrop 5.0's
  61/73/77% — note this flatters SMARTCyp more than the independent benchmark above does.

All of it is 14–15 years old, on private or undisclosed compound sets, with scoring that is not
top-k. It is useful only as a rough band (first-choice accuracy in the 60–90% range), and it is
**not** something the user guide should cite as a comparison.

## Chemistry sanity check (2026-11-18)

Before trusting any metric, the golden file's own predictions were read against known
metabolism. Grouped fold 1, top three atoms per isoform:

| molecule | top-ranked sites | the literature |
|---|---|---|
| caffeine | atoms 13, 9, 0 — all three N-methyl carbons (1A2 0.24) | N-demethylation to paraxanthine is the dominant 1A2 route; atom 2, the C8, appears for 3A4, which is the trimethyluric acid path |
| paracetamol | atom 3, the amide N, 0.87 | N-oxidation to NAPQI — the pathway that makes it hepatotoxic |
| ibuprofen | the isobutyl carbons, and **2C9 scores them highest** (0.75) | 2- and 3-hydroxylation of the isobutyl chain, and CYP2C9 is indeed ibuprofen's main metabolizer |

Three molecules, three correct primary sites, and the correct isoform preference for ibuprofen
without being told it. Caffeine's absolute numbers are low (0.24) because it is small with several
near-equivalent sites — the *ranking* is what the tool shows, and it is right.

Ibuprofen also shows the tie problem concretely: atoms 0 and 2 (the two equivalent terminal
methyls) score **identically** to every decimal, on every isoform.

## `SAGPooling` in full, for the port (read from PyG 2.6.1, 2026-11-18)

The one operation whose definition is not in upstream's own file. `SAGPooling(256, ratio=0.8)`
with the default `GraphConv` and `tanh`, feeding the FiLM graph vector. Per graph in the batch:

1. **Score each atom.** `attn = GraphConv(256 → 1)(h, edge_index)`, i.e.
   `attn_i = lin_rel(Σ_{j∈N(i)} h_j) + lin_root(h_i)`, where `lin_rel` carries a bias and
   `lin_root` does not (weights `sagpool.gnn.lin_rel.{weight,bias}`, `sagpool.gnn.lin_root.weight`;
   aggregation is a plain sum, and our edge list already holds both directions).
2. **Project and squash.** `score = tanh((attn · p) / ‖p‖₂)` with `p = sagpool.select.weight`,
   shape (1, 1) — so with one channel this reduces to `tanh(attn · sign(p))`.
3. **Keep the best 80%.** `k = ceil(0.8 · n_atoms)`, the top `k` atoms by score.
4. **Scale what is kept.** `h_kept = h[perm] * score[perm]` — the kept features are multiplied by
   their own tanh score, so the next step is a **score-weighted** mean, not a plain one.
5. **Pool.** `g = mean(h_kept)` over the kept atoms of that graph.

**Tie-breaking is a real hazard here.** `topk` sorts descending with `torch.sort`, which is not
stable by default, so when two atoms score identically the atom kept at the `k`-th boundary is not
guaranteed. Symmetry-equivalent atoms score identically often (see the float16 note above). The
port must break ties by **lowest atom index** and the golden file must include a symmetric molecule
so the choice is pinned.

## Upstream quirks (read from the source, commit `0faa38b`, 2026-11-17)

The port must reproduce these exactly, because the trained weights expect them. Where upstream has
a plain bug on a path that never runs, the port documents it and does not copy it.

- **Heavy atoms only, no hydrogens.** `dataset.py` reads the SDFs with `Chem.SDMolSupplier(path)`
  and inference uses `Chem.MolFromSmiles` — both default to `removeHs=True`, and neither ever calls
  `AddHs`. So every atom index is a **heavy-atom** index and `ATOM_LIST`'s hydrogen entry is dead
  weight the weights never see. (This is unlike `ionisation/molgpka.py`, which features `mol_h`.)
  It also means `predict_on` has an easier job than the pKa port's: no hydrogen bookkeeping.
- **No fluorine in the heteroatom count.** `hetero_neighbor_count` uses `{7, 8, 16, 17, 35, 53}`
  — N, O, S, Cl, Br, I. Fluorine (9) is **absent**, though it is in `ATOM_LIST`. Almost certainly
  an oversight upstream; reproduce it, because the weights were trained with it.
- **`bond_features(None)` returns 12 values, not 13** (it omits the "other" bond-type slot that
  `bond_one_hot(bt_idx, len(_BOND_TYPES) + 1)` adds on the live path). `mol_to_graph` never calls
  it with `None`, so the trained weights are unaffected. Do not port the broken branch.
- **Bond stereo falls back to index 0** (`STEREONONE`) for any unlisted value, rather than to an
  all-zero vector.
- **Edge order matters.** Bonds are walked in RDKit bond order and each becomes two directed edges
  appended as `(i→j, j→i)`, so the reverse of edge `2k` is `2k+1`. `rev_edge_index` encodes that,
  and is used only by the Chemprop D-MPNN encoder — the GraphCliff default ignores it.
- **Private RDKit APIs.** The logP and TPSA per-atom contributions come from
  `Crippen._GetAtomContribs` and `rdMolDescriptors._CalcTPSAContribs` (both underscore-private, and
  TPSA without the S and P contributions). Upstream pins `rdkit==2025.3.6` and chemsuite ships
  `2026.09.2`, so a drift here would silently shift every prediction.
  **Checked 2026-11-17: no drift.** Both versions give byte-identical contributions on 12 diverse
  molecules (caffeine, paracetamol, ibuprofen, a benzodiazepine, parathion for P and S, carbazole,
  a sulfoxide ylide, a CF3 aniline, thiophene, salicylic acid and a bare salt). Both APIs still
  exist and behave the same. The golden test keeps covering them as a regression guard, since they
  are private and could change without notice.
- **Isoform order is fixed** by `main.py`: `['1A2', '2A6', '2B6', '2C8', '2C9', '2C19', '2D6',
  '2E1', '3A4']`. The isoform embedding rows are indexed in exactly that order, so
  `metabolism/sites.py`'s `ISOFORMS` must match it.
- **Shapes, from `main.py`'s defaults:** hidden 256, `num_layers` **4**, `groups` 4, `mid_K` 3,
  4 attention heads, `attn_tau` 0.1, train batch 32. Inside each of the 4 `GraphCliffFilter`
  layers the projection widens to `3 × 256 = 768` and **`ShortGINE` runs at 768**, so its GIN MLP
  is 768→768→768 and its edge MLP is 13→768→768; the split into `x2, x1, v` happens *after*
  `ShortGINE`, and only `x2` (256 wide) goes through `LongPoly`.
- **`GraphFiLM`'s hidden layer is 128**, not 256: the MLP is 256→128→512 (`hidden: int = 128`
  default, and `ATTNSOM` passes only `node_dim` and `graph_dim`).
- **`LongPoly`'s own default is K=5**, but `GraphCliffFilter` passes `mid_K=3`, so there are four
  Chebyshev terms (T0…T3). The recursion starts `T0 = x`, `T1 = P·x`, `T_k = 2·P·T_{k-1} − T_{k-2}`,
  each term scaled per channel-group, then `· group_scale + group_bias`, then
  `dropout(SiLU(RMSNorm(·)))`.
- **The propagation matrix uses unit edge weights** normalized by out-degree:
  `normalize_edges` scatter-adds ones over `src`, so `Ã_e = d_src^-1/2 · d_dst^-1/2`. `ShortGINE`
  separately uses `deg = bincount(dst).clamp(min=1.0)` and weights each message by
  `1/(deg_src · deg_dst)^0.5` — two different degree conventions in one layer, both intentional.
- **The head returns logits**, with no sigmoid. Apply the sigmoid in the port, once.
- **`GraphCliffMultiRegressor` is an alias** for `ATTNSOM`; older checkpoints use that name, so
  the weight converter should accept both prefixes.

## Notes

- **No torch in the app, ever.** Training happens in WSL; the app runs the exported `.npz` in
  numpy, as `ionisation/unipka` does.
- **BLAS threads:** pin workers to one BLAS thread. numpy's multithreaded BLAS loses on matrices
  this small, and a broadcast zero-stride array makes numpy skip BLAS entirely — both recorded in
  the pKa port.
- **`som` is taken.** `chemsuite/analysis/som.py` is self-organizing maps with a `Som` class.
  This package is `chemsuite/metabolism/`, and the result type is `MetaboliteSite`.
- **It is a screening estimate**, not a measured value. Say so in the tooltips and the user guide,
  in the wording `descriptors/catalog.py`'s `PKA_NOTE` already uses — and ship the
  applicability-domain column, as the QSAR roadmap entry expects of any prediction.
