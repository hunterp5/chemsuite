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
     so the quoted 10-fold CV mean remains an honest estimate of it.

     **The production run is done (2026-11-18)**: 20 grouped folds at a 32-epoch cap, so each
     model trains on 95% of the molecules (`runs/production`, `fold1_best.pt` the candidate).

     | run | folds | train share | MCC | top-1 | top-2 | top-3 |
     |---|---|---|---|---|---|---|
     | grouped CV | 10 | 90% | 0.524 ± 0.061 | 0.618 | 0.743 | 0.811 |
     | **production** | **20** | **95%** | **0.537 ± 0.067** | **0.606** | **0.746** | **0.815** |

     The extra 5% of training data moves nothing outside seed noise, which is worth knowing in
     itself: **this architecture is not data-starved on Zaretzki**, so "train on more of it" is not
     the lever. Note the per-isoform spread stays wide (2C19 top-1 0.679 against 3A4 0.561) and 3A4
     — the isoform that matters most clinically — remains the worst, as it was in the grouped CV.
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
- **M5 is done (2026-11-18): the sites are marked on the drawings**, and the feature is the
  app's, not the tool's — a structure column can mark the atoms *any* column names, so nothing in
  `app` knows what wrote them.
  - `core/atom_names.py` (Qt-free) reads the convention the tool already wrote: comma-separated
    `C7`, optionally with a score (`C7 0.81, N3 0.05, C12 0.04`). `atoms_named` checks each name
    against the molecule beside it — the element must be that atom's element and the index inside
    the molecule — so a column written for other structures marks nothing instead of the wrong
    atoms.
  - `core/depict.py` grew `graded_colors` (the best site at full strength, each weaker *score*
    fainter, blended toward white or the dark base; equal scores share a shade, so
    symmetry-equivalent atoms do not invent a ranking), `drawn_highlight` and
    `render_png(highlight_colors=…)` → RDKit's `highlightAtomColors`. `depiction_mol` now tags
    each kept atom with its index in the molecule as given, since removing a 3D structure's
    hydrogens renumbers the rest; the tag is how a site's index survives into the drawing.
  - `app/table/highlights.py` is a new module that owns every kind of cell highlight: the old
    `Highlight` (a column-wide substructure query, Structure Search's), the new `AtomSource` (the
    setting: which column names the atoms) and `AtomHighlight` (one row's resolved atoms, whose
    `key` carries the cell text, so `StructureImages`' cache key varies per row). `table_view.py`
    came out 3 lines *shorter* than before the step; `ColumnActions` now takes the view itself
    instead of four callbacks.
  - The entry points: **Highlight Atoms ▸** in a structure column's header menu (listing only
    columns that name atoms of those structures), Clear Highlighting, and the tool itself, which
    marks the most prescribed isoform it was asked for (3A4, 2D6, 2C9, 1A2) when a run finishes —
    the ranked `Sites` column if it was written, else `Top Site`. `ToolContext.highlight_atoms`
    is the only way a tool reaches it.
  - Saved with the window layout (`atom_highlights`), as Align Drawings is, so a reopened
    document still marks them.
  - Tests: `tests/core/test_atom_names.py`, two in `tests/core/test_depict.py`,
    `tests/app/table/test_atom_highlight_ui.py` and one in `tests/tools/metabolism`. The full
    suite was run (shared rendering).
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

**The three caveats were checked against their source (2026-11-18), and all three resolve against
ATTNSOM:**

- **The metric is identical.** `src/som/metrics.py`: "top-1 / top-2 **per-instance** accuracy — the
  fraction of (molecule, isoform) instances for which at least one true site appears among the k
  highest-scoring atoms", instances without a positive skipped (there are none). That is upstream
  ATTNSOM's `topk_accuracy` definition, over the same candidate set of all heavy atoms.
- **The split is identical.** `src/som/splits.py`: "**Every split here partitions molecules.**" In
  `random` mode the group key is the SMILES itself, so **their random column is exactly our
  molecule-grouped split**. They reached our leakage conclusion independently, in the same words:
  "Splitting on instances would put the same structure on both sides."
- **The tie policy favors us, not them.** They break ties by seeded random permutation and note
  that breaking them in index order "silently flatters a model whose first atoms happen to be the
  labelled ones" — which is exactly what upstream's `argsort` does, and what produced our 0.618.

So on the same data, the same metric and the same split policy, with the tie rule tilted our way:

| approach, molecule-grouped | top-1 | scaffold top-1 | scaffold top-2 |
|---|---|---|---|
| Hückel + BDE + GBM | 0.706 | 0.674 | 0.812 |
| **Hückel physics + GBM** | **0.684** | 0.660 | 0.811 |
| message-passing GNN | 0.682 | 0.643 | 0.788 |
| SMARTCyp (untrained) | 0.629 | 0.635 | 0.770 |
| fingerprint + GBM | 0.622 | 0.601 | 0.741 |
| **our ATTNSOM** | **0.618** | – | 0.743 |
| atom-type frequency | 0.419 | 0.425 | 0.607 |

### That verdict was wrong: the labels are not the same labels (corrected 2026-11-18)

**A fourth difference, found after the above was written and committed, invalidates the whole
table.** `scripts/03_evaluate.py` takes `--label-field`, and its default — the one every committed
number in their `results/` was produced with — is **`som_all`**, where
`SOM_FIELDS = ("PRIMARY_SOM", "SECONDARY_SOM", "TERTIARY_SOM")` and `som_all` is their union.
Upstream ATTNSOM, and therefore every number of ours, trains and scores on **`PRIMARY_SOM` alone**.

Counted over the nine SDFs (`/home/hp/bench/probe/labels.py`):

| labels | positive atoms | rate |
|---|---|---|
| `PRIMARY_SOM` only — ours and upstream's | 2,853 | **0.0681** |
| all three tiers — theirs | 4,163 | **0.0994** |

Their rate of 0.0994 is exactly the 9.94% their README quotes and the 0.0989 mean `positive_rate`
in their `metrics_per_run.csv`, so this is confirmed, not inferred. **Their task has 46% more
atoms counting as a correct hit.** Top-1 asks "is the highest-scoring atom *an* annotated site",
so a 46% larger target set inflates it mechanically. Their 0.684 and our 0.618 are answers to
different questions, and the gap cannot be read as a difference in model quality.

The three differences checked earlier still hold, and still favor us on tie-breaking — but they
are now beside the point next to this one.

**So the pivot is not justified yet, and the decisive experiment is a single flag:** re-run their
harness with `--label-field som_primary`, which scores *their* models under *our* label
definition, directly against our 0.618. Their code supports it; nothing has to be written. Until
that number exists, ATTNSOM stays, and the paragraphs below are a hypothesis, not a decision.

### The experiment answers it: the pivot is dead, ATTNSOM wins (2026-11-18, seed 0)

Ran with `scripts/som/{setup_bench_env.sh,run_bench.sh}`. First seed of the `som_primary` pass,
beside their published `som_all` figures for the same models:

| model | top-1 `som_primary` | top-2 `som_primary` | their `som_all` top-1 | label effect |
|---|---|---|---|---|
| atom-type prior | 0.337 | 0.500 | 0.419 | −0.082 |
| fingerprint + GBM | 0.507 | 0.613 | 0.622 | −0.115 |
| atom/bond + GBM | 0.566 | 0.682 | – | – |
| **Hückel + GBM** | **0.566** | **0.675** | **0.684** | **−0.118** |
| **our ATTNSOM, grouped** | **0.618** | **0.743** | – | – |
| **our ATTNSOM, production** | **0.606** | **0.746** | – | – |

**Switching to the primary-only labels costs their models 0.08–0.12 top-1, uniformly**, which is
the mechanical effect predicted above and confirms the diagnosis. On the label definition our tool
actually uses, **ATTNSOM leads the best GBM by about 0.04–0.05 top-1 and by 0.07 top-2** — the
reverse of what their README appeared to say.

So: **keep ATTNSOM. M2 proceeds as planned**, with `graphnet.py`, the `SAGPooling`
reimplementation, the vendored float16 weights and the activation golden file all back in scope.
The Hückel + GBM is worth remembering as a cheap, well-calibrated baseline — and `som-benchmark`
remains the best external check we have — but it is not the better model here.

### …but that table is also biased, in our favor (checked 2026-11-18)

Running the `som_all` control exposed a second problem with the seed-0 table above. Two things
came out of comparing it against their committed numbers:

| model | their committed `som_all`, seed 0 | our `som_all`, seed 0 | drift | our `som_primary`, seed 0 |
|---|---|---|---|---|
| atom-type prior | 0.4198 | 0.4198 | +0.0000 | 0.3373 |
| fingerprint + GBM | 0.5849 | 0.5849 | +0.0000 | 0.5071 |
| atom/bond + GBM | 0.6344 | 0.6156 | −0.0189 | 0.5660 |
| Hückel + GBM | 0.6439 | 0.6462 | +0.0024 | 0.5660 |
| message-passing GNN | 0.6344 | 0.6415 | +0.0071 | 0.5307 |
| GNN + Hückel | 0.6368 | 0.6509 | +0.0142 | 0.5094 |

**Good news — the environment is sound.** `prior` and `fp_gbm` reproduce their committed values
exactly, and the trained models drift by at most 0.019, mostly under 0.01, which is the expected
cost of rdkit 2026.09.2 and numpy 2.5.5 against their macOS pins. Numbers from this environment
can be trusted.

**Bad news — seed 0 is a hard seed.** Their 10-seed `som_all` mean for the Hückel GBM is 0.684,
but seed 0 gives 0.644: about 0.04 low. So comparing a seed-0 `som_primary` figure (0.566)
against our *10-fold mean* (0.618) is not a fair comparison, and the bias runs **in our favor**.
Correcting for the seed penalty puts their Hückel GBM near 0.61 on primary-only labels, which is
a **tie** with our 0.606–0.618, not a loss.

**So the honest status is: unresolved, and much closer than either earlier table suggested.** The
10-seed `som_primary` run is what settles it, and nothing should be concluded until it lands.
`fp_gbm` is dropped from it — one seed took 96 minutes against 20 seconds for `phys_gbm`.

### The settled answer: a tie on accuracy (10 seeds, 2026-11-18)

`scripts/som/bench_table.py` builds this, and refuses to build it at all unless each result file
holds only its own run's rows. Control first — our `som_all` against their committed `som_all`,
both 10-seed means:

| model | theirs | ours | drift |
|---|---|---|---|
| atom/bond + GBM | 0.7025 | 0.7009 | −0.0017 |
| message-passing GNN | 0.6824 | 0.6835 | +0.0011 |
| GNN + Hückel | 0.6726 | 0.6710 | −0.0016 |
| Hückel + GBM | 0.6838 | 0.7025 | +0.0187 |
| atom-type frequency | 0.4192 | 0.4192 | +0.0000 |

Largest drift 0.019, three models inside 0.002 — **the environment reproduces theirs**, so its
`som_primary` figures can be trusted. And the label field alone costs **−0.100, −0.103, −0.109,
−0.116 and −0.101 top-1** across the five models: uniform, large, and exactly the mechanical
effect diagnosed above. That is the whole of the original 0.684-vs-0.618 "gap".

The comparison, on the labels our tool uses:

| model | top-1 | top-2 |
|---|---|---|
| atom/bond + GBM | 0.601 ± 0.047 | 0.730 ± 0.048 |
| **Hückel + GBM** | **0.599 ± 0.046** | **0.744 ± 0.048** |
| message-passing GNN | 0.575 ± 0.032 | 0.703 ± 0.042 |
| GNN + Hückel | 0.555 ± 0.042 | 0.688 ± 0.049 |
| atom-type frequency | 0.318 ± 0.040 | 0.493 ± 0.029 |
| **our ATTNSOM, grouped 10-fold** | **0.618** | **0.743** |
| **our ATTNSOM, production 20-fold** | **0.606** | **0.746** |

**ATTNSOM is ahead by 0.017 top-1 against a seed spread of ±0.047, and level on top-2 (0.746 vs
0.744). That is a tie.** Neither architecture has an accuracy case over the other.

**Two residual biases both run in our favor**, so if anything the tie flatters ATTNSOM:

1. **Training share.** Their split is 475/68/136 molecules, so ~70% train; our 10-fold trains on
   ~85% and the production run on 95%. We give our model 15–25% more data for the same test task.
   The production run suggests this is worth little — 95% scored *below* 90% — but it is not zero.
2. **Tie-breaking.** Our numbers come from upstream's `argsort`, which resolves ties by atom index;
   theirs use a seeded permutation. Their own docstring notes index order "silently flatters a
   model whose first atoms happen to be the labelled ones", and we have measured that ties are
   common on symmetric molecules.

**What that means for the decision.** Accuracy no longer separates the two, so the choice falls to
everything else — and there the gap is not close:

| | ATTNSOM | Hückel/atom-bond + GBM |
|---|---|---|
| runtime dependency | a hand-written numpy port of 8.5M parameters | scikit-learn, already a base dependency |
| new code | `graphnet.py`, `features.py`, a `SAGPooling` reimplementation | `features.py` plus a Hückel eigendecomposition |
| shipped weights | ~17 MB vendored float16 | kilobytes |
| test fixture | an activation golden file, layer by layer | predicted probabilities |
| calibration | focal-loss logits, uncalibrated | ECE/Brier reported, and GBMs calibrate readily |
| upstream risk | must reproduce PyG's `topk` tie behavior exactly | none |

**The one measurement that would still sharpen this** is running ATTNSOM inside their harness, on
their splits with their tie policy, which removes both residual biases. It is worth doing only if
the decision hinges on a difference smaller than 0.05 top-1 — and since the engineering case is
now what decides, it probably does not.

One more trap, worth recording because it nearly corrupted the comparison: `03_evaluate.py`
*updates* rows in an existing `metrics_per_run.csv` keyed by (model, mode, seed) instead of
rewriting the file. A short run therefore leaves their committed `som_all` rows in place for every
model, mode and seed it does not recompute, so the output file silently mixes their numbers with
ours — which is exactly what the first smoke run produced. `run_bench.sh` now clears
`metrics_per_run.csv` and `metrics_summary.csv` before each pass.

## Decided with the user, 2026-11-18: pivot to the GBM, primary-only labels

Accuracy is a tie, so the decision was taken on everything else, and **M2–M6 below are
superseded by this section**. Two further findings made the choice easier than the table above
suggests.

**The Hückel descriptors are not worth having.** Their own `docs/HUCKEL_DESCRIPTORS.md` says so,
in a box headed *"Do these features earn their place? Measured: no"*: dropping all 21 of them and
refitting costs +0.001 top-1 (p = 0.96), and the document explicitly "documents what the
descriptors are, not a claim that they help". Our own 10-seed run agrees — `basic_gbm` 0.601
against `phys_gbm` 0.599. **So the target is the plain atom-feature GBM**, and the
eigendecomposition, the Streitwieser heteroatom parameters, the degenerate-frontier averaging and
the π-fragment handling all disappear with it. What looked like the clever part of the pivot turns
out to be the part to leave out.

**And the feature set is 34 atom descriptors, all cheap RDKit calls** (read from
`basic_atom_block` and `charge_block`, 2026-11-18), plus a 9-wide isoform one-hot — 43 columns:

| block | columns |
|---|---|
| element one-hot | `el_C N O S P F Cl Br I other` (10) |
| hybridization one-hot | `hyb_sp sp2 sp3 other` (4) |
| topology and valence | `n_hydrogens`, `heavy_degree`, `formal_charge`, `is_aromatic`, `in_ring`, `min_ring_size`, `n_rings_at_atom`, `aromatic_neighbours`, `hetero_neighbours`, `has_double_bond`, `has_triple_bond`, `is_benzylic_or_allylic_CH`, `is_terminal` (13) |
| charge and lipophilicity | `gasteiger_charge`, `estate_index`, `crippen_logp`, `crippen_mr`, `neighbour_charge_max/min/mean` (7) |
| isoform | one-hot over the nine (9) |

There is **no bond block**: bond information enters through `has_double_bond`, `has_triple_bond`
and the neighbour counts. Features are computed once per unique molecule and reused across the
nine isoforms, which is also how the app should batch them.

**The estimator**, read from `_fit_gbm`: `HistGradientBoostingClassifier(max_iter=300,
early_stopping=False, l2_regularization=1.0, min_samples_leaf=20)`, with `learning_rate` ∈ {0.05,
0.10} × `max_leaf_nodes` ∈ {31, 63} chosen on the validation split by average precision.

**Labels: `PRIMARY_SOM` only** (decided with the user). It is what our measurements use, it is the
harder and more defensible question, and `SECONDARY_SOM`/`TERTIARY_SOM` would inflate every quoted
figure by about 0.10 top-1. The user guide must state the label definition beside each number, so
no reader mistakes ours for a published `som_all` figure.

### What M2 becomes

A Qt-free `chemsuite/metabolism/`, added to `QT_FREE` in `tests/test_import_boundaries.py`:

| file | contents |
|---|---|
| `features.py` | The 43 columns above, ported from `som-benchmark`'s `features.py` (MIT) with attribution. Keep its `_safe()` guard: Gasteiger charges can come back NaN or infinite. |
| `trees.py` | Numpy inference over the exported ensemble: walk `feature`, `threshold`, `left`, `right`, `is_leaf`, `value` per tree, sum the leaf values onto `baseline_prediction`, then a logistic. ~40 lines, and the only "port" left. |
| `sites.py` | `ISOFORMS`, `MetaboliteSite(atom, isoform, score, domain)`, `predict_sites(mol, isoforms)`, `predict_on(mol)`. Loads `data/som_gbm.npz` with `importlib.resources` + `functools.cache`. **Ties broken by lowest atom index** — symmetry-equivalent atoms score identically. |
| `predict.py` | The chunked, cancellable batch driver copied from `ionisation/predict.py`, unchanged by the pivot. |
| `data/som_gbm.npz` | The trees, order of 100 KB rather than 17 MB. |
| `data/LICENSE-SOMBENCH.txt` | `som-benchmark`'s MIT text, beside the other vendored licences. |

**No `weights.py`, no download, no consent prompt, no cache folder** — and no `graphnet.py`, no
`SAGPooling`, no float16 archive, no activation golden file. The golden fixture becomes predicted
probabilities for a handful of molecules, a few kilobytes of JSON.

**Deliberately not pickle.** A pickled `HistGradientBoostingClassifier` would be opaque, tied to
one scikit-learn version, and a load-time hazard; exporting the tree arrays to `.npz` and
evaluating them in numpy keeps the house pattern (`ionisation/molgpka.py`) and survives sklearn
upgrades.

**Training moves out of WSL entirely.** scikit-learn is already a base dependency, so
`scripts/som/train_gbm.py` runs in the app's own venv — no conda, no CUDA, no torch, and a fit
takes seconds. `scripts/som/export_gbm.py` writes the `.npz`. The Zaretzki data still never ships.

### What survives from the ATTNSOM work

The licence diligence; the leakage finding and its grouped-split harness
(`scripts/som/leakage_check.py`); `bench_table.py` and the benchmark environment; the per-isoform
reporting; the grouped-numbers-only rule; the chemistry sanity check as a method; and the
applicability-domain index (Morgan atom-environment fingerprints), which is now cheap beside the
model rather than dwarfed by it.

**Dropped, and `tests/metabolism/som_golden.json` with them**: `graphnet.py`, the `SAGPooling`
reimplementation, the float16 weight archive, `convert_weights.py` and `golden.py`. The trained
ATTNSOM checkpoints stay in WSL as the baseline our own GBM must match, and the sections below are
kept as the record of how the decision was reached.

### M2 is done (2026-11-18): the package, the model and its measured accuracy

`chemsuite/metabolism/` is written, trained and tested, and the whole of it took 24 seconds of
CPU rather than a GPU afternoon — which is the pivot's dividend.

| file | what it is |
|---|---|
| `__init__.py` | `ISOFORMS` (the one-hot's column order, and so part of the archive's format) and `COMMON_ISOFORMS` |
| `features.py` | the 43 columns, `atom_block`, `design_matrix`, and `atom_environments` for the domain |
| `trees.py` | `Ensemble` and `probabilities`: the exported trees walked in numpy |
| `sites.py` | `MetaboliteSite`, `predict_sites`, `top_site`; loads the archive once per process |
| `predict.py` | `predict_metabolism`, the chunked cancellable driver |
| `data/som_gbm.npz` | **250 KB**: 300 trees, 18,300 nodes, and the domain index |
| `data/LICENSE-SOMBENCH.txt` | som-benchmark's MIT text, with a note on what was reused |

**The features are the app's, and training imports them.** `scripts/som/train_gbm.py` builds its
design matrix from `chemsuite.metabolism.features`, so train and serve cannot drift. The port was
checked against `som-benchmark`'s own `basic_atom_block` + `charge_block` over 22 molecules taken across the
curated set: **identical to 0.0 in every column**, names and order included.

**The numpy evaluation is exact, not close.** `scripts/som/export_gbm.py` refuses to write an
archive whose predictions differ from scikit-learn's `predict_proba` by more than 1e-9; the
measured difference for the shipped model is **1.1e-16** over all 41,878 training rows. So
`trees.py` is not an approximation of the estimator, it *is* the estimator, and no pickle and no
scikit-learn version pin are needed at run time.

**Accuracy, molecule-grouped, `PRIMARY_SOM` only, 10 seeds** (20% test / 10% validation of
molecules, tuned on seed 0's validation split by average precision; chosen `learning_rate` 0.05,
`max_leaf_nodes` 31):

| metric | mean | sd |
|---|---|---|
| MCC | 0.475 | 0.041 |
| top-1 | **0.592** | 0.039 |
| top-2 | **0.727** | 0.049 |
| top-3 | **0.813** | 0.047 |

| isoform | 1A2 | 2A6 | 2B6 | 2C8 | 2C9 | 2C19 | 2D6 | 2E1 | 3A4 |
|---|---|---|---|---|---|---|---|---|---|
| MCC | 0.476 | 0.459 | 0.503 | 0.552 | 0.470 | 0.585 | 0.444 | 0.404 | 0.434 |
| top-1 | 0.611 | 0.614 | 0.588 | 0.618 | 0.607 | 0.683 | 0.593 | 0.574 | 0.526 |
| top-2 | 0.742 | 0.779 | 0.707 | 0.753 | 0.727 | 0.817 | 0.722 | 0.724 | 0.670 |

That lands where the benchmark said it would: `basic_gbm` 0.601 top-1 in som-benchmark's own
10-seed run, ATTNSOM 0.618 grouped. **Ours is 0.592 — within 0.7 of a seed standard deviation of
both**, so the tie the decision was taken on holds for our own implementation of it. 3A4 is again
the weakest isoform and 2C19 the strongest, as in every run so far.

**What the user guide may say**, with the label definition beside it: for a molecule the model has
not seen, the true primary site is the top-ranked atom about **59%** of the time, in the top two
**73%**, and in the top three **81%**, per isoform numbers as above. No comparison against
SMARTCyp, XenoSite or FAME, since those are published on the inflated instance-split basis.

**The production model is fitted on all 2,003 instances** with the tuned settings, which is the
point of a grouped-CV estimate: the quoted numbers describe a model trained on 90% of the
molecules, and the shipped one has seen 10% more than that.

**The applicability domain ships with it.** `atom_environments` reads each atom's Morgan
environment identifiers (radii 0–2, `AdditionalOutput.GetBitInfoMap`) and the archive carries the
sorted unique identifiers seen in training, so `MetaboliteSite.domain` is the share of an atom's
environments the training set held: 1.0 for a drug-like amide, below 1.0 for a silane. RDKit
reports an environment only where it adds atoms to the smaller one, which is why the index is
`(radius, identifier)` pairs rather than one identifier per radius.

**Three decisions taken while writing it**, each because the alternative was quietly wrong:

* **`RemoveAllHs`, not `RemoveHs`.** Training molecules came from SMILES with implicit hydrogens
  throughout. `RemoveHs` keeps a hydrogen that carries an isotope, a charge or a stereo marker,
  and a hydrogen left in the graph changes its neighbour's `heavy_degree` — a feature. Atom
  indices come back in the caller's numbering either way.
* **Ties are broken by the lower atom index**, in `sites.py` and in the training script's top-k
  metric alike, so the measured number and the displayed one rank symmetry-equivalent atoms the
  same way. (Ibuprofen's two isobutyl methyls score identically to the last bit.)
* **A blank cell is not a failure.** `predict.py`'s `is_blank` treats `None`, `b""` and
  whitespace as "nothing here" and only counts a cell that held something unparsable, which is
  what the dialog's failure count should mean.

**Chemistry spot-check, in the test suite rather than a notebook**
(`tests/metabolism/test_sites.py`): testosterone → C6, allylic to the enone, 3A4's
6β-hydroxylation (score 0.59, next atom 0.04); paracetamol → the amide nitrogen, the route to
NAPQI (0.57); caffeine → C8, oxidized to 1,3,7-trimethyluric acid (0.27); ethanol → the carbinol
carbon for 2E1 (0.74). All four are the top-ranked atom.

**35 tests** in `tests/metabolism/` cover the column contract, the split rule and the missing-value
flag of the trees, the scikit-learn parity of a freshly fitted model, the chemistry above, the
hydrogen handling, the domain, the blank and failure rows, progress and cancel.
`tests/metabolism/som_golden.json` (36 KB) pins every isoform's scores for 15 molecules, written by
`scripts/som/golden.py`.

**Still to do:** M4 (the tool), M5 (showing the sites on the structure) and M6 (the user guide).

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
