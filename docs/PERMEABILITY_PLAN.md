# Permeability: two shipped models on chemsuite's own QSAR engine

Decided and measured 2026-11-25. **Two endpoints ship** — the MDR1-MDCK efflux ratio and Caco-2
log Papp — each as a `.csqm` archive fitted by `chemsuite/qsar` and applied in numpy, from
descriptors alone. PAMPA was fitted and held back, and the 3D arm was measured and closed; see
[What the numbers decided](#what-the-numbers-decided).

## Why this is not a new modelling package

The ROADMAP listed permeability under *Tier 3: heavy add-ons*, assuming chemprop and torch —
about 2 GB of install. Both halves of that assumption have since stopped being true.

**The field moved to where chemsuite already is.** The February 2026 critical assessment of the
TDC ADMET leaderboards found data leakage in several top entries and concluded that classical
ML is the most stable family across endpoints; the April 2026 *Do Larger Models Really Win in
Drug Discovery?* benchmark agreed. The current leader on the TDC Caco-2 leaderboard —
CaliciBoost, MAE 0.256 log units — is AutoML over PaDEL/Mordred/RDKit descriptors, with no
pretraining and no graph network. Against an assay whose inter-laboratory scatter is 0.3–0.5
log units, 0.26 is already the noise floor: architecture is no longer where the gain is.
Applicability domain is. The November 2026 natural-products benchmark makes that point
bluntly — three commercial ADMET predictors correlated r = 0.05–0.10 with measurement once
they were off their training chemistry.

**chemsuite already had the engine.** `chemsuite/qsar` is a complete Auto QSAR system: a
feature pool over the descriptor catalog and fingerprints (`qsar/features.py`), a model zoo
(`qsar/zoo.py`), nested-CV search with scaffold and cluster splits (`qsar/search.py`,
`qsar/splits.py`), conformal intervals and an applicability domain (`qsar/conformal.py`), and
a `.csqm` archive applied in pure numpy (`qsar/apply.py`, `qsar/model.py`) whose export refuses
to write unless numpy reproduces scikit-learn to 1e-9.

So the state of the art for this endpoint **is** what chemsuite can already fit, and this work
is curated data, three fitted `.csqm` archives shipped in-tree, and a thin tool that applies
them. Permeability moves out of Tier 3 into the **base install**: no new runtime dependency, no
torch, nothing downloaded at run time, and no pickle ever loaded.

chemprop and ADMET-AI in an external interpreter are **dropped**, consistent with the standing
decision that chemprop and torch stay out of the app.

## The endpoints, and the data behind them

Two ship; the third was curated, fitted and held back.

| Endpoint | Ships | Task | Unit | Source | Curated rows | Licence |
|---|---|---|---|---|---|---|
| `mdck_er` | **yes** | regression | log10 (B-A/A-B) | Fang et al. 2023, Biogen public ADME set, commit `b00df00` | **2,638** | **MIT**, shipped beside the model |
| `caco2` | **yes** | regression | log10 cm/s | Wang et al. 2016, via TDC `caco2_wang` | **893** | paper supporting information redistributed by TDC; cited in the NOTICE |
| `pampa` | no | classification | permeable / not | NCATS PAMPA pH 5, via TDC `pampa_ncats` | **2,034** (85% permeable) | US government work, not under copyright |

**PAMPA is a classification, and that was not the original intent.** The plan assumed the
~6,500-compound continuous NCATS set. What TDC actually serves (`pampa_ncats`) is 2,035
compounds with a binary label, and ADME@NCATS publishes its data through a JavaScript
application with no pinnable file. Rather than scrape it, the endpoint ships as a
classification — which the engine supports natively, Mondrian class sets and all, and which
suits an 85/15 split better than a regression would. Folding in a continuous PAMPA set later
is a data change, not a code change.

**Caco-2 is the small one**, and the obvious way to grow it is the Esaki 2024 ChEMBL curation
(doi:10.1186/s13321-024-00826-z), which kept 1,661 entries of a much larger pull by excluding
ambiguous cell lines, P-gp/BCRP overexpression and basolateral-to-apical measurements. Those
exclusions are not applied today because none of the three current sources needs them: each is
a single named assay with its own column. They become mandatory the day that set is folded in,
and `scripts/permeability/curate.py` says so where the rule would go.

## Curation

`scripts/permeability/fetch.py` pins every file by SHA-256 and line count and refuses a file
that has changed upstream, as `scripts/cyp/fetch.py` does; it fetches through
`chemsuite.core.download.fetch_bytes`, so it inherits the app's certificate handling, its
HTTPS-only rule and its audit entry. `curate.py` then, per endpoint:

- standardizes with chemsuite's own standardizer — clean up, strip salts, neutralize, but
  **no** tautomer canonicalization and **no** protonation — so a feature means the same number
  the table's own descriptor columns would give;
- deduplicates by InChIKey, averages readings within the endpoint's agreement window (0.5 log,
  the assay's own reproducibility) and drops wider pairs, counting them;
- **refuses to continue** if more than 5% of a regression endpoint's values sit on one bound,
  because a censored target cannot be scored honestly as a regression.

Measured when this was written: no endpoint is censored (one value at each extreme for both
regression endpoints), 3 structures failed to standardize across all three sources, and
duplicate merging removed 1, 14 and 1 compound respectively.

Neither the snapshot nor the curated CSVs are committed; `samples/` is local test data.

## The bake-off

`scripts/permeability/train.py` fits through the app's own engine and reports the engine's
**out-of-fold** scores under a scaffold split — the only honest figure available without a
second holdout. Four arms per endpoint, same split, same seed, same budget, differing only in
the feature pool:

| Arm | Pool | The question it answers |
|---|---|---|
| `plain` | every descriptor needing no conformer or model, plus Morgan counts (2,151 columns) | the baseline, and what the evidence says should win |
| `pka` | `plain` + the pKa-derived descriptors (`cLogD 7.4` and friends) | worth seconds per molecule instead of milliseconds? Ships only if it wins by ≥ **0.05 log MAE** |
| `descriptors` | `plain` without the fingerprint block (103 columns) | do 2,048 sparse columns earn their width on sets this size? |
| `shape` | `plain` + the 3D shape descriptors over one ETKDG conformer | the cheap proxy for the 3D question, and the gate on the Uni-Mol arm |

Arms are compared at the `quick` budget and the shipped model is fitted at `standard`, because
a comparison is only meaningful when every arm is given the same search.

### Measured results

Out-of-fold, scaffold split, `quick` budget, seed 0, on this machine (2026-11-25). Every figure
is the engine's own `Fitted.scores`; the winning candidate is named because it changed between
endpoints.

**MDR1-MDCK efflux ratio** (2,638 compounds, log10 ratio):

| Arm | Candidate | MAE | RMSE | R² | Features kept | Minutes |
|---|---|---|---|---|---|---|
| `plain` | hgb | 0.380 | 0.494 | 0.484 | 200 / 2,151 | 9.2 |
| `descriptors` | hgb | 0.370 | 0.485 | 0.504 | 92 / 103 | 5.3 |
| `shape` | hgb | **0.364** | **0.483** | **0.506** | 200 / 2,162 | 11.3 |

**Caco-2 permeability** (893 compounds, log10 cm/s):

| Arm | Candidate | MAE | RMSE | R² | Features kept | Minutes |
|---|---|---|---|---|---|---|
| `plain` | rf | 0.371 | 0.475 | 0.622 | 200 / 2,151 | 4.3 |
| `descriptors` | svr | 0.372 | 0.475 | 0.622 | 88 / 103 | 2.4 |
| `shape` | rf | 0.371 | 0.475 | 0.622 | 200 / 2,162 | 3.8 |

**PAMPA** (2,034 compounds, binary, 85% permeable): `plain` / rf — ROC-AUC 0.740, balanced
accuracy 0.570, MCC 0.252, 200 / 2,151 features, 7.4 min. The other two arms were not run: the
endpoint was held back on this number, so measuring its pool further would have bought nothing.

### The models that ship

Refitted at the `standard` budget — the full zoo rather than its four cheapest entries — on the
`descriptors` pool the comparison chose. These are the numbers the user guide quotes:

| Endpoint | Candidate | MAE | RMSE | R² | Spearman | Features | Archive | Fit |
|---|---|---|---|---|---|---|---|---|
| `mdck_er` | SVR | **0.360** | 0.481 | 0.512 | 0.732 | 92 / 103 | 3.7 MB | 27.9 min |
| `caco2` | random forest | **0.354** | 0.454 | 0.655 | — | 88 / 103 | 8.8 MB | 20.0 min |

The bigger budget was worth it: MDCK improved 0.370 → 0.360 and Caco-2 0.372 → 0.354, and in
both cases a different estimator won than at `quick`. MDCK's shipped model also beats the
`shape` arm's best (0.364), so closing the 3D arm cost nothing at all.

Measured end to end through the tool's own path, one process and no worker pool:
**13 ms per row**, so a two-thousand-row table is about half a minute. That is the figure to
compare against the `pka` arm's ~60 CPU-seconds per molecule below.

### What the numbers decided

**Fingerprints are out.** On MDCK the 103 descriptors alone beat descriptors plus 2,048 Morgan
counts (0.370 against 0.380) in 60% of the time; on Caco-2 the two are identical to three
decimals. The 2,048 sparse columns give the selection cascade more noise to discard and leave it
with 200 features that are no better than the 92 it would otherwise keep. Two endpoints and
three different winning estimators agree, so this is not a fluke of one fit.

**The 3D arm is closed.** `shape` minus `plain` — the same pool plus 11 shape descriptors — is
+0.016 log on MDCK and 0.000 on Caco-2. Against an assay that scatters 0.3–0.5 log that is
nothing, and it is nowhere near the 15.7% MAE reduction CaliciBoost reported from 3D
descriptors. The gate in the last section therefore does not open: **no Uni-Mol head, no second
84 MB of weights, no WSL training run.** `shape` also needs an ETKDG conformer per row at
predict time, which would turn a millisecond-per-row tool into a second-per-row one for a tie.

**PAMPA is held back.** ROC-AUC 0.740 says the ranking carries signal, but balanced accuracy
0.570 against a 0.5 baseline, on data that is 85% permeable, says the model mostly predicts the
majority class. A per-compound call would need more caveat than it is worth. `curate.py` and
`train.py` still cover the endpoint and the classification path is still exercised by
`tests/permeability/test_predict.py`, so it is one line in `ENDPOINTS` to restore — the right
trigger being a continuous PAMPA set rather than a better model on this one.

**The `pka` arm was abandoned on measured cost, with its accuracy never obtained.** This is the
one place the bake-off did not run to a number, so it is stated plainly rather than tidied up.
The arm was started twice — on MDCK (2,638 molecules) and then on Caco-2 (893) — and neither
finished building its feature pool: MDCK consumed about 67,000 CPU-seconds in two hours and
Caco-2 about 50,000 in one, both still short of a complete pool on twelve cores. That puts the
pKa descriptors at **upwards of 60 CPU-seconds per molecule**, measured consistently on two
independent sets.

That number, not an MAE, decides the gate. The question this arm was set up to answer assumed
the cost was "seconds rather than milliseconds" per molecule; it is nearer a minute of CPU per
molecule, and **the same pass runs at predict time**, so a two-thousand-row table would cost
hours on a twelve-core machine for a column whose neighbors in Chemistry ▸ Calculate take
milliseconds. No MAE gain available on a 893-row set could make that shippable as the default,
so the comparison stopped being load-bearing and the compute was better spent on the models
that ship.

**The honest caveat:** this reverses the "let the bake-off decide" rule on cost grounds rather
than accuracy grounds, and it is a judgment rather than a measurement of what was asked. The
accuracy question is still open and still cheap to settle on an idle machine overnight:

```sh
.venv/Scripts/python.exe scripts/permeability/train.py --arm pka --endpoint caco2
```

If `pka` turns out to clear +0.05 log MAE, the right outcome is not to change the default but
to offer it as a second, slower model — which is new scope, not a tweak. Before rerunning it,
`train.py` should pass a `progress` callback into `features.build` so the pool's rate is
visible; its absence is why two runs had to be abandoned blind.

**Caco-2 is data-limited, not model-limited.** Three pools and three different winning
estimators give identical scores at the `quick` budget: nothing about the features moves it. A
bigger search does help a little (0.372 → 0.354 at `standard`), which bounds what tuning is
worth — the remaining distance to the published TDC leader on this same Wang data (0.256) is
not another pool or another budget away. The lever is more data, so folding in the Esaki 2024
ChEMBL curation (+1,661 entries, with its exclusions applied) is the recommended next step for
this endpoint. Part of the gap is also split and curation rather than model: ours deduplicates
to 893 compounds and uses its own scaffold split, and the February 2026 audit found leakage
inflating several leaderboard entries.

## Reproducing it

Every step below has been run on this machine, in this order:

```sh
.venv/Scripts/python.exe scripts/permeability/fetch.py
.venv/Scripts/python.exe scripts/permeability/curate.py
# the comparison, quick budget, same split and seed for every arm
.venv/Scripts/python.exe scripts/permeability/train.py --arm plain --arm descriptors --arm shape
# the models that ship, standard budget, descriptors only
.venv/Scripts/python.exe scripts/permeability/train.py --ship descriptors \
    --endpoint mdck_er --endpoint caco2
# the pinned predictions the golden test checks
.venv/Scripts/python.exe scripts/permeability/golden.py
```

`--ship` prints the `Pin(size, sha256)` line for each archive; those go into `PINS` in
`chemsuite/permeability/models.py` by hand, so a refit is never a silent reship.
`tests/permeability/test_models.py` then checks the shipped copies against them, and
`tests/permeability/test_golden.py` pins what they predict to 1e-7. Both skip in a checkout
without the archives, so a missing model is visible rather than silent.

The `pampa` endpoint is still curated and trainable; it is simply not in `ENDPOINTS`.

## The 3D arm: measured, and closed

`chemsuite/ionisation/unipka/model.py` is already a complete numpy forward pass of the Uni-Mol
encoder, shipped for pKa, so a permeability head would have been cheap in *code* — but not in
*bytes*: fine-tuning moves the backbone, so a second endpoint means a second ~84 MB file unless
the encoder is frozen and only a head trains. The published evidence was mixed (3D descriptors
cut MAE 15.7% in CaliciBoost; the 2026 PROTAC preprint saw no gain), so the plan made it
conditional rather than guessing.

**The gate did not open.** The `shape` arm — the same pool plus 11 shape descriptors over an
ETKDG conformer, which is the cheap proxy for what a 3D model would contribute — gained
**+0.016 log MAE on MDCK and 0.000 on Caco-2**. Against an assay that scatters 0.3–0.5 log
between laboratories that is indistinguishable from nothing, and it is nowhere near the 15.7%
that would have justified a transformer. On Caco-2 the three pools are identical to three
decimals, which says the limit there is 893 training rows, not the representation.

So: **no Uni-Mol head, no second 84 MB of weights, no WSL training run.** Reopening it would
need a new reason — most plausibly a much larger Caco-2 set after the Esaki merge, where a
representation could start to matter. Until then this is closed, not deferred.

## Sources

- Fang et al., *J. Chem. Inf. Model.* 2023, **63**, 3263–3274, doi:10.1021/acs.jcim.3c00160;
  data at `github.com/molecularinformatics/Computational-ADME` (MIT).
- Wang et al. 2016, Caco-2 permeability, via Therapeutics Data Commons `caco2_wang`.
- NCATS PAMPA pH 5, via Therapeutics Data Commons `pampa_ncats`.
- Esaki, Yonezawa & Ikeda, *J. Cheminform.* 2024, **16**, 32, doi:10.1186/s13321-024-00826-z.
- CaliciBoost, *J. Cheminform.* 2025, **17**, 184, doi:10.1186/s13321-025-01137-7.
- *Critical Assessment of ML models for ADMET Prediction in TDC leaderboards*, bioRxiv,
  February 2026.
- *Do Larger Models Really Win in Drug Discovery?*, bioRxiv, April 2026.
- *The Challenges of PROTAC Permeability Prediction*, arXiv 2610.07998.
