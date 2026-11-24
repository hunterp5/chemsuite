# Auto QSAR/QSPR (plan approved 2026-11-23)

chemsuite builds its own QSAR/QSPR models: a large feature pool, automatic feature selection
and model selection inside nested cross-validation, honest prediction intervals, and a model
file that applies in another session and refines when new data arrives. It runs on numpy,
scipy, scikit-learn and RDKit, which are already base dependencies. **No new dependency, and
no torch.**

`docs/ROADMAP.md` §6 Tier 1 carries "QSAR, Auto QSAR" (added 2026-11-17) as a port of
MolManager's `qsar_models.py` and `auto_qsar*.py` (~2.9k lines). This plan replaces that port
with chemsuite's own design, on the pattern the pKa and site-of-metabolism models already
proved: fit with scikit-learn, **apply in numpy**, never load a pickle.

## Why

- **It is the tool that closes the loop.** chemsuite computes descriptors, fingerprints,
  conformers and pKa, and it plots and clusters data. "Given what we have measured, what
  should we make next" is the question all of that is for, and nothing in the app answers it.
- **Automatic is the honest kind.** A chemist choosing features and a model by hand, scored on
  the same data, overfits. Selection and tuning belong inside the cross-validation folds, and
  that is something code does correctly or not at all.
- **A model is only useful if it leaves the session.** Scoring an idea next month, or a
  colleague's library, means a model file — and refitting it as the project's data grows.

## Alternatives, with their licenses

| Tool | License | Verdict |
|---|---|---|
| **QSARtuna** (MolecularAI, formerly Qptuna; J. Cheminform. 2024) | Apache-2.0 | **Design reference, not a dependency.** Its `pyproject.toml` pins `chemprop 1.5.2` (and so torch) and `scikit-learn 1.4.0`, against our `>=1.5`. It is built around config files and a CLI, with mlflow and apischema, not a library a GUI drives. Its structure — a search over a model zoo, nested CV, uncertainty and applicability domain as first-class outputs — is the right shape, and is what this plan follows. |
| **QSPRpred** (CDDLeiden; J. Cheminform. 2024) | open source; RDKit + scikit-learn core, torch optional | Rejected. Lighter than QSARtuna, but it brings its own dataset, model and on-disk-store abstractions that duplicate `core/table`, `ToolContext` and document parts, and it would own our model file format. |
| **chemprop 2.3.1** | MIT | **Not in the app, for now** (decided with the user 2026-11-23). It requires `torch>=2.1`, `lightning>=2.0`, `pandas`, `astartes[molecules]`, `descriptastorus`, `rich`, `ConfigArgParse`, `cuik_molmaker_pin`, `myerson`. torch 2.15.0 wheels: 127 MB win_amd64, 539 MB manylinux x86-64 **plus** pinned `cuda-toolkit==13.2.2`, cudnn, cusparselt, nccl, nvshmem and `triton~=3.9`, and 139 MB macOS arm64 — multiple GB installed on Linux. Unlike the pKa and SoM models, the user fits this model on their own machine, so there is no train-outside-and-ship-numpy escape hatch: torch would have to live in the app's interpreter. Qt and torch also ship competing OpenMP runtimes on Windows, a known import-time crash. And `docs/SOM_PLAN.md` already measured a graph network at 8.5M parameters as **level** with a GBM over 34 atom descriptors. |
| **Optuna** | MIT | Rejected as a dependency. TPE buys little over random search plus halving at the budgets a desktop fit allows, and we need our own loop anyway for progress and cancellation (below). About six extra wheels for a second code path. |
| **ONNX via skl2onnx** | Apache-2.0 | Not now. A safe-to-load interchange format, worth revisiting only if someone needs to score outside Python. |

**Revisit chemprop as an optional external engine** — a path in Settings, like
`CHEMSUITE_OPENMM_PYTHON` and AiZynthFinder, run by subprocess over a JSON contract and shown
in the same leaderboard — only if real project data shows the in-app zoo falling short.

## Decided with the user, 2026-11-23

| Question | Decision |
|---|---|
| QSARtuna or QSPRpred as a dependency | No; QSARtuna's paper is the design reference. |
| chemprop in the app | Not at all, for now. |
| Hyperparameter search | Our own: random sampling of curated spaces plus successive halving. |
| Model storage | A `.csqm` numpy archive and a document part, plus a **write-only** pickle export. The app never unpickles a model. |
| First phase | Regression and classification together, the task detected from the target column. |

---

## The engine: `chemsuite/qsar/`

Qt-free, and it never imports the app (`tests/test_import_boundaries.py`). Mirrored by
`tests/qsar/`. Every module stays well under the 1,000-line budget.

### `recipe.py` — what makes a model reproducible

A `Recipe` dataclass of plain data with a JSON round-trip: the structure column, the
descriptor keys, the fingerprint variants, the extra numeric table columns, the preprocessing
choices, the split, the seed and the budget. A fit stores its recipe; apply and refine replay
it. Nothing about a model is implicit.

### `features.py` — the pool

Builds the design matrix from a `Recipe`, reusing what exists:

- descriptors through `chemsuite/descriptors/catalog.py` (`CATALOG`, `GROUPS`, `COMMON`,
  `in_group`) and `compute_descriptors` (`descriptors/compute.py`);
- fingerprints through `compute_fingerprints` and `compute_counts`
  (`analysis/fingerprints.py`), with `variant()` for radius and size — Morgan bits and counts,
  FCFP, RDKit, atom pairs, torsions, MACCS;
- any numeric column the user already has in the table.

It returns `X`, the feature names, and the mask of rows it could featurize. Text-valued
descriptors (`Descriptor.is_text`) are skipped; the pKa descriptors (`PKA_KEYS`) are available
but off by default, since each one costs a model run per row.

The default pool is deliberately large — Morgan counts at 2048 bits plus every numeric
descriptor group is of order 2,500 columns on a few hundred rows. Selection, not restraint, is
what makes that safe.

### `splits.py` — how the data is divided

Random, **scaffold** (`descriptors/scaffold.murcko`), **cluster** (Butina through
`analysis/clustering.butina`), and **by column** (a date or series column, for a
train-on-the-past split). Stratified variants for classification. It generates nested outer and
inner folds, and never puts the same scaffold or cluster on both sides.

Scaffold is the default: a random split on a congeneric series flatters every model.

### `select.py` — the cascade

Pure functions, each **fitted on a training fold only**, returning the kept indices plus the
statistics to report:

1. drop constant and near-zero-variance columns;
2. prune correlated pairs above a threshold, keeping the one more correlated with y;
3. a univariate filter — F-test for regression, mutual information as an option;
4. embedded importance — an L1 model or a forest's importances — down to a target count.

### `zoo.py` — the candidates

scikit-learn only. Each entry declares its search space, whether it needs scaling, and a cost
class the budget uses.

| Candidate | Space | Scaled |
|---|---|---|
| Ridge | `alpha` log-uniform 1e-3…1e3 | yes |
| ElasticNet | `alpha` log 1e-4…1e1, `l1_ratio` 0.05…0.95 | yes |
| PLS (`PLSRegression`) | `n_components` 1…min(20, features, n−1) | yes |
| SVR / SVC (RBF) | `C` log 1e-2…1e3, `gamma` log 1e-4…1e0 or `scale`, `epsilon` 0.01…0.5 | yes |
| RandomForest | `max_features` {0.1, 0.3, sqrt, 1.0}, `min_samples_leaf` {1, 2, 5}, `max_depth` {None, 10, 20}; 500 trees fixed | no |
| HistGradientBoosting | `learning_rate` log 0.01…0.3, `max_leaf_nodes` {15, 31, 63}, `min_samples_leaf` {5, 10, 20}, `l2_regularization` 0…1; `max_iter` 500 with early stopping | no |
| kNN | `n_neighbors` 1…15, `weights` {uniform, distance}; Jaccard on fingerprint bits, Euclidean otherwise | numeric only |
| MLP | `hidden_layer_sizes` {(64,), (128,), (256, 64)}, `alpha` log 1e-5…1e-1, `learning_rate_init` log, `early_stopping` | yes |

Classification mirrors these, with `LogisticRegression` in place of Ridge and ElasticNet. SVC
uses `decision_function`, not Platt scaling: conformal calibration needs a score, not a
calibrated probability, and `probability=True` costs an internal 5-fold fit.

### `search.py` — the budget

Our own nested loop, not `HalvingRandomSearchCV` (which needs
`sklearn.experimental.enable_halving_search_cv` and gives us no hook for progress or cancel):

- **outer folds** give the honest estimate. For each outer fold, the whole search runs on its
  training part only, the winner is fitted there, and it is scored on the held-out part.
- **inner folds** choose. Per candidate, sample N configurations from its space, evaluate them
  on a fraction of the inner folds, keep the best half, widen, repeat. Preprocessing and the
  whole of `select.py` are fitted inside each inner training fold.
- **the shipped model** is the same search run once over all the data and fitted on all of it.
  The reported scores are the outer-fold ones, never the shipped model's own fit.

Three budgets in the dialog — Quick, Standard, Thorough — set the candidate list, the
configurations per candidate and the fold counts. Seeded throughout: the same seed, data and
recipe give the same model. Progress is reported per candidate configuration, and cancel is
checked between them.

### `metrics.py`

Regression: R², Q² (out-of-fold), RMSE, MAE, Spearman. Classification: ROC-AUC, PR-AUC, MCC,
balanced accuracy, calibration error. Default objective: RMSE for regression, ROC-AUC for two
classes, MCC above two.

### `conformal.py` — intervals and the applicability domain

- **Regression:** cross-conformal intervals from the outer folds' out-of-fold absolute
  residuals, normalized by a difficulty estimate (a kNN residual model), so the interval is
  wider where the model is weaker. At a chosen confidence, 90% by default.
- **Classification:** Mondrian prediction sets — a per-class quantile of out-of-fold
  nonconformity, so a row's set is empty, one class, or several.
- **Applicability domain:** the mean Tanimoto distance to the k nearest training compounds
  (fingerprints), or the leverage on the selected descriptor block, with "outside" set at the
  95th percentile of the training set's own distances. One number and one flag per row.

A conformal interval is the honest answer to "how much do I trust this prediction", and it
costs no extra fitting: the out-of-fold predictions the outer loop already produced are exactly
its calibration set.

### `trees.py` — its own walker, and why not the one that exists

A general walker over flat node arrays: forests as well as boosting, regression as well as
classification, several values in a leaf, and a per-tree output column for multiclass
boosting.

**`chemsuite/metabolism/trees.py` stays where it is** — decided with the user on 2026-11-24,
against this plan's first instinct of moving it. That walker is specialized to the
site-of-metabolism model: a binary classifier, a scalar baseline, one log-odds value per leaf
and a logistic on top. It is also pinned by `scripts/som/golden.py` and the shipped
`som_gbm.npz`, neither of which gains anything from being generalized, and making
`metabolism` import from `qsar` would couple two unrelated features the wrong way round. The
price is about 25 lines of node descent living in two places, and each docstring points at the
other.

### `export.py` and `apply.py` — fit with scikit-learn, apply in numpy

`export.py` turns a fitted estimator into flat arrays, by family:

| Family | Arrays |
|---|---|
| Ridge, ElasticNet, LogisticRegression, PLS | coefficients and intercept |
| SVR, SVC | support vectors, dual coefficients, intercept, gamma |
| RandomForest, HistGradientBoosting | node arrays, as `trees.py` already reads |
| kNN | the training matrix and its targets |
| MLP | the weight matrices, the biases and the activation name |

**`export.py` refuses to write when `apply.py`'s prediction differs from scikit-learn's by more
than 1e-9** — the rule `scripts/som/export_gbm.py` already follows, and what makes the numpy
path the estimator rather than an approximation of it. `apply.py` imports numpy and nothing
else: scikit-learn is never imported to apply a model.

### `model.py` — the `.csqm` archive

An uncompressed `.npz` on the `macromolecule/structure/store.py` pattern: a `_format` integer
array, a `_meta` JSON blob as `uint8`, and plain numeric or fixed-width text arrays beside
them. Read with `np.load(..., allow_pickle=False)`, and a `ModelError` for anything else. The
manifest holds:

| Key | Holds |
|---|---|
| `format`, `created`, `chemsuite` | the format number, the date, the version that wrote it |
| `task`, `target` | regression or classification; the target column and its class labels |
| `recipe` | the `Recipe`, verbatim |
| `preprocess` | imputation values, centers, scales |
| `features` | the names, and the indices the cascade kept |
| `estimator` | family, hyperparameters, the names of its arrays |
| `scores` | the outer-fold metrics: the leaderboard, and the winner's |
| `conformal` | confidence, the calibration quantiles, the difficulty model |
| `domain` | kind, k, the threshold |
| `training` | n, the row keys (InChIKey), y, and the weights |
| `lineage` | one entry per fit or refit: date, n, scores, note |

Training keys and y are stored because refinement needs them: a model refines by refitting its
recipe on the old rows plus the new ones, not by mutating a frozen object.

### `refine.py`

Replay the recipe on the old rows plus the new, either re-tuning lightly around the chosen
hyperparameters or re-running the full search, then append a lineage entry. A row already in
the training set is matched by key and its y updated, not duplicated.

---

## The tool: `chemsuite/tools/qsar/`

Data ▸ QSAR ▸ **Auto QSAR…**, **Apply Model…**, **Refine Model…**, with the model panel
reopening through its `restore`, following `chemsuite/tools/space/tool.py`.

- `tool.py` — the `Tool` entries, `requires=("rdkit", "sklearn", "scipy")`.
- `dialog.py` — a kit `ToolDialog` (`app/widgets/dialogs/kit.py`), built like
  `tools/metabolism/dialog.py`: **Input** — structures, the target column, `RowsChoice`;
  **Method** — the detected task shown, the split, the budget; **More Options** (collapsed) —
  the feature pool, the confidence, the seed, the selection thresholds; **Output** — the
  columns to write, as checkboxes. No helper text anywhere: explanations are
  `Section.row(..., help=...)` tooltips, and the only text shown is a blocking `problem()`
  ("There is no numeric column to model", "Fewer than 20 rows have a value").
- `panel.py` — the model panel: the leaderboard, observed against predicted and residuals in
  pyqtgraph, feature importance, the lineage, and the actions Apply to Table, Refine with New
  Data…, Save Model…, Export for Python….
- `part.py` — the document part under key `qsar_models` (`ToolContext.document_part`,
  `app/document/workspace.py`): each fitted model's archive bytes and the panel state, with a
  `revision` so the document knows it needs saving. A model survives closing the panel.
- `apply_dialog.py`, `refine_dialog.py` — choose a model (a `.csqm` file, or one already in the
  document), the rows, and the columns.

Output columns, one set per run, suffixed `(2)`, `(3)` when a name is taken: **Prediction**,
**Residual** (where y is known), **Interval Low** and **Interval High** — or **Predicted
Class** and **Class Set** — **Domain Distance**, **In Domain**, and **Fold** (which outer fold
held the row out, so a plot can show honest predictions).

The fit runs through `ctx.run` as a background job. It is a whole-dataset result, so **no**
`keep_partial` on cancel.

### The pickle export

Export for Python… writes `model.joblib` (joblib is already a scikit-learn dependency) plus a
sibling manifest naming the Python, numpy and scikit-learn versions it was written under, for
use in a notebook. **The app never reads it**: Apply Model… takes `.csqm` only, and says so.
Both formats embed training data — the user guide says that in one sentence.

---

## What the code settled (2026-11-24, built)

- **One basis per fold.** Filling and selection depend on the training rows and the recipe,
  not on the candidate or its settings, so they are computed once per fold
  (`search.basis_for`) and shared by every configuration tried there. Done naively, a Quick
  run refits a 100-tree selection forest for every configuration; this changes no number,
  because the basis is still fitted on that fold's training rows alone.
- **The conformal quantile is an order statistic**, the `ceil((n + 1) · confidence)`-th
  smallest nonconformity — not `np.quantile`, whose interpolation gives the finite-sample
  correction away. The tests pin both small-sample cases.
- **A row is never its own neighbor** in the local difficulty or the domain distance
  (`_nearest(..., skip=1)`), or every training row would sit at distance zero.
- **Fewer than `MIN_ROWS` (20) rows with a value is refused** by `search.run` with a reason,
  rather than failing later inside a fold.
- **The leaderboard's settings come from the final search over all the rows; its scores are
  out of fold.** The headline score is always the procedure's, because the leaderboard's
  winner won a comparison on those same folds.
- **Cost, measured:** a Quick run (4 candidates, 8 configurations, 3 outer × 3 inner folds) on
  160 rows and 43 columns takes **28 s**, of which the 500-tree random forest is 22 s; without
  the forest and the neural network, 6 s. Standard will be minutes. **Next:** run the folds on
  `ctx.process_pool` so a fit uses the machine it is on — the fits are independent, so this is
  a scheduling change, not a numerical one.
- **`SVC` is built without naming `probability`:** scikit-learn deprecated the parameter in
  1.9, and passing it even as False warns. Conformal calibration wants the decision value.

Then, from step 4:

- **Distances are taken directly, never by the expanded formula.** `‖a‖² − 2a·b + ‖b‖²` is one
  BLAS pass and cancels catastrophically: a row measured against itself comes out around 1e-8
  instead of 0, and scikit-learn's inverse-distance weighting gives a *coincident* training row
  all of the weight, so 1e-8 silently takes the other branch. The kNN export missed
  scikit-learn by 1.9e-07 until this was fixed — caught by the 1e-9 check, which is exactly
  what it is for. `chemsuite/qsar/neighbors.py` now holds one tiled, accurate implementation,
  used by both the neighborhood model and the conformal/domain distances.
- **PLS is exported as a linear model.** `PLSRegression.predict` only *centers* its columns —
  "the coefficients are already scaled" — so the centering folds into the intercept and what
  is left is `x @ coefᵀ + b`, exact to 9e-16. There is no separate PLS family to apply.
- **A random forest's leaf values are normalized on the way out**, because scikit-learn has
  changed between storing class counts and class fractions.
- **A support vector machine is exported for two classes only**, and `export` refuses above
  that rather than approximating: scikit-learn's multiclass margin is a vote between pairs.
  The zoo must leave `svc` out when a target has more than two classes.
- **No tree in the zoo carries a missing-value rule**, because the pool's holes are filled
  before any estimator sees them, so `missing_left` is written as zeros for forests.
- **`fit.py` is the only module that knows all the others.** The apply path — `model`,
  `apply`, `trees`, `prepare`, `neighbors`, `conformal` — never imports the fitting side, and
  a test applies a saved model in a fresh interpreter with `sklearn` blocked at the import
  hook to prove it.

Then, from steps 5 and 6 — the tool:

- **A model is a document asset, not part state.** `tools/qsar/part.py` keeps only asset ids
  and the panel's layout; each model's bytes go in the document's asset store, where proteins
  live. That buys three things the part could not give on its own: adding a model is an op, so
  a fit undoes with its columns in one step; an undone fit's entry disappears from the list
  because its asset has gone; and `Model` satisfies `core.assets.Packable`, so the store hands
  back the decoded model and a session unpacks it once rather than once per prediction.
- **The columns a fit writes are out of fold, and that is a design decision, not a detail.**
  `Fitted.held` carries the predictions the search already made, each from a fold that did not
  contain its own row. Re-predicting with the finished model would have been one line shorter
  and would have quietly replaced honest numbers with flattering ones, so `fit.py` builds the
  `Prediction` and the tool never computes one for training rows.
- **`ToolContext.texts()` was added** rather than widening the ratchet. Reading a label column
  needed `table.column(...).displays(table.slots_of(...))`, which
  `tests/test_architecture.py::test_tools_use_the_table_through_the_context` refuses in a tool
  — correctly. The helper is the general one the rule asks for: cell texts for oids, beside
  `plotting.data.numbers` for floats.
- **The zoo is asked what it can export before the search, not after.** `for_task(...,
  classes=n)` drops a candidate whose arrays cannot be written for that many classes, so a
  three-class target is never offered support vectors. Without it, a Thorough fit could search
  for minutes and fail at the final export, which is the worst possible moment.
- **Asset ids are plain ASCII.** A target column can be called `IC50 (µM) / 24 h`, and that id
  becomes a name inside the saved document.
- **Identity is the InChIKey**, computed by `features.keys_for`. A refit matches on it, so a
  remeasured compound replaces its old value instead of being learned twice; a different
  tautomer or salt is a different compound, which is the conservative answer.
- **A refit is a refit.** `refine.py` relearns the fill values, the features, the settings and
  the calibration on the union, because a model nudged by new rows keeps a selection chosen on
  the old ones and a score measured against folds the new rows were never in. After a few
  rounds nobody could say what the number meant. One lineage line per round, not two.
- **Export for Python is write-only and session-only.** The pickle needs the scikit-learn
  object, which only exists while the document that fitted it is open; a reopened model is
  arrays. The dialog says so rather than refitting something slightly different behind the
  user's back, and nothing in chemsuite ever reads a pickle.

## Steps

1. **This document**, and the decisions in CLAUDE.md's **Decided** section in short form (no
   chemprop or torch in the app; our own search; a model archive is never unpickled), since the
   project moves between two machines and memory does not travel.
2. **Engine core** — `recipe.py`, `features.py`, `splits.py`, `select.py`, `zoo.py`,
   `search.py`, `metrics.py`, with their tests.
3. **Uncertainty and domain** — `conformal.py`, with coverage tests.
4. **Archive and apply** — `trees.py` (its own, general), then `export.py`, `apply.py`,
   `model.py`, the 1e-9 verification and the round trip.
5. ✅ **The tool** — the package above, the document part and the panel, with
   `docs/USER_GUIDE.md` and the `ROADMAP.md` entry in the same commit.
6. ✅ **Refinement** — `refine.py`, the lineage, the Refine dialog.

CoMSIA (ROADMAP §6) builds on this later. chemprop as an external engine stays deferred.

## Tests

Each one fails without the code it covers.

- **Recovery:** y built as a known function of a few descriptors — the cascade selects those
  features and the outer folds reach high R².
- **Leakage:** y shuffled — outer-fold R² ≈ 0, ROC-AUC ≈ 0.5. This is the test that proves
  selection and tuning sit inside the folds. It mirrors `scripts/som/leakage_check.py`.
- **Export exactness:** every family, scikit-learn against `apply.py`, within 1e-9, and
  `export.py` refusing to write when it is not.
- **Round trip:** fit → `.csqm` → a fresh process → identical predictions, with
  `allow_pickle=False` throughout; and a `.csqm` carrying a pickled object raises `ModelError`.
- **Conformal coverage:** a 90% interval covers about 90% of held-out rows.
- **Splits:** no scaffold and no cluster on both sides.
- **Document:** fit, save, reopen — the panel, the leaderboard and the model come back.
- **Refine:** added rows keep the recipe, extend the lineage, and do not duplicate a row
  already trained on.

Targeted runs per CLAUDE.md: `pytest tests/qsar`, `pytest tests/tools/qsar`, plus
`tests/test_architecture.py` and `tests/test_import_boundaries.py` between steps. The full
suite runs once, before the `dev`→`main` merge.
