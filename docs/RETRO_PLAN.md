# Plan: retrosynthesis with AiZynthFinder (drafted 2026-11-06)

## 1. Goal and ground rules

A chemist picks rows and asks for routes: chemsuite searches each target with a **local**
AiZynthFinder install, writes the search statistics onto the row, and shows the routes as a
tree of precursor cards that can be read, copied and drawn in the sketcher.

- **Local only.** The policy models (ONNX), the templates and the stock file all come from a
  `config.yml` on this machine. A config that sets `use_remote_models: true`, or names a
  MongoDB stock on a host other than localhost, is refused before any process starts: nothing
  about a structure leaves the machine.
- **Its own interpreter.** AiZynthFinder 4.4.1 needs RDKit 2023.x and NumPy 1.26; chemsuite
  runs a much newer RDKit. So it is never imported into chemsuite's process. It runs in the
  interpreter named by `tool_paths.aizynthfinder_python`, already in Settings ▸ External Tools
  together with `aizynthfinder_config`.
- **Nothing vendored.** No AiZynthFinder code is copied into chemsuite. We talk to it through
  a small worker script and its own public API (`AiZynthFinder`, `prepare_tree`,
  `tree_search`, `build_routes`, `extract_statistics`, `stock_info`).
- **Presentation kept.** MolManager's route view is the one to match: a target list, a route
  table, a precursor tree whose cards are framed green (in stock) or amber (still needs a
  route), numbered Step buttons on the branch lines, and step metadata beside the canvas.

## 2. Where the speed comes from

MolManager spawned `aizynthcli` once per job. Every job therefore paid the model load
(the USPTO policy plus the 1.3 GB ZINC stock) before the first target, and paid it again for
the next job. Measured on this machine (M-series, three targets, 15 s limit each): **33.5 s
warm against 77.7 s** with a fresh process per job — the ~15 s load is paid once, not three
times, and a second job on the same config starts searching immediately.

chemsuite keeps **one worker process per interpreter+config, alive across jobs** and speaks a
line protocol to it:

| | MolManager | chemsuite |
|---|---|---|
| Models loaded | once per job (~15 s each time) | once per session (first search) |
| Progress | parsed from the checkpoint file | one reply per target |
| Cancel | kill the process, read the checkpoint | the finished targets are already in hand; the process stays warm |
| Several algorithms | a second script, its own process | the same warm worker, one request per algorithm |

The worker is a single `-c` script (no file written next to untrusted data), reads one JSON
request per line on stdin and writes one JSON reply per line on stdout, so stdout stays clean
even though AiZynthFinder logs to stderr. Requests: `hello` (version, policy/filter/stock
names), `search` (one target), `quit`.

## 3. Phases

| Phase | Work | Done when |
|---|---|---|
| A1 ✅ | Qt-free core, `chemsuite/retro/`: resolve the interpreter and config, refuse non-local configs, read the policy/filter/stock catalogue from `config.yml`, the worker script and its protocol, `SearchOptions`, parse replies into table columns and compact routes (score, solved, steps, nested tree, stock names, step metadata). | Done 2026-11-07: `tests/retro/test_retro_core.py` (40 tests). The parsers run on recorded AiZynthFinder output; the local-only rules on crafted configs; a fake interpreter that honours `-c` runs chemsuite's own worker script against a stub `aizynthfinder`, covering the protocol, the options reaching AiZynthFinder's settings, a target that fails mid-job, a library printing to stdout, a crash, a warm worker over several targets, and cancel keeping the finished targets. |
| A2 ✅ | Chemistry ▸ Retrosynthesis…: the dialog (rows from a structure column, policies/filters/stocks from the config, time and iteration limits, max routes, algorithm, exploration C, stop at first, cluster), the job, the `AiZynth …` columns, partial results on cancel, rows that already have values skipped. | Done 2026-11-07: `tests/app/test_retro_ui.py` drives the dialog against the fake interpreter: the config's policies and stocks offered, a missing install and a remote config blocking the run, the columns written, the options reaching the search, rows already searched skipped, and one structure failing without stopping the rest. |
| A3 ✅ | The Routes pane: target list, route table, precursor tree with green/amber cards, numbered step buttons whose metadata opens under the canvas, pan and zoom, copy a card or add it to the table, draw a step in the sketcher. Saved with the document (`tool_id`/`state`/`restore`). | Done 2026-11-07: the pane draws the target to the right of its precursors with green/amber cards, numbers the steps from the earliest reaction, opens a step's metadata, copies a card, and comes back from a saved session with the same target and route selected. |
| A4 ✅ | Settings ▸ External Tools: Download… for the public data (Zenodo: USPTO policy, ringbreaker, filter, ZINC stock, written `config.yml`), and Test that reports the AiZynthFinder version from the worker. User guide section. | Done 2026-11-07: `chemsuite/retro/download.py` fetches into `.part` files and renames them when complete, skips finished files, resumes with a range request, keeps the `.part` on cancel, refuses non-HTTPS, writes each fetch to the audit log, and never overwrites a `config.yml` it did not write. Download… runs it as a queued job and sets the config path when it finishes; policy `network_allowed: false` greys it out. `tests/retro/test_retro_download.py` (stub opener, no network) and `tests/app/test_retro_ui.py` (the button, a failure on the row, the policy). Also done: the Routes canvas redraws structures for the zoom they are shown at (steps of √2 × pixel ratio, sharper picture made once zooming pauses), fixing structures drawn at half size on Retina screens; scroll pans, ⌘/Ctrl-scroll and pinch zoom with limits; refits on resize until moved; double-click fits. |
| Later | Route clustering shown as groups in the route table; several algorithms in one pass; a Routes column so routes travel with a row. | — |

## 4. Open questions

- **Install.** Done on this machine on 2026-11-07: the `aizynth` conda environment (Python
  3.12, AiZynthFinder 4.4.1, RDKit 2023.9.6, NumPy 1.26.4, PyArrow 25) and the public USPTO
  models plus the ZINC stock in `~/aizynth-data`, with a `config.yml` beside them. The tests
  still use the fake interpreter, so they need none of it.
- **Where routes live.** They are kept in the Routes panel's saved state, so they reopen with
  the document. A dedicated route column is only worth it if routes should survive copy-paste
  between documents.
