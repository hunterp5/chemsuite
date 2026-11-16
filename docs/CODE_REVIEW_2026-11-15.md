# chemsuite code review, 15 November 2026

Reviewed at `dev` commit `8f1dc77` (about 85,000 lines in `chemsuite/`, 734 tracked files).
Nothing in the code was changed by this review.

## Verdict

The architecture is in good shape and is held there by tests, not by intention: layer
boundaries, the tool API, the ratchets and the process registry all fail a check when broken.
The weak points are at the edges, where data enters and leaves the app, and in a few places
where a single shared resource (the undo history, the worker pool, the molecule cache) is
assumed to behave.

- **2 Critical**, both silent data corruption on import or export, both reproduced.
- **12 Major**: one measured performance bug that also drops undo steps, one export that
  rewrites cell text as Excel formulas, and reliability, integrity and maintainability issues.
- **22 Minor**.

Fix C1, C2, M1 and M2 first. Each is small, self-contained and reproduced below.

**Status, 15 November:** C1, C2, M1 and the Excel half of M2 are fixed in the working tree,
each with a test that fails without the fix. After M1, one-cell edits on the 200,000-row
table take 0.1–0.4 ms and all 60 undo steps are kept. CSV export is unchanged (a policy
decision, see M2). Also done: a damaged document's `IndexError` is now reported as a damaged
file (part of M6), and native crashes write `crash.log` in the log folder (part of M12).
The full suite passes (2,639 passed). Everything else is open.

## How this review was done, and what it did not cover

| Depth | Areas |
|---|---|
| Read in full | `core/` (container, document, download, settings, mols, history, jobs), `core/table/table.py`, `io/` (load, delimited, sdf, formats, xlsx), `app/shell/jobs.py`, `job_host.py`, `main.py`, `app/document/tools.py`, `session.py`, open/save/close in `main_window.py`, `table_model.py`, the image cache in `table_view.py`, `protein/dock/vina.py`, `protein/structure/{components,fetch,store}.py`, `retro/{runner,download}.py`, `tests/test_architecture.py`, CI |
| Skimmed or searched only | `core/table/columns.py`, `analysis/`, `sar/`, `descriptors/`, `conformers/`, `ionisation/`, the tool dialogs, `app/plots`, `app/viewer`, `protein/ui`, `protein/workspace` |
| Not reviewed | `sketch/` and `app/sketcher` (about 6,000 lines), `protein/prepare` (hydrogens, energy, sites), `protein/look`, `molview`, the shaders and the QRhi drawing code, the scientific correctness of any model or descriptor, the test suite's quality |

Each finding says how it was established:

- **Reproduced**: I ran it against this checkout and quote the output.
- **Read**: found by reading the code; not executed. Treat as very likely, not proven.

No full test run was made (nothing changed). No profiler was used apart from the timing in M1,
so the performance ranking outside M1 is from reading the code.

---

## Critical

### C1. Text imports silently replace undecodable characters with `�`  (Reproduced)

`io/formats.py:194` `detect_encoding` decides between UTF-8 and Windows-1252 from the first
64 KB only, and `open_text` (`formats.py:216`) then decodes the whole file with
`errors="replace"`. A CSV whose first non-ASCII character comes after 64 KB is read as UTF-8,
and every such character becomes U+FFFD with no warning.

Reproduced: an 8,001-row Windows-1252 CSV (95 KB) whose last row has the unit `µM` loads with
that cell as `'�M'`. The micro sign is gone, the file opens "successfully", and the
damage is saved into the document.

SD, Mol2 and RD files have no fallback at all: `sdf.py:95` (and the same line in `mol2.py`,
`rdf.py`) always decodes as UTF-8 with `errors="replace"`, so any Latin-1 SD file from an
older vendor tool loses its `µ`, `°`, `Å` and accented names the same way.

Fix: decode strictly. Try UTF-8 over the whole stream (an incremental decoder costs nothing
extra); on the first failure fall back to cp1252 for the file, and use that same helper for
the block formats. Never use `errors="replace"` on user data without counting and reporting
the replacements.

### C2. SD export writes fields that do not read back  (Reproduced)

`io/sdf.py:224` writes each field as `> <name>\n{value}\n\n` with the cell text as it is. The
SD format ends a field at the first blank line and a record at a `$$$$` line, so:

- a cell containing a blank line is truncated at it on re-import;
- a cell containing a line `$$$$` splits the record in two.

Reproduced: a two-row table with notes `"line one\n\nline three"` and `"a\n$$$$\nb"` exports
and reads back as **three** rows with notes `['line one', 'a', '']`.

A column name containing `>` or a newline, and a `Title` cell (only its first line is kept),
have the same class of problem. Multi-line text cells are ordinary (pasted notes, assay
comments), so the first case will happen in normal use.

Fix: on write, replace blank lines inside a value with a single space line and escape or
refuse a `$$$$` line; sanitise field names. Add a round-trip property test (hypothesis is
already a dev dependency) for arbitrary text cells through every export format.

---

## Major

### M1. Undo accounting makes every edit slow and silently drops undo steps  (Reproduced)

`core/history.py:139` recomputes `sum(c.nbytes for c in self._undo)` on every push. A
`PutCells` op counts `old.nbytes`, and for a text column `TextColumn.take` shares the live
column's dictionary, so `TextColumn.nbytes` (`columns.py:440`, `_Dictionary.nbytes` at
`:335`) walks **every string in the whole column** with `sys.getsizeof`, once per command in
the history, on every edit. `ParsedColumn.nbytes` (`columns.py:593`) loops over every blob
the same way.

Two consequences, both measured on a 200,000-row table with a high-cardinality text column
(compound IDs), editing one cell at a time:

```
edit 1: 0.043 s    edit 30: 0.920 s    edit 60: 1.001 s
undo steps kept: 34 of 60; reported history bytes: 524 MB
```

1. A single-cell edit takes about one second once the history has a few dozen entries.
2. The shared dictionary is counted in full for each command, so the 512 MB cap is hit after
   about 34 one-cell edits and older undo steps are discarded. The user sees no message; the
   real memory held is a few kilobytes.

Fix: compute an op's size once when it is created and store it; count only what the op owns
(codes, and the strings it added), never a shared dictionary; keep a running total in
`History` instead of re-summing.

### M2. Excel export turns text starting with `=` into a formula  (Reproduced)

`io/xlsx.py:142` appends cell text to an openpyxl sheet. openpyxl stores any string that
starts with `=` as a formula (`data_type 'f'`); confirmed here for `=1+1` and
`=HYPERLINK(...)`. So a text cell such as `=CH2` or an ID beginning with `=` is exported as
a broken formula, and a value imported from an untrusted SD field (`=HYPERLINK(...)`, DDE
payloads) becomes a live formula in the spreadsheet a colleague opens. CSV export
(`delimited.py`, `write_delimited`) has the same injection exposure when the file is opened
in Excel (`=`, `+`, `-`, `@` prefixes).

Fix: for `.xlsx`, write text through a cell with `data_type = "s"` set explicitly (or
`quotePrefix`). For CSV, decide a policy (prefix `'` on export of text cells that start with
a formula character, or an export option) and document it.

### M3. One crashed or hung worker process breaks every later job  (Read)

- `app/shell/jobs.py:268–280` creates the shared `ProcessPoolExecutor` once and never
  replaces it. If a worker dies (a native crash in RDKit or CDPKit on one molecule, the OOM
  killer), the executor is broken for good: the running job fails with `BrokenProcessPool`,
  and so does every later job that uses the pool, until the app restarts. No code in the
  package mentions `BrokenProcessPool`.
- The job that hit the crash loses all rows already computed. The "keep finished rows" rule
  covers `Cancelled` only.
- `ordered_map` (`core/jobs/parallel.py:162`) cannot stop a chunk that is running. A
  molecule that makes a chunk run for hours (tautomer enumeration, a pathological conformer
  search) keeps its worker after the job is cancelled; enough of them and the pool has no
  free workers. Only Vina has a per-item timeout.
- `ToolContext.in_process` (`app/document/tools.py:264`) uses a one-worker pool. A cancelled
  UMAP or t-SNE keeps running there and the next map waits behind it with no indication.

Fix: catch `BrokenProcessPool` in the runner, discard and rebuild the pool, and surface
"a worker crashed on rows a–b" while keeping the finished rows. For `in_process` and for
cancelled long chunks, replace the worker (terminate and restart) instead of waiting.

### M4. One molecule cache is shared by the window and by background jobs  (Read)

`ParsedColumn.obj` (`core/table/columns.py:543–566`) reads and updates a class-level
`OrderedDict` (`StructureColumn._cache`, 1,024 entries) with no lock. It is called by the
GUI thread while painting (`table_view.py:141`, `table_model.py:184`) and by job threads
exporting a snapshot (`io/sdf.py:215/232`, `mol2.py:314`, `molecules.py:251`). A snapshot
shares the class, so it shares the cache.

- `cache.get(blob)` followed by `cache.move_to_end(blob)` raises `KeyError` if the other
  thread evicted the entry in between. The window is small, but an export of 200,000 rows
  evicts constantly, so scrolling during an export is the case that hits it.
- The same `Chem.Mol` object is handed to two threads. RDKit molecules are not safe for
  concurrent use when one side triggers lazy state (ring info, conformer access).
- An export pushes all 200,000 molecules through the 1,024-entry cache, evicting everything
  the table view had cached.

Fix: jobs should load molecules with `mol_from_blob` directly and never touch the display
cache; make the cache an instance owned by the GUI side, or guard it with a lock.

### M5. Downloads are not verified, and resume can corrupt a file  (Read)

`core/download.py` and `retro/download.py`:

- `retro/download.py:56–76` fetches about 1.4 GB of ONNX models, templates and an HDF5 stock
  that a worker interpreter then loads. There is no pinned SHA-256 or size; integrity rests
  on TLS alone. (The Vina binaries, by contrast, are checksum-pinned in
  `protein/dock/bin/checksums.json`; do the same here.)
- `download.py:97` checks only that the first URL is `https://`. `urllib` follows redirects,
  including to `http://`, and both Zenodo and figshare redirect.
- Resume (`download.py:147–162`) sends `Range: bytes=N-` with no `If-Range`/ETag. If the
  remote file changed, or the `.part` is from another URL, the result is a silent splice of
  two files; `expected = have + length` always matches, so the check passes.
- A `.part` that is already complete (crash between the last write and `os.replace`) gets
  HTTP 416 on every retry and the download can never finish until the user deletes it by
  hand.
- `download_file` has no size cap (only `fetch_bytes` has), and neither sets a read deadline
  beyond the 60 s socket timeout, so a server that trickles one byte a minute holds a job
  forever (cancel still works).

Fix: add `sha256` and `size` to `DataFile`, verify before the rename, delete the `.part` on
mismatch; install a redirect handler that refuses non-HTTPS targets; send `If-Range`; treat
416 as "check the hash and finish".

### M6. Opening an untrusted file trusts its contents too far  (Reproduced in part)

A `.csd` or a compound file from a collaborator is untrusted input. The readers are careful
in places (`np.load(allow_pickle=False)`, validated PDB and component IDs, JSON not pickle
for metadata) but not throughout:

- **Reproduced:** a document whose text-column codes are out of range raises a bare
  `IndexError` from `core/document.py:250`. `read_container` (`container.py:131`) converts
  `KeyError/ValueError/TypeError/zlib.error/struct.error` to "damaged file" but not
  `IndexError`, `OverflowError` or `MemoryError`, so the user gets a traceback instead of
  "this document is damaged".
- Structure blobs are RDKit binary molecules, loaded with `Chem.Mol(blob)`
  (`core/mols.py:56`). That parser is C++ and is not designed for hostile input; a crafted
  blob can crash the process. Blobs are never validated on load, only when first drawn or
  used, so the crash comes later and looks unrelated.
- `zlib.decompress` on every chunk, asset and the oids (`document.py:206/217`,
  `container.py:180`) and `gzip` on imports have no output limit: a small file can expand to
  fill memory.
- `container.py:126` opens the SQLite file read-write ("read only by convention"). Open it
  with `file:…?mode=ro` (and `PRAGMA trusted_schema=OFF`); as written, opening a document on
  a read-only share can fail, and SQLite may create journal files beside it.
- `io/xlsx.py:60`: openpyxl parses XML with the standard library unless `defusedxml` is
  installed; it is not a dependency and is absent from this environment, so entity-expansion
  ("billion laughs") workbooks are not blocked.
- Saved regular expressions in filters and searches are restored and run on the GUI thread
  when a document opens (`app/document/session.py`), so a document can carry a pattern that
  freezes the window (catastrophic backtracking).
- `protein/structure/store.py:unpack` accepts arrays of any shape and dtype; mismatched
  lengths surface later as index errors in drawing code. A missing `_meta` raises `KeyError`
  rather than `StoreError`.

Fix: widen `read_container`'s conversion to `Exception` subclasses that mean "bad data"
(`LookupError`, `ArithmeticError`, `MemoryError`); validate codes and array shapes at load;
use `zlib.decompressobj().decompress(data, max_length)`; add `defusedxml` to the
dependencies; consider validating blobs in the worker pool at open so a bad one is flagged as
an invalid cell.

### M7. A retrosynthesis worker that dies mid-run is reported as success  (Read)

`retro/runner.py:302–307`: if the worker stops after some targets, `search` logs the
exception and returns the replies so far. `tools/retro/dialog.py:finish` then writes those
rows and shows "Searched N structures…". The remaining rows stay blank and nothing tells the
user the run was cut short; only the log has it. Raise a distinct result (for example
`Cancelled`-like `Stopped(partial, reason)`) so the dialog can write the finished rows *and*
say "stopped after N of M: worker exited".

Also in this file, a small race: `start()` replaces `self._replies` while the previous
process's reader thread may still be about to `put(None)` into whatever `self._replies` then
is, which would make the new worker look dead at once. Give each process its own queue
object captured by its reader thread.

### M8. Structure cells are drawn on the GUI thread, through a PNG round trip  (Read)

`app/table/table_view.py:123–187`: on a cache miss the delegate's paint path unpickles the
molecule, computes 2D coordinates, renders with Cairo, encodes a PNG, and decodes it again
with `QImage.fromData` (`:180`). That is several milliseconds per cell, per new cell, on the
thread that must also scroll. With 10–20 structure cells per screen this is the main reason
fast scrolling through a large table will stutter, and a pending (unparsed) cell is also
parsed there. See "Performance" for the proposed change.

### M9. The per-row job driver is copy-pasted across thirteen modules  (Read)

The project rule "per-row tools keep finished rows on cancel" is implemented by repeating the
same loop: make chunks, `ordered_map`, extend the collected lists, report progress, catch
`Cancelled` and re-raise with the partial. It appears in `descriptors/compute.py`,
`ionisation/predict.py`, `analysis/{fingerprints,reactions (×4),reagents,standardize,
structure_search,rgroup_library}.py`, `conformers/etkdg.py` and
`sar/{central_ring,constellation,fragments,mmp}.py`. `_chunked(total, executor)` exists
three times and `_chunks(...)` five times.

This is the largest duplication in the codebase and it is the one a rule depends on. It is
also why M3 has no single place to be fixed. One `core/jobs` helper
(`map_rows(fn, cells, *, executor, cancel, progress, combine)`) would carry the partial-on-
cancel rule, crash handling and progress for every tool.

### M10. `main_window.py` is at its budget and still the hub  (Read)

`app/document/main_window.py` is 1,582 lines against a ratcheted budget of 1,596, with 112
methods and 54 import lines. It owns menus, file open/save/export, drag and drop, background
parsing, clipboard actions for structures and reactions, panel restoration and the close
protocol. The ratchet has stopped it growing, which is the point, but the next feature that
touches the window has nowhere to go. The cheapest cuts, each a clean seam already:
document I/O (`open_path`, `_save_to`, `export_dialog`, `_start_parsing`, about 180 lines)
into a `DocumentFiles` object, and the structure/reaction copy-paste actions (around
`:1185–1320`) next to `structure_copy.py`.

Related: `ToolContext.structure_image` reaches `window._images` (`tools.py:468`), the one
private reach from the tool API into the window.

### M11. The IT policy fails open  (Read)

`core/settings.py`: an unreadable or malformed policy file is logged and ignored
(`:272`), which leaves `network_allowed` at its default of true and unlocks every locked
setting. `machine_policy_path` (`:246`) lets the `CHEMSUITE_POLICY` environment variable
replace the machine path, so any user can point it at an empty file. If the policy is meant
to be a control (the audit log and "downloads are turned off" suggest it is), it should fail
closed when a policy file exists but cannot be read, and the override should be honoured only
when no machine policy exists (or only in development builds).

### M12. No protection against losing work outside a clean save  (Read)

- There is no autosave or recovery file. A native crash (RDKit, the GPU driver, numba) loses
  everything since the last save, and `faulthandler` is not enabled, so such a crash leaves
  no trace in the log either.
- Nothing detects that a document changed on disk since it was opened, and there is no lock
  file. Two windows, two machines or a sync client writing the same `.csd` is last-writer-
  wins with no warning. The project already moves between two machines.

Fix: enable `faulthandler` into the log folder at start-up (one line); write a periodic
recovery snapshot using the existing `snapshot()` + `save_document`; record the file's
mtime and size at open and compare before saving.

---

## Minor

| # | Finding | Where |
|---|---|---|
| m1 | Files written through `mkstemp` + `os.replace` end up mode `0600` on macOS and Linux and lose the previous file's permissions and ACLs: a document saved on a shared drive becomes unreadable to the group. Copy the old mode (or apply the umask) before the rename. Applies to documents, every export and settings. | `core/container.py:66`, `io/formats.py:230`, `core/settings.py:341` |
| m2 | Atomic write is implemented five times with different guarantees: only `io.formats.atomic_write` fsyncs; `SettingsStore.save` leaks the temp file on a non-`OSError`; `CcdCache._write` and the download `.part` use fixed names. One helper in `core`. | as m1, plus `protein/structure/components.py:239`, `core/download.py:146` |
| m3 | `iter_blocks` does `buffer += chunk` and rescans the whole buffer until it finds `$$$$`. A large file with no record terminator (a wrong file renamed `.sdf`, one huge record) is quadratic in time and held entirely in memory. Cap the buffer and fail with a clear message. | `io/sdf.py:70–90`; same shape in `mol2.py:56`, `rdf.py:88` |
| m4 | `csv` field limit is the default 131,072 characters; a cell larger than that (a molfile or sequence in a CSV) aborts the import with `_csv.Error`. Raise the limit. | `io/delimited.py:112` |
| m5 | The CSV header is read with `readline()`, so a quoted header containing a newline is split. | `io/delimited.py:86` |
| m6 | CI tests Python 3.12 only; the package claims 3.11–3.13. `stop_pool` and `ordered_map` use `ProcessPoolExecutor` privates (`_processes`, `_executor_manager_thread`, `_max_workers`) that are free to change between versions. Add 3.11 and 3.13 to the matrix on Linux at least. | `.github/workflows/ci.yml`, `core/jobs/parallel.py:144–145, 193` |
| m7 | Dependencies have lower bounds only and there is no lock or constraints file, so two machines (and CI) can resolve different RDKit, numpy, numba or scikit-learn versions. Results that depend on them (descriptors, conformers, maps) are then not reproducible between the Mac and the PC. | `pyproject.toml` |
| m8 | No static type check in CI although the code is fully annotated and `# type: ignore` is budgeted. | CI |
| m9 | On Unix, `_stop_group` sends `SIGKILL` to the process group after the leader has been reaped; if nothing else was in the group, that id can already belong to an unrelated new process group. Rare. Skip the final kill when `killpg(group, 0)` says the group is empty, or keep the leader unreaped until after it. | `core/jobs/children.py:116` |
| m10 | On Windows a child is assigned to its Job Object after it has started, so anything it spawns in the first moments escapes kill-on-close. Start suspended, assign, resume. | `core/jobs/children.py:60–61` |
| m11 | The Vina fallback copies the binary to a predictable folder in the system temp directory and trusts an existing file of the same size. On a shared Linux machine another user can plant it. The app passes a per-user cache folder, so this only bites library use; remove the fallback or use a `0700` per-user folder. | `protein/dock/vina.py:112` |
| m12 | Any change to cells, rows, columns or styles drops every filter's cached mask, including filters on other columns. Background parsing applies 2,000 rows at a time, so with a filter active a 200k-row import recomputes each mask up to 100 times. Drop only the masks whose `columns` were touched (and the all-column quick search). | `core/table/table.py:810` |
| m13 | Removing or converting a column drops the filters and sort keys that used it, and undo restores the column but not those. | `core/table/table.py:_remove_column`, `_replace_column` |
| m14 | A table listener that raises stops the remaining listeners from hearing `Changed`, leaving panels stale. Isolate listeners (log and continue). | `core/table/table.py:_emit` |
| m15 | Two running instances share `settings.json` (last writer wins: recent files and tool paths can be lost) and the same rotating log files (rotation fails on Windows when another process holds the file). | `core/settings.py`, `core/logbook.py:228` |
| m16 | `.missing` markers in the CCD cache never expire, so a component added to the CCD later is never found. Cached `.cif` files are never revalidated. | `protein/structure/components.py:_fetched` |
| m17 | `gzip.decompress` of a fetched PDB entry has no size limit (up to 256 MB compressed is accepted). | `protein/structure/fetch.py:70` |
| m18 | `_status` builds a one-element numpy array per call and is asked for every painted structure cell; `image()` does the same for its cache key. Index `_blobs[slot]` directly. | `app/table/table_model.py:323`, `table_view.py:144` |
| m19 | `app/shell`, described as the layer every window uses, imports one domain feature directly: `chemsuite.retro` in `settings_dialog.py:71–72, 508, 547, 562` and `main.py:96`. A registry for "external tool tests" and "things to close at exit" would keep the shell free of features, as `Tool` does for menus. | `app/shell` |
| m20 | `core/document.py` serialises columns by reading their private arrays (8 uses of `_values`, `_codes`, `_dict`, `_blobs`, …). Encoding and decoding belong on the column classes. | `core/document.py:156–174, 238–264` |
| m21 | `wait_for` runs a nested event loop for Save. While it spins, the window is modal-blocked by the progress dialog, but timers, job callbacks and other windows still run, so a second save or a close of another window can interleave. Worth a test for "save, then close arrives during the wait". | `app/shell/job_host.py:145–166` |
| m22 | Module-level mutable state makes tests order-dependent and windows implicitly coupled: `_windows`, `_WORKERS`, `children._live`, `notify._tray`, `medchem._sascorer`, the class-level molecule caches, `settings_hub()`. Not wrong for a desktop app, but each is a place a second document can affect the first. | various |

---

## Duplicate logic, redundant utilities, unused code

### Duplicates

| What | Copies | Where | Suggested home |
|---|---|---|---|
| `_mol(cell)`: bytes → `mol_from_blob`, str → `mol_from_text` | 10 identical, plus `mol_of` | `analysis/{fingerprints,reagents,standardize}.py`, `conformers/etkdg.py`, `descriptors/compute.py`, `ionisation/predict.py`, `sar/{central_ring,constellation,fragments,mmp}.py`, `tools/protein/rows.py:81` | `core.mols.mol_from_cell` |
| Chunk → `ordered_map` → collect → `Cancelled(partial)` driver | about 13 | see M9 | one helper in `core/jobs` |
| Atomic write (temp file + replace) | 5 | see m2 | `core` |
| Block reader: `iter_blocks` / `parse_block` / read loop | 3 | `io/sdf.py`, `io/mol2.py`, `io/rdf.py` | a shared record splitter parameterised by terminator (not diffed line by line; same shape) |
| Header handling and the "decide reaction columns from the first rows" closure, and `_batch` | 2 | `io/delimited.py:86–107`, `io/xlsx.py:66–83` | one function taking rows |
| oid → slot lookup with the clip-and-compare idiom | 3 | `Table.slots_of`, `Table.has_oid`, `TableModel.rows_of_oids` | a `Table.find_slots(oids)` returning −1 for missing |
| "parse a number or None" | 4 | `core/table/{criteria,search,split_join}.py`, `protein/prepare/templates.py` | `core/table/numbers.py` |
| `CHUNK` and `TIMEOUT_S` | redefined, unused | `retro/download.py:45–46` (the module delegates to `core.download`) | delete |

### Unused

Found by a name scan (a definition whose name appears nowhere else in `chemsuite`, `tests`
or `scripts`), then checked by search. Dynamic dispatch can hide a use, so confirm before
deleting.

No references anywhere:
`analysis/selectivity.py:100 off_target_candidates`, `app/shell/theme.py:71 fusion_palette`,
`app/sketcher/canvas.py:377 selection_atoms`, `app/table/filters/panel.py:185
set_all_enabled`, `app/table/advanced_search.py:492 remove_row`,
`core/table/columns.py:394 unique_values` and `:435 unused_entries`,
`plotting/data.py:134 is_plottable_axis`, `plotting/orbit.py:106 eye_direction`,
`protein/structure/select.py:126 _per_atom`, `protein/ui/sequence_view.py:369 scroll_to`,
`sketch/edits.py:238 last_edit`, `ionisation/unipka/__init__.py:65 Predictor` (probably
work in progress for the Uni-pKa model, per HANDOFF).

Used only by tests (production code never calls them), a sample:
`Table.clear_filters`, `Table.remove_column` (the app uses `remove_columns`),
`ParsedColumn.invalid_mask`, `core/structure_input.py is_valid` and `mol_weight`,
`protein/dock/box.py from_residue` and `volume`, `protein/chem/interactions.py
ligand_query`, `sketch/layout.py best_layout`, `unipka.macro_pkas`. Either they are API
waiting for a caller or the tests are keeping dead code alive.

The scan also flagged `_Constant`, `_BinOp`, … in `analysis/formula.py` and several Qt
overrides; those are dispatched by name or called by Qt and are in use.

---

## Assumptions that can fail silently

**Inputs**

- Text is UTF-8 if the first 64 KB is (C1). SD field text has no blank lines (C2). Cell text
  does not start with `=` (M2). A CSV cell is under 128 KB (m4). A file ending in `.sdf`
  contains `$$$$` (m3).
- A `.csd` is well-formed apart from the error types `read_container` lists (M6). Structure
  blobs are valid RDKit binaries. Compressed sections are small when expanded.
- Row counts fit in memory several times over: `Table.snapshot()` for every save and export
  copies all numeric arrays and the object-pointer arrays; import builds all batches before
  concatenating. Fine at the 200k-row design point; there is no guard or message at 5M.
- `NumberColumn` is float64: integers above 2^53 (some registry numbers, long barcodes)
  lose digits if a column is inferred as numeric. I did not check `column_from_values`'s
  inference rule for this; worth a test.

**Network**

- The first URL's scheme describes the whole exchange (redirects, M5). A resumed file is the
  same file (M5). The server stays responsive (no overall deadline). A 404 is permanent
  (m16). Downloads are complete if the rename happened (true) and correct if complete (not
  checked).
- `network_allowed()` re-reads the policy file on each download and treats "cannot read" as
  "allowed" (M11).

**Concurrency**

- The table is single-threaded by convention, and the tools respect it: they take inputs on
  the GUI thread and compute on copies. The exceptions are the shared molecule cache (M4)
  and lazily parsing a pending cell inside `obj()`, which writes into a column array from
  whichever thread asked.
- Worker processes do not die and chunks finish (M3).
- A job's `on_done` runs after the document may have changed; `guard` and
  `apply_results`'s deleted-row filter handle that well. `apply_parsed` re-checks the cell
  text. This part is solid.
- One app instance per user (m15). One writer per document (M12).

**Memory**

- Undo memory is estimated, and the estimate is wrong by orders of magnitude for shared
  dictionaries (M1). The image cache is bounded by bytes (192 MB), the molecule cache by
  count (1,024, about right), the log by entries. Decompression is unbounded (M6).
- `TextColumn` dictionaries are append-only; a column edited heavily (find and replace over
  200k unique strings) keeps every old string until `compacted()` runs at save. The live
  dictionary can therefore hold several generations of a column.

---

## Module boundaries and coupling

What is sound, and enforced: Qt-free packages, `core` below the domain packages, the app
never importing a tool, tools going through `ToolContext`, only the bridge knowing both the
table and the workbench, external programs through one registry. I found no violation the
tests miss in those rules.

Where the structure leaks:

1. **`app/shell` → `chemsuite.retro`** (m19). The lowest app layer names one feature.
2. **`core/document.py` ↔ column internals** (m20). The file format is coupled to private
   attributes of four classes; adding a column kind means editing the serialiser.
3. **`core/jobs/parallel.py` ↔ `concurrent.futures` privates** (m6).
4. **`ToolContext` → `window._images`** and `StructureImages` living in
   `app/table/table_view.py` while `app/shell/main.py`, the viewer and `ToolContext` all use
   it. It is a shared service; move it to `app/widgets` or `app/shell` and hand it to the
   context.
5. **`protein/ui/panel.py` and `protein/workspace/session.py` still read the `Table`**: a
   known ratchet in `test_import_boundaries.py`, noted here so it stays visible.
6. **Tools still reading column internals**: 40 budgeted uses across nine files
   (`COLUMN_INTERNALS_IN_TOOLS`), the largest being `tools/retro/dialog.py` (9) and
   `tools/utilities/calculator.py` (7).
7. **The main window** (M10) is the most coupled component: it is the only class that knows
   the runner, the workspace, the view, the panels, the session and the file formats at
   once.
8. **The class-level molecule cache** couples every document and every thread through one
   dictionary (M4).

The `getattr`/`hasattr` budget (54) is full, as CLAUDE.md says. `core/download.py` alone
uses four on the response object; a tiny `Response` protocol would free them.

---

## Performance: where the gains are

Ranked by expected effect on what a user feels. Only the first is measured.

1. **Undo accounting (M1).** Measured cost today: one-cell edits climb from 43 ms to about
   1 s at 200k rows as the history fills. With sizes stored per op the cost should stay near
   the first edit's (not measured, since nothing was changed). Small change, large effect,
   and it fixes the lost undo steps.
2. **Structure drawing off the GUI thread (M8).** Paint a placeholder on a cache miss, render
   in the quick lane (or the process pool for big screens), and repaint the cell when the
   image lands; drop results for cells that scrolled away. Skip the PNG encode/decode by
   avoiding the Cairo-to-PNG path (options to evaluate: `MolDraw2DSVG` with Qt's SVG
   renderer, or drawing with `QPainter` from RDKit's coordinates as the sketcher's renderer
   does; I did not prototype either). Pre-render one screen ahead in the
   scroll direction.
3. **Filter masks (m12).** Per-column invalidation makes editing and background parsing with
   filters active proportional to the filters touched, not to all of them.
4. **Row insert and delete reset the model.** `TableModel._resync` calls `beginResetModel`
   whenever the row count changes, which throws away the view's layout and repaints
   everything. Appending rows (the common case: enumeration results, pasted rows) could use
   `beginInsertRows` when the new oids sort after the old ones and no sort is active.
5. **Exports bypass the display cache (M4)** and stop evicting the table's molecules.
6. **A dedicated, restartable worker for GIL-bound work (M3)** so a cancelled map does not
   delay the next one.
7. **Save**: single-threaded zlib level 1 per 64k-row chunk, column by column. Compressing
   chunks in the thread pool (zlib releases the GIL) would scale with cores; the SQLite
   writes stay on one thread. Worth doing only if save time is a complaint.
8. **Open**: every column of a `.csd` is decoded up front. The format already stores columns
   separately "so a later version can read single columns lazily"; that later version is the
   biggest available win for opening very wide documents.

Not measured and worth a profile before any work: start-up import time (RDKit, scikit-learn
and pyqtgraph are heavy; several are already imported lazily), and the 3D views.

---

## Blind-spot audit

Things this review would otherwise have let pass without comment.

**Areas I did not read and you should not assume are clean**

- The sketcher (`sketch/`, `app/sketcher`, about 6,000 lines including the 1,341-line
  canvas) and all of `protein/prepare`. Both parse or build chemistry by hand and are where
  subtle correctness bugs would live.
- The GPU code. I did not check buffer lifetimes, device-loss handling or what happens when
  a QRhi resource fails to create mid-session (start-up is checked by `graphics.py`).
- Scientific validity: descriptor definitions, pKa models, interaction cut-offs, MMP and
  Free–Wilson maths. A code review does not establish that a number is right.
- The tests themselves. About 2,050 tests exist and CI runs them on three systems, but I
  did not assess what they assert. The four reproduced bugs above each sit in a module with
  a test file, which says the tests cover behaviour more than edge cases. Hypothesis is
  installed; round-trip properties for every import/export format and for `.csd` would have
  caught C1, C2 and M2.

**Risks outside the code**

- **Distribution.** There is no frozen build, installer, code signing or notarisation yet
  (`freeze_support()` is the only trace). Shipping to chemists on managed Windows and macOS
  machines is a project of its own: unsigned binaries, the bundled `vina` executables, the
  spawn-based process pool under a frozen app, and Qt plugin discovery all tend to fail
  there first.
- **Licences.** GPL-3.0 with RDKit (BSD), PySide6 (LGPL), CDPKit (LGPL-2.1), Vina
  (Apache-2.0), pacmap/faiss, umap/numba, openpyxl, gemmi (MPL-2.0), meeko (LGPL-2.1),
  propka. I did not audit these or the bundled data; CLAUDE.md records licence decisions
  for the pKa data, and the `security-auditor` agent in this repo is set up for exactly this
  check.
- **Dependency drift** (m7) and no vulnerability scan (`pip-audit`) in CI. GitHub Actions
  are pinned by tag, not commit.
- **Privacy of logs.** The log and the optional audit log record full file paths and every
  URL fetched (PDB IDs and ligand codes reveal what someone is working on). That is the
  purpose of an audit log, but the ordinary log does it too, and logs get attached to bug
  reports.
- **Accessibility and localisation.** Structure cells are images with no accessible text;
  number parsing and formatting appear to assume `.` as the decimal separator. I did not
  check what a German-locale CSV with `1,5` becomes on import; test it.
- **Scale ceiling.** Everything is in memory by design. There is no message when a file is
  too big; the app will simply get slow and then fail to allocate. DATABASE_PLAN.md exists,
  so this is known, but the failure mode today is unannounced.

**Where my own findings are weakest**

- M3, M4, M5, M7 and M8 are from reading. M4's race needs two threads in a narrow window; I
  am confident the code permits it, not that anyone has hit it.
- The "unused" list is a name scan.
- Severity is judged for a desktop tool used by a small group with files from colleagues.
  If documents will come from outside the organisation, M6 moves up to Critical.

---

## Suggested order

1. C1, C2, M2: import/export correctness. Add round-trip property tests with them.
2. M1: store op sizes; running total.
3. M3 + M9 together: one row-mapping helper that also rebuilds a broken pool.
4. M5: hashes and redirect check for downloads. M6: error conversion, `mode=ro`,
   `defusedxml`, decompression limits.
5. M12 (first step): `faulthandler`, and an mtime check before save.
6. M4, M8: molecule cache ownership, then asynchronous cell drawing.
7. The duplicates table and m1/m2 as routine clean-up; M10 before the next feature that
   needs the main window.

Items 3 and 6 are structural (shared infrastructure, the renderer path) and by the project's
own rule should start in plan mode and end with a full test run.
