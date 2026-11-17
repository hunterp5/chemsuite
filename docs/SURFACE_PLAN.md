# Protein surface options (Macromolecule Viewer; plan approved 2026-11-17)

## Context

Asked for 2026-11-15 (HANDOFF, "Next (protein)"). Today each protein has one surface: a
Gaussian "blob" (`macromolecule/look/surface.py`, VMD QuickSurf-like) over the atoms the user
showed it on, always colored by hydrophobicity (`Display.surface_scheme` is never set from the
UI), with one transparency for every protein (Display ▸ Surface Transparency). Focus draws a
pocket skin with the same fixed look. Chemists expect what PyMOL and ChimeraX give: the
molecular (SES) surface by default, other kinds, color by any scheme or electrostatics,
opacity per surface, mesh and dot styles, carving around a ligand, smoothing.

Decided with the user (2026-11-17):
- **Default kind: SES** (probe 1.4 Å) for new surfaces; documents saved before load as Blob.
- **Coulombic charges:** Amber ff14SB-style partial charges per atom, with the protonation
  state a prepared protein has (read from which hydrogens are present, as ChimeraX does);
  united-atom defaults (hydrogen charges folded onto their heavy atom) for an unprepared
  protein; Gasteiger for ligands; ions their formal charge; waters left out.
  Note: Prepare itself stores only formal charges and hydrogens, not partial charges. So
  "the charges Prepare assigned" means the ff14SB charges of the states Prepare chose.

This file is the approved plan, linked from ROADMAP §7.0c and HANDOFF.

## Design

### Qt-free geometry (`chemsuite/macromolecule/look/`)

- **`surface.py`**: `build(structure, atoms, model, options)` where `options` is a frozen
  `SurfaceOptions(kind, probe, quality, smoothing)`.
  - Kinds `ses` (default), `sas`, `vdw`, `blob` (today's, unchanged).
  - vdW/SAS: a signed distance field `min_i(|x − cᵢ| − (rᵢ + p))` on the grid, computed per
    radius class with `scipy.ndimage.distance_transform_edt(..., return_indices=True)` on the
    voxelized centers and then the exact distance to the nearest atom's real position (no
    Python loop over atoms; few radius classes).
  - SES: the SAS field's inside, its EDT (distance to the solvent region), contoured at the
    probe radius (EDTSurf's method). Exact enough at ≤ 0.5 Å, and smoothed (below).
  - The same `_surface_nets` contouring, padding (closed meshes), gradient normals and
    nearest-atom owners as today. Quality picks the spacing (Coarse 0.8, Normal 0.5, Fine 0.35
    Å), still widened past `MAX_VOXELS`.
  - **Smoothing** (Off / Low / High): Taubin λ|μ smoothing of the vertices over the mesh's
    edges (no shrinkage), normals recomputed from the faces; Blob ignores it.
- **`mesh.py`**: `edges(mesh)` (unique edges, for the Mesh style); `dots(mesh, spacing)` (small
  octahedra on a vertex subsample, for the Dots style, owners kept); `carve(mesh, points,
  distance)` keeps triangles whose vertices all lie within *distance* of the carve points
  (cKDTree), so a carved surface keeps the whole protein's shape. The carve points are the
  positions of the atoms selected when carving, taken from every visible protein (a ligand
  or pose can belong to another one), and are stored as points. Edits never have to remap
  them.
- **`electrostatics.py`** (new): `charges(structure, ligand_charges)` per atom (the table
  below, variants chosen by the hydrogens present, united atom without them);
  `potential(points, centres, charges)` = Σ 332·q/(ε r) with ε = 4r, a cKDTree cutoff of
  15 Å; colored red (−10) – white – blue (+10) kcal/(mol·e), as ChimeraX's default.
- **Charge table**: `look/data/amber_charges.json` (residue → atom → charge, including the
  His/Asp/Glu/Lys/Cys/Tyr/terminal variants), written once by `scripts/build_amber_charges.py`
  from a permissively licensed source. **Step 0 checks the licence** (security-auditor):
  PDB2PQR's `AMBER.DAT` (BSD-3) first, OpenMM's ff14SB XML (MIT) second; the source and paper
  (Maier et al. 2015) are cited in the file. AMBER's own leap files do not ship.

### What a display holds (`look/display.py`)

- New per-protein fields, saved in `Display.state()` (old documents: `kind` missing → Blob):
  `surface_options` (kind, probe, quality, smoothing), `surface_style` (solid | mesh | dots),
  `surface_scheme` (any of `colours.SCHEMES`, plus `coulombic` and `flat`), `surface_rgb`
  (for flat), `transparency` (exists), `carve` (points, (k, 3)) and `carve_distance`.
- **Edits keep them:** the display rebuilt after a structure edit (today it copies `visible`,
  `scheme`, `surface_scheme` and remaps `surface`) also carries over the options, style, flat
  color, transparency and carve. This gets a test (an Inspect apply or residue edit must not
  reset a surface).
- The `_rgba` cache key adds the flat color.
- The surface part's build gets the options. Carving and the style are applied after the build,
  so a carve or style change rebuilds nothing heavy.
- **`potential`, a new deferred part:**
  - It is built in the background only when the scheme is Coulombic.
  - Its version is tied to the surface's version and the charges. It is handed out only once
    that surface version has been accepted.
  - It is computed on the uncarved mesh, so changing the carve costs nothing.
- The pocket skin gets its own transparency and style (session state, as `pocket_surface` is)
  and uses the same kind as its protein's surface.
- `scene.compose` and **`relook`**: per-vertex colors for Coulombic. `meshes.look` takes
  per-vertex colors as an alternative to per-owner ones, and `relook` keeps them on a selection
  change. The layers carry their style.

### Drawing (`macromolecule/ui/view.py`)

- `set_mesh(..., style=)`: Solid as today. **Dots** is the octahedra mesh through the
  existing triangle pipelines. **Mesh** draws the mesh's edges with a `Lines` topology built
  from the existing `pro_mesh` shaders (no new shader, no `.qsb` rebuild).
  - Mesh is more than a new pipeline. `meshes.chunks`, `_upload_chunks` and `_draw_mesh`
    assume triangles (`mesh.triangles`, `.size` counts), so a layer's index array and topology
    are generalized, with a pipeline for each mode and topology.
  - Check first that `RhiView.pipeline()` accepts a Lines topology.
  - Both styles work opaque, transparent, and in chunks without 32-bit indices.

### Menus (`macromolecule/ui/menus.py`, `menubar.py`)

Display ▸ **Surface** replaces Display ▸ Surface Transparency. It acts on the proteins with a
selection, else every visible one (as Color does), and checks show the first one's value:
- **Kind:** Molecular (SES) · Solvent Accessible (SAS) · Van der Waals · Gaussian
- **Style:** Solid · Mesh · Dots
- **Color:** By Hydrophobicity, By Element, By Chain, … (every scheme), **By Electrostatics**,
  **Single Color…**
- **Transparency:** Opaque · 20% … 80%
- **Carve Around Selection** (3 · 5 · 8 · 12 Å) and **Show Whole Surface**
- **Probe Radius** (1.0–2.0 Å; SES and SAS only), **Quality** (Coarse · Normal · Fine),
  **Smoothing** (Off · Low · High)

View ▸ Ligand ▸ **Pocket Surface** becomes a submenu: On, Transparency, Style.
Tooltips explain (no helper text). New text is American English ("Color"); the Display ▸
Colour menu is renamed Color on the way.

## Phases (each a commit on `dev`; targeted tests + architecture/boundaries per step)

0. **Licence check** of the charge source (security-auditor): PDB2PQR's `AMBER.DAT` (BSD-3)
   first, then OpenMM's ff14SB XML (MIT).
   - The generator parses the file directly and never imports OpenMM or pdb2pqr. It is data
     only, so "pdb2pqr is dropped" still holds.
   - **Checkpoint:** report the source, its licence and the citation to the user, and commit
     the JSON only after their OK. The decision then goes in CLAUDE.md's Decided section.
1. **Geometry:** kinds, quality, smoothing, carve, edges, dots in `look/`. `surface.build`
   always names its kind. Tests in `tests/macromolecule/look/`:
   - every kind's mesh is closed (`open_edges == 0`);
   - one atom's vdW and SAS areas are near 4πr², with the tolerance set by measuring what
     surface nets gives at 0.5 Å first (no code tuned to a guessed number);
   - SES encloses less than SAS and is smooth over the cleft between two spheres;
   - smoothing keeps the mesh closed and moves vertices less than 0.3 Å;
   - carve keeps only the triangles near the points;
   - Blob is unchanged.

   A benchmark (`-m bench`) measures SES **time and peak memory** on a test protein (target:
   under 2 s at Normal for about 5,000 heavy atoms). `distance_transform_edt`'s index arrays
   near `MAX_VOXELS`, once per radius class, are watched.
   - With SES the new default, grep the tests and benchmarks that build surfaces through
     `Display` (`show_surface`, `surface_mesh`, `test_bench_protein.py`). Keep them on the kind
     they meant; `test_a_selection_change_is_drawn_within_a_frame` is already timing-sensitive.
2. **Display, session, saving, menus:** the fields above, `Display.state` round-trip and old
   state loading as Blob, session setters on the selected-or-all proteins, the Surface menu and
   the pocket submenu. Tests in `tests/macromolecule/look/test_display*.py`,
   `tests/macromolecule/ui/`.
3. **Styles in the view:** Mesh (lines pipeline) and Dots. GUI tests check the layers and their
   buffers offscreen; screenshots in the app.
4. **Electrostatics:** `electrostatics.py`, the `potential` part, By Electrostatics. The
   15 Å neighbor sum over about 10⁵ vertices runs in chunks and goes in the benchmark. Tests:
   - a lone Lys NZ makes the nearby surface blue, and an Asp's red;
   - His charges differ between prepared and unprepared;
   - a ligand gets Gasteiger charges;
   - the potential follows the 4r dielectric;
   - a selection change keeps the Coulombic colors.
5. **Docs:** `docs/USER_GUIDE.md` (the viewer's surface section), ROADMAP §7.0c, HANDOFF.

## Verification

- Per phase: `pytest tests/macromolecule/look tests/macromolecule/ui
  tests/tools/macromolecule tests/test_architecture.py tests/test_import_boundaries.py`,
  `ruff check` / `ruff format` on touched paths.
- The full suite once at the end (the renderer and document save/load are touched), and say so.
- In the app (run skill), on a sample protein (e.g. 1ATP from `tests/data/protein`): SES vs
  Blob side by side, each style, Coulombic before and after Prepare, carve around the ligand,
  pocket transparency; save, reopen, and an older document opening as Blob. Before/after
  screenshots.
