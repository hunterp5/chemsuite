# This file is part of chemsuite.
# Copyright (C) 2026 Hunter Picard
#
# chemsuite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# chemsuite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with chemsuite.  If not, see <https://www.gnu.org/licenses/>.
"""Minimisation in the workbench: a prepared protein relaxed in OpenMM's own interpreter
(:mod:`chemsuite.simulate`), kept as a **new protein** beside the one it was made from.

- Protein ▸ Prepare… with Minimise ticked leaves the protein as it was and adds two:
  ``X-prep`` (prepared) and ``X-min`` (prepared, then minimised), in one undoable step.
- Protein ▸ Minimise… adds ``X-min`` for a protein already prepared.

Each copy names the protein it was made from (``derived_from`` in its meta), and the
workbench shows the copy in its place (:meth:`WorkbenchSession._swap_sources`). Each copy
keeps its own :class:`~chemsuite.protein.workspace.prepared.PreparedRecord` and report.

A minimisation always starts from the *prepared* model. Its restraints and its RMSD are
measured from the coordinates passed in, so minimising a minimised model again would drift
away from the deposit. That is why a minimised protein cannot be prepared or minimised again.
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Mapping
from dataclasses import dataclass, replace

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.shell.settings_hub import settings_hub
from chemsuite.core.assets import Asset, PutAsset
from chemsuite.core.jobs import CancelToken
from chemsuite.protein.structure import store
from chemsuite.protein.structure.components import CcdCache
from chemsuite.protein.structure.model import WATER, Structure
from chemsuite.protein.workspace.prepared import Minimisation, PreparedRecord, atoms_digest
from chemsuite.protein.workspace.reports import PrepareReport, Row
from chemsuite.protein.workspace.session import DERIVED_FROM
from chemsuite.simulate import build, config, runner
from chemsuite.simulate.minimise import Done, minimise
from chemsuite.simulate.result import Contact
from chemsuite.simulate.spec import STAGES, MinimiseSettings, SpecError

from .tool import KIND, ccd, protein_workspace

#: an asset's meta key saying what made it, and its values
STAGE = "stage"
PREPARED, MINIMISED = "prepared", "minimised"
SUFFIX = {PREPARED: "-prep", MINIMISED: "-min"}
#: what stops a minimisation without it being a bug (said in the report, not a crash box)
FAILURES = (build.BuildError, runner.SimulateError, SpecError)

_last = MinimiseSettings()  # remembered for the session


def last_settings() -> MinimiseSettings:
    return _last


def remember(settings: MinimiseSettings) -> None:
    global _last
    _last = settings


# -- the worker ----------------------------------------------------------------------------------
def ready() -> str:
    """Why OpenMM cannot run here, or ``""``. Starts no process."""
    return config.ready_problem(settings_hub().settings)


def simulator() -> runner.Simulator:
    """The warm OpenMM worker (check :func:`ready` first)."""
    where = config.interpreter(settings_hub().settings)
    if where is None:
        raise runner.SimulateError(ready())
    return runner.worker_for(where)


def platforms() -> tuple[str, ...]:
    """The platforms the worker has said work, this session (none until one has run)."""
    where = config.interpreter(settings_hub().settings)
    return runner.platforms_seen(where) if where is not None else ()


# -- the proteins --------------------------------------------------------------------------------
def stage_of(ctx: ToolContext, protein: str) -> str:
    """What made *protein*: PREPARED, MINIMISED, or "" for one opened or fetched."""
    asset = ctx.table.assets.get(protein)
    return str(asset.meta.get(STAGE, "")) if asset is not None else ""


def root_of(ctx: ToolContext, protein: str) -> str:
    """The protein *protein* was first made from (itself, for one opened or fetched): the
    name its copies are called after, so a copy of ``X-prep`` is ``X-min``, not
    ``X-prep-min``."""
    seen = {protein}
    while True:
        asset = ctx.table.assets.get(protein)
        source = str(asset.meta.get(DERIVED_FROM, "")) if asset is not None else ""
        if not source or source in seen:
            break
        seen.add(source)
        protein = source
    for suffix in SUFFIX.values():
        if protein.endswith(suffix) and stage_of(ctx, protein):
            return protein[: -len(suffix)]
    return protein


def copy_id(ctx: ToolContext, protein: str, stage: str, taken: Collection[str] = ()) -> str:
    """A new id for *protein*'s copy at *stage* (``X-min``, ``X-min-2``…); *taken* holds ids
    claimed in the same step but not yet added."""
    base = root_of(ctx, protein) + SUFFIX[stage]
    used = {*ctx.table.assets, *taken}
    found, n = base, 2
    while found in used:
        found, n = f"{base}-{n}", n + 1
    return found


def copy_op(new: str, source: str, s: Structure, stage: str, meta: Mapping[str, str]) -> PutAsset:
    """Adding *s* as protein *new*, made from *source* at *stage*; *meta* is the original's
    (its file, its ligands' ⚠)."""
    meta = {**meta, DERIVED_FROM: source, STAGE: stage}
    return PutAsset(new, None, Asset(new, KIND, new, s, meta))


def minimised_already(ctx: ToolContext, protein: str) -> str:
    """The ⚠ for a protein that is minimised already ("" when it is not)."""
    record = protein_workspace(ctx).prepared(protein)
    if stage_of(ctx, protein) != MINIMISED and (record is None or record.minimised is None):
        return ""
    asset = ctx.table.assets.get(protein)
    source = str(asset.meta.get(DERIVED_FROM, "")) if asset is not None else ""
    instead = f": minimise {source} instead" if source in ctx.table.assets else ""
    return f"{protein} is minimised already{instead}"


def has_waters(s: Structure) -> bool:
    return bool((s.atom_kind == WATER).any())


def problem(ctx: ToolContext, protein: str, settings: MinimiseSettings) -> str:
    """Why Protein ▸ Minimise… cannot run on *protein*, or ``""``."""
    if why := ready():
        return why
    if why := minimised_already(ctx, protein):
        return why
    asset = ctx.table.assets.get(protein)
    if asset is None or asset.kind != KIND:
        return "There is no such protein"
    record = protein_workspace(ctx).prepared(protein)
    if record is None:
        return f"Prepare {protein} first"
    s = store.as_structure(asset.payload)
    if not record.matches(s):
        return f"{protein} was edited after it was prepared: prepare it again"
    return settings.problem(has_waters(s))


# -- the work ------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Relaxed:
    """A minimisation done: the structure it gave (one model), and what it was."""

    structure: Structure
    done: Done
    settings: MinimiseSettings


def relax(
    worker: runner.Simulator,
    s: Structure,
    record: PreparedRecord,
    settings: MinimiseSettings,
    ccd: CcdCache | None = None,
    model: int = 0,
    report: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> Relaxed:
    """*s*'s *model* minimised as *record* says it was prepared (on a worker thread). Raises
    one of :data:`FAILURES` with what stopped it."""
    system = build.system_spec(s, record.states_on(s), ccd=ccd, cancel=cancel)
    names = [stage.name for stage in STAGES]

    def progress(stage: str, _iteration: int, _energy: float) -> None:
        if report is not None and stage in names:
            report(names.index(stage), len(names))

    done = minimise(worker, s, system, settings.forcefield(), settings.protocol(),
                    built=record.built_in(s), model=model, platform=settings.platform,
                    progress=progress, cancel=cancel)  # fmt: skip
    out = s.with_coords(done.result.coordinates(s, model))
    return Relaxed(out, done, settings)


def record_of(new: str, source: str, record: PreparedRecord, relaxed: Relaxed) -> PreparedRecord:
    """The minimised protein's record: the preparation's (the same atoms), and how it was
    minimised. The coordinates are the protein's own, so the result keeps the numbers only."""
    result = relaxed.done.result.to_json()
    for heavy in ("uid", "xyz"):
        result.pop(heavy, None)
    minimised = Minimisation(source, relaxed.settings.forcefield().to_json(),
                             relaxed.done.system.to_json(), result)  # fmt: skip
    return replace(record, protein=new, atoms=atoms_digest(relaxed.structure),
                   minimised=minimised)  # fmt: skip


def _contact(c: Contact | None) -> str:
    return "none" if c is None else f"{c.a} – {c.b}, {c.distance:.2f} Å"


def report_of(new: str, relaxed: Relaxed, ph: float) -> PrepareReport:
    """What the minimisation did, for the report panel: the energies, how far the deposited
    atoms moved, contacts, the force left, the residues that moved most (each a line that
    zooms to it), and what it was run with."""
    s, result, settings = relaxed.structure, relaxed.done.result, relaxed.settings
    where = {s.residue_label(r): r for r in range(s.n_residues)}
    rows = [Row("Minimised", f"Energy {result.energy_before:,.0f} → "
                f"{result.energy_after:,.0f} kJ/mol")]  # fmt: skip
    rows += [Row("Minimised", text, True) for text in result.warnings(settings.rmsd_limit)]
    moved = f"Heavy atoms moved {result.rmsd:.2f} Å (RMSD)"
    if result.rmsd_pocket is not None:
        moved += f", pocket {result.rmsd_pocket:.2f} Å"
    if result.rmsd_ligand is not None:
        moved += f", ligand {result.rmsd_ligand:.2f} Å"
    rows.append(Row("Minimised", moved))
    rows.append(Row("Closest contact", f"Before: {_contact(result.contact_before)}; after: "
                    f"{_contact(result.contact_after)}"))  # fmt: skip
    if result.max_force_atom:
        residue = result.max_force_atom.rsplit("/", 1)[0]
        rows.append(Row(residue, f"Largest force left: {result.max_force:,.0f} kJ/mol/nm on "
                        f"{result.max_force_atom.rsplit('/', 1)[-1]}",
                        residue=where.get(residue, -1)))  # fmt: skip
    rows += [Row(label, f"Moved {distance:.2f} Å, the most of any residue"
                 if i == 0 else f"Moved {distance:.2f} Å", residue=where.get(label, -1))
             for i, (label, distance) in enumerate(result.moved_most)]  # fmt: skip
    ff = settings.forcefield()
    charges = sorted({model for _label, model, _q in result.charges if model})
    used = f"{ff.protein} + {ff.water} + {ff.small_molecule}, {ff.solvent}"
    if charges:
        used += f"; ligand charges {', '.join(charges)}"
    rows.append(Row("Force field", used))
    versions = dict(result.versions)
    run = f"{result.platform}, {sum(seconds for _n, seconds in result.seconds):.1f} s"
    if versions.get("openmm"):
        run += f", OpenMM {versions['openmm']}"
    rows.append(Row("Run", run))
    return PrepareReport(new, tuple(rows), ph=ph, net_charge=int(s.charge.sum()))


def failed_row(error: BaseException) -> Row:
    """The ⚠ line a preparation's report gets when its minimisation could not be run."""
    text = str(error).splitlines()[0] if str(error) else type(error).__name__
    return Row("Minimised", f"Not minimised: {text}", True)


# -- Protein ▸ Minimise… -------------------------------------------------------------------------
def minimise_protein(ctx: ToolContext, protein: str, settings: MinimiseSettings):
    """Minimise *protein* (prepared) in the background and add the result as a new protein,
    with its report. Returns the job, or None when it cannot run (:func:`problem`)."""
    from .prepare_dialog import show_report

    if problem(ctx, protein, settings):
        return None
    workspace = protein_workspace(ctx)
    record = workspace.prepared(protein)
    s = store.as_structure(ctx.table.assets[protein].payload)
    shown = workspace.session().display(protein)
    model = shown.model if shown is not None else 0
    if s.n_models > 1:
        s = s.with_coords(s.coords[model])
        model = 0
    worker = simulator()
    cache = ccd(ctx)
    remember(settings)

    def work(job) -> Relaxed:
        return relax(worker, s, record, settings, cache, model, job.report, job.cancel_token)

    def done(relaxed: Relaxed) -> None:
        if protein not in ctx.table.assets:
            ctx.show_message(f"{protein} is no longer in this document")
            return
        new = copy_id(ctx, protein, MINIMISED)
        meta = ctx.table.assets[protein].meta
        ctx.edit(f"Minimise {protein}", [copy_op(new, protein, relaxed.structure, MINIMISED, meta)])
        workspace.set_prepared(record_of(new, protein, record, relaxed))
        report = report_of(new, relaxed, record.ph)
        workspace.set_report(report)
        show_report(ctx, report)
        result = relaxed.done.result
        warned = f", {report.warnings} to check" if report.warnings else ""
        ctx.show_message(f"Minimised {protein} as {new}: heavy atoms moved {result.rmsd:.2f} Å"
                         f"{warned}")  # fmt: skip

    return ctx.run(f"Minimising {protein}…", work, done)
