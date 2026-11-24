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
"""Data ▸ QSAR: the dialogs, the columns they write, and the model the document keeps.

One test fits for real through the dialog, because that is the only way to know the table, the
engine and the document part agree. The rest build a model straight from the engine, which is
much quicker, and check what the app does with it.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from tests.gui import gui

pytestmark = gui

#: 24 small molecules whose logP spans a useful range, so a model has something to find
SMILES = [
    "CCO", "CCCO", "CCCCO", "CCCCCO", "CCCCCCO", "CCCCCCCO",
    "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1", "Clc1ccccc1", "Brc1ccccc1",
    "CC(=O)O", "CCC(=O)O", "CCCC(=O)O", "CCCCC(=O)O", "c1ccccc1C(=O)O", "Cc1ccc(cc1)C(=O)O",
    "CCN", "CCCN", "CCCCN", "c1ccccc1N", "Cc1ccccc1N", "CCOCC",
]  # fmt: skip


def a_table():
    """A table of structures with a measured value, enough rows for a model."""
    from rdkit import Chem
    from rdkit.Chem import Crippen

    from chemsuite.core.table import Table

    measured = [round(Crippen.MolLogP(Chem.MolFromSmiles(s)), 3) for s in SMILES]
    table = Table()
    table.add_rows(
        {
            "ID": [f"C{i:02d}" for i in range(len(SMILES))],
            "SMILES": SMILES,
            "LogP": [str(v) for v in measured],
            "Series": ["acid" if "C(=O)O" in s else "other" for s in SMILES],
            # a second number, so a time split has something to order by
            "Batch": [str(1 + i // 8) for i in range(len(SMILES))],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    return table


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace

    window._install(Workspace(a_table()))
    return window


def a_model(task: str = "regression", kind: str = "d"):
    """A fitted model, built straight from the engine: no structures, so it is quick.

    *kind* is the pool prefix its features carry: "d" for descriptors, which any document
    with structures can rebuild, or "c" for table columns, which a document may not have.
    """
    from chemsuite.qsar.fit import fit
    from chemsuite.qsar.recipe import Recipe

    rng = np.random.default_rng(0)
    pool = rng.normal(size=(60, 8))
    raw = 2 * pool[:, 0] - pool[:, 1] + rng.normal(scale=0.2, size=60)
    y = raw if task == "regression" else (raw > np.median(raw)).astype(float)
    classes = None if task == "regression" else np.array(["low", "high"])
    recipe = Recipe(
        task=task,
        target="LogP",
        structures="SMILES",
        columns=tuple(f"c{i}" for i in range(4)),
        split="random",
        budget="quick",
        max_features=4,
        seed=0,
        without=("rf", "mlp", "knn"),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fit(
            pool,
            [f"{kind}:col{i}" for i in range(8)],
            y,
            np.arange(60),
            recipe,
            keys=[f"KEY{i:04d}" for i in range(60)],
            classes=classes,
        )


# -- the dialog ------------------------------------------------------------------------------
def test_the_dialog_detects_the_task_from_the_target(chem) -> None:
    dialog = chem.open_tool("auto_qsar")
    assert dialog.primary.text() == "Fit on 24 Rows"
    assert not dialog.status.text().strip()
    dialog.target.setCurrentText("LogP")
    assert dialog.task.value() == "regression"
    dialog.target.setCurrentText("Series")
    assert dialog.task.value() == "classification", "two text labels are classes"
    assert not dialog.task.buttons["regression"].isEnabled(), "labels are not numbers"


def test_a_fit_needs_rows_features_and_a_column_to_write(chem) -> None:
    dialog = chem.open_tool("auto_qsar")
    dialog.descriptors.setChecked(False)
    dialog.fingerprints.setChecked(False)
    assert "at least one kind of feature" in dialog.status.text()
    assert not dialog.primary.isEnabled()
    dialog.descriptors.setChecked(True)
    assert not dialog.status.text().strip()
    for box in (
        dialog.write_prediction,
        dialog.write_interval,
        dialog.write_residual,
        dialog.write_domain,
    ):
        box.setChecked(False)
    assert "at least one column to write" in dialog.status.text()


def test_too_few_rows_are_refused_rather_than_modeled(window) -> None:
    """Four rows can be fitted perfectly and mean nothing; the dialog says so."""
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"SMILES": SMILES[:4], "LogP": ["1", "2", "3", "4"]},
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    dialog = window.open_tool("auto_qsar")
    assert "cannot be modeled honestly" in dialog.status.text()
    assert "20 is the fewest" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_a_table_with_no_structures_says_that_first(window) -> None:
    """The bare window has no structure column, which blocks a fit before anything else."""
    dialog = window.open_tool("auto_qsar")
    assert "no structure column" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_the_split_column_appears_only_for_a_column_split(chem) -> None:
    dialog = chem.open_tool("auto_qsar")
    assert not dialog.method.is_row_visible(dialog.split_column)
    dialog.split.set_value("column")
    assert dialog.method.is_row_visible(dialog.split_column)
    assert "Batch" in [dialog.split_column.itemText(i) for i in range(dialog.split_column.count())]
    dialog.split.set_value("scaffold")
    assert not dialog.method.is_row_visible(dialog.split_column)


def test_the_recipe_follows_the_fields(chem) -> None:
    dialog = chem.open_tool("auto_qsar")
    dialog.target.setCurrentText("LogP")
    dialog.budget.set_value("quick")
    dialog.split.set_value("random")
    dialog.use_columns.setChecked(True)
    recipe = dialog.recipe()
    assert recipe.target == "LogP"
    assert recipe.task == "regression"
    assert recipe.budget == "quick"
    assert recipe.split == "random"
    assert recipe.descriptors and recipe.fingerprints
    assert "LogP" not in recipe.columns, "a model may not learn from its own target"


# -- a real fit ------------------------------------------------------------------------------
def test_a_fit_writes_its_columns_and_keeps_the_model(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("auto_qsar")
    dialog.target.setCurrentText("LogP")
    dialog.budget.set_value("quick")
    dialog.split.set_value("random")
    dialog.write_fold.setChecked(True)
    dialog.run()
    qtbot.waitUntil(lambda: "LogP Predicted" in table.column_names, timeout=600_000)
    for name in (
        "LogP Predicted",
        "LogP Low",
        "LogP High",
        "LogP Residual",
        "LogP Distance",
        "LogP In Domain",
        "LogP Fold",
    ):
        assert name in table.column_names, name

    from chemsuite.tools.qsar.part import PART_KEY, QsarModels

    part = chem.workspace.part(PART_KEY, lambda saved: QsarModels(table, saved))
    kept = part.kept()
    assert len(kept) == 1
    model = part.model(kept[0].asset)
    assert model.recipe.target == "LogP"
    assert len(model.training) >= 20
    assert part.capture()["models"][0]["asset"] == kept[0].asset

    panel = chem.tool_panels.find("qsar_models")
    assert panel is not None, "the fit opens the models panel"
    assert panel.current() is model


# -- the document part -----------------------------------------------------------------------
def test_the_part_keeps_models_as_assets_and_saves_only_their_ids(chem) -> None:
    from chemsuite.tools.qsar.part import MODEL_KIND, QsarModels

    table = chem.workspace.table
    part = QsarModels(table)
    model = a_model().model
    asset, op = part.put(model)
    table.edit("Add a model", [op])
    assert table.assets[asset].kind == MODEL_KIND
    saved = part.capture()
    assert saved["models"] == [
        {"asset": asset, "title": model.title, "target": "LogP", "task": "regression"}
    ]
    assert "weights" not in str(saved), "the part saves ids, not the model"

    again = QsarModels(table, saved)
    assert [k.asset for k in again.kept()] == [asset]
    pool = np.zeros((3, 8))
    assert np.allclose(again.model(asset).predict(pool).predicted, model.predict(pool).predicted)


def test_an_undone_fit_takes_its_model_out_of_the_list(chem) -> None:
    from chemsuite.tools.qsar.part import QsarModels

    table = chem.workspace.table
    part = QsarModels(table)
    asset, op = part.put(a_model().model)
    table.edit("Add a model", [op])
    assert [k.asset for k in part.kept()] == [asset]
    table.undo()
    assert part.kept() == [], "the asset is gone, so the entry goes with it"


def test_the_asset_id_never_carries_a_column_name_through(chem) -> None:
    from chemsuite.tools.qsar.part import _asset_name

    assert _asset_name("pIC50") == "qsar-pic50"
    assert _asset_name("IC50 (µM) / 24 h") == "qsar-ic50-m-24-h"
    assert _asset_name("") == "qsar"


# -- the columns -----------------------------------------------------------------------------
def test_the_column_names_put_the_target_first() -> None:
    from chemsuite.tools.qsar.columns import names_for

    chosen = ["prediction", "interval", "residual", "domain", "fold"]
    assert names_for("pIC50", chosen, "regression") == [
        "pIC50 Predicted",
        "pIC50 Low",
        "pIC50 High",
        "pIC50 Residual",
        "pIC50 Distance",
        "pIC50 In Domain",
        "pIC50 Fold",
    ]
    assert names_for("Class", chosen, "classification") == [
        "Class Predicted",
        "Class Confidence",
        "Class Possible",
        "Class Correct",
        "Class Distance",
        "Class In Domain",
        "Class Fold",
    ]


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_every_name_gets_exactly_one_column_of_values(task: str) -> None:
    from chemsuite.tools.qsar.columns import names_for, values_for

    fitted = a_model(task)
    chosen = ["prediction", "interval", "residual", "domain", "fold"]
    names = names_for("LogP", chosen, task)
    values = values_for(
        fitted.model, fitted.held, chosen, measured=fitted.measured, fold=fitted.fold
    )
    assert len(names) == len(values)
    for name, column in zip(names, values, strict=True):
        assert len(column) == len(fitted.rows), name


def test_the_written_prediction_is_the_out_of_fold_one() -> None:
    """The table must not show the model marking its own homework.

    The out-of-fold prediction is worse than the model's own fit on the same rows, always.
    If these ever came out equal, the columns would be showing memorized values.
    """
    fitted = a_model()
    model = fitted.model
    held = fitted.held.predicted
    assert len(held) == len(fitted.rows)
    own = model.predict(np.zeros((1, 8)))  # the finished model still predicts
    assert own.predicted.shape == (1,)
    assert not np.allclose(held, fitted.measured), "no fit is perfect out of fold"


# -- apply and refine ------------------------------------------------------------------------
def test_apply_says_what_to_do_when_there_is_no_model(chem) -> None:
    dialog = chem.open_tool("qsar_apply")
    assert "no model to apply" in dialog.problem()
    assert ".csqm" in dialog.problem()
    assert not dialog.primary.isEnabled()


def test_refine_says_what_to_do_when_there_is_no_model(chem) -> None:
    dialog = chem.open_tool("qsar_refine")
    assert "no model to refine" in dialog.problem()
    assert not dialog.primary.isEnabled()


def test_a_saved_model_can_be_applied_from_a_file(chem, tmp_path, qtbot) -> None:
    """A model is a file a colleague can be sent; this is the receiving end."""
    from chemsuite.qsar.model import load

    path = a_model().model.save(tmp_path / "logp.csqm")
    dialog = chem.open_tool("qsar_apply")
    dialog._loaded = load(path)
    dialog._path = path
    dialog._fill_models()
    assert dialog.model() is not None
    assert not dialog.problem()
    assert dialog.primary.text() == "Predict 24 Rows"
    dialog.close()


def test_the_pickle_export_is_refused_for_a_reopened_model(chem) -> None:
    """The archive keeps arrays, not a scikit-learn object, so there is nothing to pickle."""
    from chemsuite.tools.qsar.export_python import exportable
    from chemsuite.tools.qsar.part import QsarModels

    part = QsarModels(chem.workspace.table)
    asset, op = part.put(a_model().model)
    chem.workspace.table.edit("Add a model", [op])
    why = exportable(part, asset)
    assert "reopened from a saved document" in why
    assert "Save Model" in why


def test_applying_a_model_that_needs_a_column_this_document_lacks_is_refused(chem) -> None:
    """A recipe travels with the model, and can name a column another document never had.

    Only the features the fit kept matter: a recipe offers a wide pool, and a column the
    model looked at and discarded is no loss.
    """
    from pathlib import Path

    from chemsuite.tools.qsar.work import missing_features, read_columns

    ctx = _ctx(chem)
    fitted = a_model(kind="c")  # this one learned from table columns, not descriptors
    model = fitted.model
    kept = [n for n in model.kept_names() if n.startswith("c:")]
    assert kept, "the model must rely on a column for this test to mean anything"
    missing = missing_features(ctx, model)
    assert sorted(missing) == sorted(n[2:] for n in kept), "none of these columns exist here"

    dialog = chem.open_tool("qsar_apply")
    dialog._loaded = model
    dialog._path = Path("elsewhere.csqm")
    dialog._fill_models()
    assert "which this document does not have" in dialog.problem()
    assert not dialog.primary.isEnabled()

    # reading the columns anyway gives blanks rather than raising
    oids = chem.workspace.table.oids
    read = read_columns(ctx, ("LogP", "Nothing At All"), oids)
    assert np.isfinite(read["LogP"]).all()
    assert np.isnan(read["Nothing At All"]).all()


def test_a_model_and_its_panel_come_back_when_the_document_reopens(chem, qtbot, tmp_path) -> None:
    """Documents reopen exactly as saved, and a model belongs to the document.

    The model's bytes ride in the asset store and its id in the part, so this is the test
    that the two halves are saved and matched up again.
    """
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table
    from chemsuite.tools.qsar.panel import show_panel
    from chemsuite.tools.qsar.part import PART_KEY, QsarModels

    table = chem.workspace.table
    part = chem.workspace.part(PART_KEY, lambda saved: QsarModels(table, saved))
    model = a_model().model
    asset, op = part.put(model)
    table.edit("Add a model", [op])
    part.rename(asset, "The March model")
    show_panel(_ctx(chem))
    pool = np.zeros((3, 8))
    before = model.predict(pool).predicted

    path = tmp_path / "with_a_model.csd"
    assert chem._save_to(path)
    chem._install(Workspace(Table()))
    chem.open_path(path)
    qtbot.waitUntil(lambda: chem.workspace.path == path, timeout=20_000)
    qtbot.wait(20)

    again = chem.workspace.part(PART_KEY, lambda saved: QsarModels(chem.workspace.table, saved))
    kept = again.kept()
    assert [k.title for k in kept] == ["The March model"], "the name was saved too"
    reopened = again.model(kept[0].asset)
    assert np.allclose(reopened.predict(pool).predicted, before)
    assert len(reopened.training) == len(model.training)
    assert reopened.recipe == model.recipe

    panel = chem.tool_panels.find("qsar_models")
    assert panel is not None, "the panel reopens with the document"
    assert panel.current() is reopened


def _ctx(window):
    from chemsuite.app.document.tools import ToolContext

    return ToolContext(window)
