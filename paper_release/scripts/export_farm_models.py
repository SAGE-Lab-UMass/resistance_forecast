"""Serialize the manuscript models into loadable artifacts for the ``farm`` package.

The training logic is not reimplemented here. This script reads
``model_essential_nonessential_combined_mutation_level.ipynb`` -- the notebook that
produces every number reported in the manuscript -- pulls the cells that define the
label construction and ``run_experiment``, and executes them. The fitted estimators
that come out are the same objects the notebook holds in ``trained_models``.

After fitting, every exported model is checked against the prediction CSVs the
notebook already wrote to ``results/training_predictions/``. The export fails if any
holdout probability differs, so an artifact can never silently drift from the paper.

Usage:
    python scripts/export_farm_models.py [--out ../farm/artifacts]
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

RELEASE_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = RELEASE_ROOT / "model_essential_nonessential_combined_mutation_level.ipynb"

# dataset -> model family selected in the manuscript, and the notebook key it lands under
SELECTED = {
    "Essential": "rf",
    "Nonessential": "logistic",
    "Combined": "rf",
}

# tolerance for the reproduction check; predictions should match to float precision
ATOL = 1e-12


def _code_cells(notebook: Path) -> list[str]:
    nb = json.loads(notebook.read_text())
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


def _cell_containing(cells: list[str], marker: str) -> str:
    for src in cells:
        if marker in src:
            return src
    raise LookupError(f"no notebook cell contains {marker!r}")


def run_notebook_training() -> dict:
    """Execute the notebook's own training path and return its namespace."""
    cells = _code_cells(NOTEBOOK)
    ns: dict = {"__name__": "__farm_export__"}

    # The import cell chdirs to the release root and resolves FEATURES_ROOT.
    exec(_cell_containing(cells, "def find_run_root"), ns)

    # Cell 3 equivalent: load and filter the 2021 feature table.
    exec(_cell_containing(cells, "catalog_data = pd.read_csv"), ns)

    # create_label_data + the call that builds X_ess/X_noness/X_combined/label_data.
    exec(_cell_containing(cells, "def create_label_data"), ns)

    # run_experiment, then the dataset dict, then the training loop.
    exec(_cell_containing(cells, "def run_experiment"), ns)
    exec(_cell_containing(cells, '"Essential": (' ), ns)

    # The training loop writes prediction CSVs relative to the cwd. Redirect them to a
    # scratch directory so the repository's committed results are never overwritten.
    train_cell = _cell_containing(cells, "trained_models = {}")
    with tempfile.TemporaryDirectory() as tmp:
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            exec(train_cell, ns)
        finally:
            os.chdir(cwd)
    return ns


def verify(ns: dict, dataset: str, family: str) -> dict:
    """Re-score the notebook's saved holdout rows and require an exact match."""
    reference = RELEASE_ROOT / "results" / "training_predictions" / f"{family}_{dataset}_holdout_predictions.csv"
    if not reference.exists():
        raise FileNotFoundError(f"missing notebook reference predictions: {reference}")

    ref = pd.read_csv(reference)
    model = ns["trained_models"][f"{dataset}_{family}"]
    features = list(ns["numeric_columns"])

    probs = model.predict_proba(ref[features])[:, 1]
    delta = np.abs(probs - ref["pred_prob"].to_numpy())
    if not np.all(delta <= ATOL):
        raise AssertionError(
            f"{dataset}/{family}: exported model disagrees with the notebook "
            f"(max |delta| = {delta.max():.3e} over {len(delta)} holdout rows)"
        )
    return {"n_holdout_checked": int(len(delta)), "max_abs_prob_delta": float(delta.max())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=RELEASE_ROOT.parent / "farm" / "artifacts",
        help="directory to write the serialized models into",
    )
    args = parser.parse_args()

    import joblib  # imported late so --help works without it

    print(f"executing training cells from {NOTEBOOK.name} ...")
    ns = run_notebook_training()

    features = list(ns["numeric_columns"])
    summary = pd.DataFrame(ns["training_results"])
    label_data = ns["label_data"]

    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict = {
        "created": date.today().isoformat(),
        "source_notebook": NOTEBOOK.name,
        "feature_columns": features,
        "n_features": len(features),
        "models": {},
    }

    for dataset, family in SELECTED.items():
        model = ns["trained_models"][f"{dataset}_{family}"]
        row = summary[
            (summary["dataset"] == dataset)
            & (summary["model"] == family)
            & (summary["mode"] == "holdout")
        ]
        if len(row) != 1:
            raise AssertionError(f"expected exactly one holdout row for {dataset}/{family}, got {len(row)}")
        row = row.iloc[0]

        check = verify(ns, dataset, family)
        print(f"  {dataset:<13} {family:<9} threshold={row['threshold']:.3f}  "
              f"holdout AUC={row['auc']:.3f}  verified on {check['n_holdout_checked']} rows")

        key = dataset.lower()
        joblib.dump(model, out_dir / f"{key}.joblib", compress=3)

        if dataset == "Combined":
            n_train = int(len(label_data))
        else:
            n_train = int((label_data["essentiality"] == key).sum())

        manifest["models"][key] = {
            "dataset": dataset,
            "family": "RandomForestClassifier" if family == "rf" else "LogisticRegression",
            "artifact": f"{key}.joblib",
            "threshold": float(row["threshold"]),
            "holdout": {
                "auc": float(row["auc"]),
                "sensitivity": float(row["sensitivity"]),
                "specificity": float(row["specificity"]),
                "f1_resistant": float(row["f1_res"]),
                "n_test": int(row["n_test"]),
                "n_train": int(row["n_train"]),
            },
            "n_labeled_mutations": n_train,
            "verification": check,
        }

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nwrote {len(manifest['models'])} models + manifest.json to {out_dir}")


if __name__ == "__main__":
    main()
