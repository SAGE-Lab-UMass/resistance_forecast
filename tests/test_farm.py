"""Tests for the ``farm`` inference package.

The important test is :func:`test_matches_notebook_holdout_predictions`: it re-scores
the exact rows the manuscript notebook held out and requires the packaged model to
reproduce its probabilities. If that passes, the package and the paper agree.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import farm

REPO = Path(__file__).resolve().parents[1]
PREDICTIONS = REPO / "paper_release" / "results" / "training_predictions"

# (model name, notebook prediction file, expected holdout AUC from the manuscript)
CASES = [
    ("combined", "rf_Combined_holdout_predictions.csv", 0.920),
    ("essential", "rf_Essential_holdout_predictions.csv", 0.843),
    ("nonessential", "logistic_Nonessential_holdout_predictions.csv", 0.943),
]


@pytest.fixture(scope="module")
def sample() -> pd.DataFrame:
    return pd.read_csv(PREDICTIONS / "rf_Combined_holdout_predictions.csv")


def test_feature_contract():
    assert len(farm.FEATURE_COLUMNS) == 25
    flat = [c for cols in farm.FEATURE_GROUPS.values() for c in cols]
    assert flat == list(farm.FEATURE_COLUMNS), "groups must tile the feature vector in order"


def test_available_models():
    assert set(farm.available_models()) == {"combined", "essential", "nonessential"}


@pytest.mark.parametrize("name,filename,expected_auc", CASES)
def test_matches_notebook_holdout_predictions(name, filename, expected_auc):
    """The packaged model must reproduce the notebook's holdout probabilities."""
    ref = pd.read_csv(PREDICTIONS / filename)
    model = farm.load_model(name)

    probs = model.predict_proba(ref)
    np.testing.assert_allclose(probs, ref["pred_prob"].to_numpy(), atol=1e-12)

    calls = model.predict(ref)
    np.testing.assert_array_equal(calls, ref["pred_label"].to_numpy())

    from sklearn.metrics import roc_auc_score

    assert round(roc_auc_score(ref["true_label"], probs), 3) == expected_auc


@pytest.mark.parametrize("name,filename,_auc", CASES)
def test_threshold_matches_manifest(name, filename, _auc):
    model = farm.load_model(name)
    assert model.threshold == pytest.approx(model.metadata["threshold"])
    assert 0.0 < model.threshold < 1.0


def test_predict_frame_shape_and_columns(sample):
    out = farm.predict(sample)
    assert len(out) == len(sample)
    for col in ("farm_model", "resistance_score", "threshold", "predicted_resistant"):
        assert col in out.columns
    assert out["resistance_score"].between(0, 1).all()
    assert set(out["predicted_resistant"].unique()) <= {0, 1}
    # identifier columns present on the input are carried through
    assert "gene" in out.columns


def test_predict_accepts_a_path(tmp_path, sample):
    path = tmp_path / "variants.csv"
    sample.to_csv(path, index=False)
    from_path = farm.predict(path)
    from_frame = farm.predict(sample)
    pd.testing.assert_frame_equal(from_path, from_frame)


def test_predict_returns_array_when_asked(sample):
    arr = farm.predict(sample, return_frame=False)
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (len(sample),)


def test_threshold_override_changes_calls(sample):
    strict = farm.predict(sample, threshold=0.99)["predicted_resistant"].sum()
    loose = farm.predict(sample, threshold=0.01)["predicted_resistant"].sum()
    assert strict < loose == len(sample)


def test_column_order_is_enforced(sample):
    """A reordered input must not change the answer."""
    shuffled = sample[list(reversed(sample.columns))]
    np.testing.assert_allclose(
        farm.predict(shuffled, return_frame=False),
        farm.predict(sample, return_frame=False),
        atol=1e-12,
    )


def test_missing_features_raise_a_useful_error(sample):
    broken = sample.drop(columns=["Prox_3D_zeroed", "Rosetta_ddG"])
    with pytest.raises(farm.MissingFeaturesError) as exc:
        farm.predict(broken)
    assert "Prox_3D_zeroed" in str(exc.value)
    assert set(exc.value.missing) == {"Prox_3D_zeroed", "Rosetta_ddG"}


def test_nan_features_are_zero_filled(sample):
    with_nan = sample.copy()
    with_nan.loc[0, "Rosetta_ddG"] = np.nan
    zeroed = sample.copy()
    zeroed.loc[0, "Rosetta_ddG"] = 0.0
    np.testing.assert_allclose(
        farm.predict(with_nan, return_frame=False),
        farm.predict(zeroed, return_frame=False),
        atol=1e-12,
    )


def test_unknown_model_name():
    with pytest.raises(KeyError):
        farm.load_model("does-not-exist")


def test_cli_predict(tmp_path, sample):
    from farm.cli import main

    src, dst = tmp_path / "in.csv", tmp_path / "out.csv"
    sample.to_csv(src, index=False)
    assert main(["predict", str(src), "-o", str(dst), "-m", "combined"]) == 0

    out = pd.read_csv(dst)
    assert len(out) == len(sample)
    np.testing.assert_allclose(
        out["resistance_score"].to_numpy(), sample["pred_prob"].to_numpy(), atol=1e-12
    )


def test_cli_template_round_trips(tmp_path):
    from farm.cli import main

    dst = tmp_path / "template.csv"
    assert main(["template", "-o", str(dst)]) == 0
    cols = pd.read_csv(dst).columns.tolist()
    assert all(c in cols for c in farm.FEATURE_COLUMNS)


def test_reproduces_the_deployment_forecast():
    """End-to-end: the packaged model must reproduce the manuscript's 696 of 4,525."""
    path = (REPO / "paper_release" / "results" / "forecast_2023_unclassified"
            / "Combined_on_full_2023_uncertain.csv")
    if not path.exists():
        pytest.skip("deployment forecast table not present in this checkout")

    ref = pd.read_csv(path)
    out = farm.predict(ref, model="combined")

    assert len(out) == 4525
    np.testing.assert_allclose(
        out["resistance_score"].to_numpy(), ref["pred_prob_res"].to_numpy(), atol=1e-12
    )
    np.testing.assert_array_equal(
        out["predicted_resistant"].to_numpy(), ref["pred_label"].to_numpy()
    )
    assert int(out["predicted_resistant"].sum()) == 696
