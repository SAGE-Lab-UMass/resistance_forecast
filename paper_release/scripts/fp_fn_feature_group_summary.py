#!/usr/bin/env python3
"""Summarize FP/FN feature-group shifts using the real training feature set.

This script computes model-level summaries comparing:
  - FN vs TP
  - FP vs TN

for the four biological groups actually used in training:
  - Forcefield_Energy: Rosetta_* columns
  - Structural_Context: Prox_3D_zeroed
  - Expanded_LLR: LLR_dim* columns only
  - Mutational_AAIndex: AAIndex_mut* columns only

Important: this excludes LLR_score and AAIndex_delta* because they were not
included in the supervised training matrix.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


RESULTS_DIR = (
    Path(__file__).resolve().parents[1]
    / "results"
)
ERROR_DIR = RESULTS_DIR / "error_analysis"

PREDICTION_FILES = {
    "Combined": RESULTS_DIR / "Combined_RF_predictions_on_category3.csv",
    "Essential": RESULTS_DIR / "Essential_RF_predictions_on_category3.csv",
    "Nonessential": RESULTS_DIR / "Nonessential_LOGISTIC_predictions_on_category3.csv",
}

GROUPS = {
    "Forcefield_Energy": lambda c: c.startswith("Rosetta_"),
    "Structural_Context": lambda c: c == "Prox_3D_zeroed",
    "Expanded_LLR": lambda c: c.startswith("LLR_dim"),
    "Mutational_AAIndex": lambda c: c.startswith("AAIndex_mut"),
}


def build_group_features(df: pd.DataFrame) -> dict[str, list[str]]:
    group_features: dict[str, list[str]] = {}
    for group, rule in GROUPS.items():
        cols = [c for c in df.columns if rule(c)]
        if not cols:
            raise ValueError(f"No columns found for group {group}")
        group_features[group] = cols
    return group_features


def add_outcome_labels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    y_true = out["true_label"].astype(int)
    y_pred = out["model_prediction"].astype(int)

    conditions = [
        (y_true == 1) & (y_pred == 1),
        (y_true == 1) & (y_pred == 0),
        (y_true == 0) & (y_pred == 0),
        (y_true == 0) & (y_pred == 1),
    ]
    labels = ["TP", "FN", "TN", "FP"]
    out["outcome"] = np.select(conditions, labels, default="UNKNOWN")
    return out


def zscore_columns(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    x = df[cols].apply(pd.to_numeric, errors="coerce")
    mu = x.mean(axis=0)
    sigma = x.std(axis=0, ddof=0).replace(0, np.nan)
    return (x - mu) / sigma


def summarize_model(model_label: str, path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path)
    df = add_outcome_labels(df)

    group_features = build_group_features(df)

    # Raw structural context retains its physical interpretation.
    df["Structural_Context_raw"] = pd.to_numeric(df["Prox_3D_zeroed"], errors="coerce")

    # Standardized grouped scores for the remaining training features.
    summary_cols = []
    for group, cols in group_features.items():
        if group == "Structural_Context":
            continue
        z = zscore_columns(df, cols)
        df[f"{group}_zscore_mean_abs"] = z.abs().mean(axis=1, skipna=True)
        summary_cols.append(f"{group}_zscore_mean_abs")

    rows = []
    for group in ["Structural_Context", "Forcefield_Energy", "Expanded_LLR", "Mutational_AAIndex"]:
        value_col = (
            "Structural_Context_raw"
            if group == "Structural_Context"
            else f"{group}_zscore_mean_abs"
        )

        medians = (
            df.groupby("outcome")[value_col]
            .median()
            .reindex(["TP", "FN", "TN", "FP"])
        )
        counts = (
            df["outcome"]
            .value_counts()
            .reindex(["TP", "FN", "TN", "FP"])
            .fillna(0)
            .astype(int)
        )

        rows.append(
            {
                "model": model_label,
                "group": group,
                "value_type": (
                    "raw_median" if group == "Structural_Context" else "standardized_mean_abs_z_median"
                ),
                "n_TP": counts["TP"],
                "n_FN": counts["FN"],
                "n_TN": counts["TN"],
                "n_FP": counts["FP"],
                "median_TP": medians["TP"],
                "median_FN": medians["FN"],
                "delta_FN_minus_TP": medians["FN"] - medians["TP"],
                "median_TN": medians["TN"],
                "median_FP": medians["FP"],
                "delta_FP_minus_TN": medians["FP"] - medians["TN"],
            }
        )

    summary = pd.DataFrame(rows)

    detail_cols = [
        "gene",
        "drug",
        "mutation_oneletter",
        "true_label",
        "model_prediction",
        "outcome",
        "Structural_Context_raw",
    ] + summary_cols
    detail = df[detail_cols].copy()
    detail.insert(0, "model", model_label)
    return summary, detail


def main() -> None:
    ERROR_DIR.mkdir(parents=True, exist_ok=True)

    summaries = []
    details = []
    for model_label, path in PREDICTION_FILES.items():
        summary, detail = summarize_model(model_label, path)
        summaries.append(summary)
        details.append(detail)

    summary_df = pd.concat(summaries, ignore_index=True)
    detail_df = pd.concat(details, ignore_index=True)

    summary_path = ERROR_DIR / "fp_fn_feature_group_summary.csv"
    detail_path = ERROR_DIR / "fp_fn_feature_group_variant_scores.csv"

    summary_df.to_csv(summary_path, index=False)
    detail_df.to_csv(detail_path, index=False)

    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nSaved: {summary_path}")
    print(f"Saved: {detail_path}")


if __name__ == "__main__":
    main()
