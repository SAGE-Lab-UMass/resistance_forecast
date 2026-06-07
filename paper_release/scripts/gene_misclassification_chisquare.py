#!/usr/bin/env python3
"""Gene-level chi-square tests for Category 3 forecast misclassifications.

Tests whether error status (correct vs incorrect) is independent of gene for each
selected model's reclassified Category 3 prediction file. Also reports error
enrichment: observed incorrect calls divided by the expected incorrect calls
under proportional allocation by gene sample size.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


RESULTS_DIR = (
    Path(__file__).resolve().parents[1]
    / "results"
)

PREDICTION_FILES = {
    "Combined RF": RESULTS_DIR / "Combined_RF_predictions_on_category3.csv",
    "Nonessential Logistic": RESULTS_DIR / "Nonessential_LOGISTIC_predictions_on_category3.csv",
    "Essential RF": RESULTS_DIR / "Essential_RF_predictions_on_category3.csv",
}


def monte_carlo_pvalue(obs: np.ndarray, n_iter: int = 200_000, seed: int = 1) -> float:
    """Fixed-margin Monte Carlo p-value for chi-square statistic.

    This preserves gene row totals and total incorrect calls, randomly allocating
    incorrect labels across all variants under the null that errors are
    proportional to gene sample size.
    """
    rng = np.random.default_rng(seed)
    row_totals = obs.sum(axis=1)
    n_total = int(row_totals.sum())
    n_incorrect = int(obs[:, 1].sum())

    observed_chi2 = chi2_contingency(obs, correction=False)[0]
    gene_ids = np.repeat(np.arange(len(row_totals)), row_totals)

    n_ge = 0
    for _ in range(n_iter):
        error_indices = rng.choice(n_total, size=n_incorrect, replace=False)
        simulated_incorrect = np.bincount(gene_ids[error_indices], minlength=len(row_totals))
        simulated_obs = np.column_stack([row_totals - simulated_incorrect, simulated_incorrect])
        simulated_chi2 = chi2_contingency(simulated_obs, correction=False)[0]
        if simulated_chi2 >= observed_chi2 - 1e-12:
            n_ge += 1

    return (n_ge + 1) / (n_iter + 1)


def analyze_model(model_name: str, path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "correct" not in df.columns:
        raise ValueError(f"Missing required 'correct' column in {path}")
    if "gene" not in df.columns:
        raise ValueError(f"Missing required 'gene' column in {path}")

    df["Incorrect"] = ~df["correct"].astype(bool)

    rows = []
    for gene, group in df.groupby("gene", sort=True):
        incorrect = int(group["Incorrect"].sum())
        total = int(len(group))
        correct = total - incorrect
        rows.append(
            {
                "Gene": gene,
                "Correct": correct,
                "Incorrect": incorrect,
                "Total": total,
                "ErrorRate": incorrect / total,
            }
        )

    table = pd.DataFrame(rows).set_index("Gene")
    observed = table[["Correct", "Incorrect"]].to_numpy()

    chi2, p_asymptotic, dof, expected = chi2_contingency(observed, correction=False)
    expected_df = pd.DataFrame(
        expected,
        index=table.index,
        columns=["ExpectedCorrect", "ExpectedIncorrect"],
    )

    table["ExpectedIncorrect"] = expected_df["ExpectedIncorrect"]
    table["Enrichment"] = table["Incorrect"] / table["ExpectedIncorrect"].replace(0, np.nan)

    p_monte_carlo = monte_carlo_pvalue(observed)

    print(f"\n=== {model_name} ===")
    print(table.to_string(float_format=lambda x: f"{x:.3f}"))
    print(
        f"chi2={chi2:.4f}, dof={dof}, "
        f"asymptotic_p={p_asymptotic:.5f}, "
        f"monte_carlo_p={p_monte_carlo:.5f}, "
        f"min_expected={expected.min():.3f}"
    )

    out = table.reset_index()
    out.insert(0, "Model", model_name)
    out["Chi2"] = chi2
    out["Dof"] = dof
    out["AsymptoticP"] = p_asymptotic
    out["MonteCarloP"] = p_monte_carlo
    out["MinExpected"] = expected.min()
    return out


def main() -> None:
    all_rows = []
    for model_name, path in PREDICTION_FILES.items():
        if not path.exists():
            raise FileNotFoundError(path)
        all_rows.append(analyze_model(model_name, path))

    output_path = RESULTS_DIR / "gene_misclassification_chisquare.csv"
    pd.concat(all_rows, ignore_index=True).to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    main()
