"""The 25-feature representation the FARM models expect."""

from __future__ import annotations

from typing import Iterable

import pandas as pd

#: Feature columns in the exact order the models were fitted on. Order matters:
#: scikit-learn estimators are positional, so a reordered frame silently produces
#: wrong predictions rather than an error.
FEATURE_COLUMNS: tuple[str, ...] = (
    "Rosetta_fa_atr",
    "Rosetta_fa_rep",
    "Rosetta_fa_sol",
    "Rosetta_fa_elec",
    "Rosetta_fa_dun",
    "Rosetta_ddG",
    "Prox_3D_zeroed",
    "LLR_dim33",
    "LLR_dim46",
    "LLR_dim62",
    "LLR_dim70",
    "LLR_dim124",
    "LLR_dim192",
    "LLR_dim207",
    "LLR_dim258",
    "LLR_dim267",
    "LLR_dim315",
    "AAIndex_mut1",
    "AAIndex_mut2",
    "AAIndex_mut3",
    "AAIndex_mut4",
    "AAIndex_mut5",
    "AAIndex_mut6",
    "AAIndex_mut7",
    "AAIndex_mut8",
)

#: Which modality each feature belongs to, for reporting and ablation.
FEATURE_GROUPS: dict[str, tuple[str, ...]] = {
    "biophysical_energy": FEATURE_COLUMNS[0:6],
    "structural_context": FEATURE_COLUMNS[6:7],
    "protein_language_model": FEATURE_COLUMNS[7:17],
    "aaindex": FEATURE_COLUMNS[17:25],
}


class MissingFeaturesError(ValueError):
    """Raised when an input frame does not carry the full feature representation."""

    def __init__(self, missing: Iterable[str]) -> None:
        self.missing = list(missing)
        shown = ", ".join(self.missing[:8])
        if len(self.missing) > 8:
            shown += f", ... (+{len(self.missing) - 8} more)"
        super().__init__(
            f"input is missing {len(self.missing)} of the {len(FEATURE_COLUMNS)} required "
            f"feature columns: {shown}. These are produced by the feature pipeline in "
            f"paper_release/forecast_data_preparation_combined.ipynb; FARM scores precomputed "
            f"features and does not compute them from a mutation string."
        )


def prepare_features(data: pd.DataFrame) -> pd.DataFrame:
    """Select and order the feature columns, filling missing values with zero.

    The zero fill matches the manuscript's deployment path, where variants keep a
    feature vector even when an individual descriptor could not be computed. A frame
    that is missing whole columns is an error rather than something to fill, because
    that means the feature pipeline was not run.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError(f"expected a pandas DataFrame, got {type(data).__name__}")

    missing = [c for c in FEATURE_COLUMNS if c not in data.columns]
    if missing:
        raise MissingFeaturesError(missing)

    return data.loc[:, list(FEATURE_COLUMNS)].astype(float).fillna(0.0)
