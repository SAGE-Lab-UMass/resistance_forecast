"""Loading and scoring with the serialized FARM models."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .features import FEATURE_COLUMNS, prepare_features

ARTIFACT_DIR = Path(__file__).resolve().parent / "artifacts"
MANIFEST_PATH = ARTIFACT_DIR / "manifest.json"

#: Accepted spellings for each model, so callers can use the manuscript's wording.
_ALIASES = {
    "combined": "combined",
    "all": "combined",
    "essential": "essential",
    "ess": "essential",
    "nonessential": "nonessential",
    "non-essential": "nonessential",
    "noness": "nonessential",
}

#: Columns carried through to the output frame when present on the input.
_PASSTHROUGH = ("gene", "mutation_oneletter", "mutation_pos", "drug", "confidence", "phenotype")


@lru_cache(maxsize=1)
def _manifest() -> dict:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"no model manifest at {MANIFEST_PATH}. Regenerate the artifacts with "
            f"`python paper_release/scripts/export_farm_models.py`."
        )
    return json.loads(MANIFEST_PATH.read_text())


def available_models() -> list[str]:
    """Names that :func:`load_model` accepts."""
    return sorted(_manifest()["models"])


def _resolve(name: str) -> str:
    key = _ALIASES.get(str(name).strip().lower().replace("_", "-"))
    if key is None:
        key = _ALIASES.get(str(name).strip().lower())
    if key is None:
        raise KeyError(f"unknown model {name!r}; available: {', '.join(available_models())}")
    return key


@dataclass(frozen=True)
class FarmModel:
    """A fitted FARM classifier together with the threshold it was calibrated with.

    The threshold is the median Youden index over repeated stratified cross-validation
    *within the training split only*, exactly as in the manuscript, so scoring new
    variants never consults their labels.
    """

    name: str
    estimator: object
    threshold: float
    metadata: dict

    @property
    def feature_columns(self) -> tuple[str, ...]:
        return FEATURE_COLUMNS

    def predict_proba(self, data: pd.DataFrame) -> np.ndarray:
        """Resistance-association probability for each row."""
        X = prepare_features(data)
        return self.estimator.predict_proba(X)[:, 1]

    def predict(self, data: pd.DataFrame, threshold: float | None = None) -> np.ndarray:
        """Binary resistance call at the calibrated threshold (1 = resistance-associated)."""
        cut = self.threshold if threshold is None else float(threshold)
        return (self.predict_proba(data) >= cut).astype(int)

    def predict_frame(self, data: pd.DataFrame, threshold: float | None = None) -> pd.DataFrame:
        """Predictions as a frame, carrying identifier columns through when present."""
        cut = self.threshold if threshold is None else float(threshold)
        probs = self.predict_proba(data)

        carried = [c for c in _PASSTHROUGH if c in data.columns]
        out = data.loc[:, carried].copy() if carried else pd.DataFrame(index=data.index)
        out["farm_model"] = self.name
        out["resistance_score"] = probs
        out["threshold"] = cut
        out["predicted_resistant"] = (probs >= cut).astype(int)
        return out.reset_index(drop=True)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        h = self.metadata.get("holdout", {})
        return (
            f"FarmModel(name={self.name!r}, family={self.metadata.get('family')}, "
            f"threshold={self.threshold:.3f}, holdout_auc={h.get('auc')})"
        )


@lru_cache(maxsize=8)
def load_model(name: str = "combined") -> FarmModel:
    """Load one of the manuscript's fitted models.

    ``combined`` is the deployment model used to forecast the 2023 uncertain-significance
    variants; ``essential`` and ``nonessential`` are the gene-stratum models.
    """
    import joblib

    key = _resolve(name)
    meta = _manifest()["models"][key]
    path = ARTIFACT_DIR / meta["artifact"]
    if not path.exists():
        raise FileNotFoundError(
            f"missing artifact {path}. Regenerate with "
            f"`python paper_release/scripts/export_farm_models.py`."
        )
    return FarmModel(
        name=key,
        estimator=joblib.load(path),
        threshold=float(meta["threshold"]),
        metadata=meta,
    )
