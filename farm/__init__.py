"""FARM -- Forecasting Antibiotic Resistance in *Mycobacterium tuberculosis*.

Inference interface for the models described in the FARM manuscript. The fitted
estimators shipped here are serialized directly from
``paper_release/model_essential_nonessential_combined_mutation_level.ipynb`` and are
verified against that notebook's own holdout predictions at export time.

    >>> from farm import predict
    >>> scored = predict("variants_with_features.csv")        # doctest: +SKIP

Scoring expects the 25-feature multimodal representation (see
:data:`farm.FEATURE_COLUMNS`) to have been computed already by the feature pipeline in
``paper_release/forecast_data_preparation_combined.ipynb``. FARM does not build
features from a bare mutation string, because the structural context and Rosetta terms
require per-gene structures and distance maps.
"""

from .features import FEATURE_COLUMNS, FEATURE_GROUPS, MissingFeaturesError, prepare_features
from .model import FarmModel, available_models, load_model
from .predict import predict

__all__ = [
    "predict",
    "load_model",
    "available_models",
    "FarmModel",
    "FEATURE_COLUMNS",
    "FEATURE_GROUPS",
    "prepare_features",
    "MissingFeaturesError",
    "__version__",
]

__version__ = "1.0.0"
