# TB Resistance Forecasting Code Release

This folder contains the final analysis workflow used for the manuscript-level forecasting results. 

## Included files

- `forecast_data_preparation_combined.ipynb`
  - Main notebook for dataset construction and multimodal feature assembly from the WHO 2021 and 2023 catalogues.
- `model_essential_nonessential_combined.ipynb`
  - Main notebook for model training, holdout evaluation, temporal reclassification analysis, 2023 forecasting, and figure generation.
- `forecast_utils.py`
  - Shared utility functions used by the notebooks.
- `scripts/fp_fn_feature_group_summary.py`
  - Post-processing summary of false-negative and false-positive feature-group shifts.
- `scripts/gene_misclassification_chisquare.py`
  - Gene-level enrichment analysis for residual errors.
- `requirements.txt`
  - Python dependencies used for the final workflow.

## Expected repository layout

This release assumes the project data directory is available as `data/` at the repository root. The code can be run either:

1. directly from this folder when it sits inside the repository root, or
2. from a repository root that contains both this code and the required `data/` directory.

Derived features are written to:

- `data/derived_features/` when a local writable `data/` tree is available inside the code root, or
- the parent project `data/derived_features/` directory when the code is nested within a larger repository layout.

Model outputs, evaluation summaries, figures, and forecast tables are written to:

- `results/`

## Run order

Install dependencies:

```bash
pip install -r requirements.txt
```

Then run the notebooks in this order:

1. `forecast_data_preparation_combined.ipynb`
2. `model_essential_nonessential_combined.ipynb`

## Notes on the final workflow

The release version preserves the key fixes introduced in the isolated review run:

- holdout and Category-3 evaluation thresholds are derived from training data only
- cross-validation no longer tunes thresholds on fold test labels
- path detection is portable and no longer depends on hardcoded absolute paths
- essential/nonessential evaluation subsets are derived from training-supported genes
- the ESM contribution utilities include the residue-indexing and empty-array guards used in the final analysis
