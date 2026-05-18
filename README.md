# Resistance Forecasting in *Mycobacterium tuberculosis*

This repository contains the data preparation, multimodal feature engineering, and model-training workflows used to forecast antibiotic resistance-associated missense variants in *Mycobacterium tuberculosis*.

The central goal of the project is to prioritize WHO catalogue variants of uncertain significance by asking a temporally realistic question: can a model trained on the 2021 WHO mutation catalogue identify uncertain variants that later become recognized as resistance-associated in the 2023 catalogue, and can that same model be deployed to prioritize the remaining 2023 uncertain variants for follow-up?

The workflow combines structural context, Rosetta-derived energetics, AAIndex-based physicochemical descriptors, and ESM-2 sequence likelihood features to model distinct mutational regimes across essential, nonessential, and combined gene sets.

## Project overview

The repository currently contains two layers:

1. the broader project history (`Comparison_Model`, `legacy_files`, earlier exploratory assets), and
2. a cleaned release package under `paper_release/` containing the final manuscript-facing code and data resources.

If you are using this repository for reproduction of the manuscript analyses, start with `paper_release/`.

## Scientific objective

The primary analysis focuses on WHO resistance confidence categories:

| Category | Description |
|----------|-------------|
| **1** | Associated with resistance |
| **2** | Associated with resistance -- interim |
| **3** | Uncertain significance |
| **4** | Not associated with resistance -- interim |
| **5** | Not associated with resistance |

The main forecasting workflow:
- trains on 2021 known-effect variants (Categories 1/2 as resistant and 4/5 as susceptible),
- evaluates retrospectively on the subset of 2021 Category 3 variants that were reclassified in 2023, and
- deploys the frozen Combined model to all feature-complete 2023 uncertain-significance variants.

## Feature representation

Each missense variant is represented using a multimodal feature set that combines:

- **Structural context**: three-dimensional proximity to known resistance-associated positions
- **Rosetta energetics**: mutation-associated energy and stability terms
- **Expanded ESM-2 likelihood features**: sequence-context plausibility from a protein language model
- **AAIndex descriptors**: reduced physicochemical encoding of the mutant residue

These features are used to capture different biological resistance regimes, including structurally constrained target-gene mutations and broader disruption-like mutational patterns in nonessential or accessory genes.

## Repository structure

```text
resistance_forecast/
├── Comparison_Model/                # Alternative statistical-model analyses
├── data/                            # Project-level source data and legacy resource store
├── paper_release/                   # Final manuscript-facing code and curated data package
│   ├── forecast_data_preparation_combined.ipynb
│   ├── model_essential_nonessential_combined.ipynb
│   ├── forecast_utils.py
│   ├── requirements.txt
│   ├── scripts/
│   ├── source_data/
│   │   ├── catalog/
│   │   ├── derived_features/
│   │   └── distmaps/
│   └── supplementary_data/
├── legacy_files/                    # Earlier notebooks and superseded scripts
├── mycobrowser.py                   # Helper for Mycobrowser-derived metadata
├── requirements.txt
└── README.md
```

## Recommended entry point

Use the `paper_release/` folder for the final workflow.

### Main notebooks

- `paper_release/forecast_data_preparation_combined.ipynb`
  - Builds the 2021 and 2023 feature-complete variant tables.
- `paper_release/model_essential_nonessential_combined.ipynb`
  - Trains the selected models, evaluates holdout and temporal forecasting performance, and generates the main manuscript outputs.

### Shared module

- `paper_release/forecast_utils.py`
  - Utility code used by both notebooks for preprocessing, feature assembly, and model-support functions.

### Post-processing scripts

- `paper_release/scripts/fp_fn_feature_group_summary.py`
- `paper_release/scripts/gene_misclassification_chisquare.py`

## Curated release data

The curated manuscript-facing data package lives under `paper_release/source_data/`.

### Included source data

- `paper_release/source_data/catalog/`
  - WHO 2021 and WHO 2023 catalogues, protein metadata, and AAIndex reference tables
- `paper_release/source_data/distmaps/`
  - structural distance-map resources used for proximity features
- `paper_release/source_data/derived_features/`
  - final 2021 and 2023 derived feature tables used by the modeling notebook

Additional structural resources are retained in the project-level `data/` directory.

## Supplementary deployment data

The full set of 5,215 Combined-model forecasts for feature-complete WHO 2023 uncertain-significance variants is provided in:

- `paper_release/supplementary_data/Supplementary_Data_1_Combined_2023_uncertain_variant_forecasts.csv`

This file contains the ranked deployment-style forecast output used in the manuscript.

## Quick start

Clone the repository and move into the project root:

```bash
git clone https://github.com/SAGE-Lab-UMass/resistance_forecast.git
cd resistance_forecast
```

Install dependencies for the release workflow:

```bash
pip install -r paper_release/requirements.txt
```

Then run the notebooks in order:

```text
1. paper_release/forecast_data_preparation_combined.ipynb
2. paper_release/model_essential_nonessential_combined.ipynb
```

## Main outputs

The release workflow produces:
- 2021 and 2023 final feature tables
- holdout and temporal evaluation summaries
- per-model prediction tables for reclassified Category 3 variants
- 2023 deployment forecasts for uncertain-significance variants
- figure-ready performance and interpretability outputs

## Contact

For questions or collaborations, please contact:
**[mtasmin@umass.edu](mailto:mtasmin@umass.edu)**
