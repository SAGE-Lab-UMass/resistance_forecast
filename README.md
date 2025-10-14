# Resistance Forecast in *Mycobacterium tuberculosis* Antibiotics

This repository contains data preparation and machine learning workflows for predicting the confidence level of antibiotic resistance–associated mutations in *Mycobacterium tuberculosis* (M. tuberculosis). By integrating biological, structural, and protein language model–derived features, this project aims to classify uncertain mutations based on their likelihood of conferring drug resistance.

---

##  Objective

The primary goal is to predict the **confidence category** for mutations in *M. tuberculosis* genes labeled as having “uncertain significance” in the WHO 2021 catalogue and verify with the updated confidence categories from the WHO 2023 catalogue. These categories reflect the strength of evidence linking a mutation to antibiotic resistance:

| Category | Description |
|-----------|--------------|
| **1** | Associated with resistance |
| **2** | Associated with resistance – interim |
| **3** | Uncertain significance |
| **4** | Not associated with resistance – interim |
| **5** | Not associated with resistance |



##  Dataset Overview

- Mutations from **17 essential genes** involved in *M. tuberculosis* resistance mechanisms  
- Total # mutations (WHO 2021): **4,709**
- Category 3 (“Uncertain”) mutations dominate the dataset and form the **test set** for model inference.



##  Features Used

Each mutation is represented by a set of biologically and structurally meaningful features:

| Feature | Description |
|----------|--------------|
| **Proximity (3D)** | Structural distance to nearest resistance-conferring mutation |
| **Thermostability (Rosetta)** | Predicted change in protein stability upon mutation |
| **Log-Likelihood Ratio (LLR)** | Change in per-residue log probability from ESM2 |
| **AAIndex (PCA Distance)** | Biochemical/physicochemical distance between wildtype and mutant amino acids |

These features are derived from multi-source biological data including Rosetta simulations, WHO catalogs, and pretrained ESM2 embeddings.

##  Repository Structure

```

RESISTANCE_FORECAST/
│
├── data/                      # Raw and processed datasets (WHO, embeddings, Rosetta outputs)
├── results/                   # Model outputs, predictions, plots
├── Comparison_Model/           # Stat models 
├── legacy_files/               # Previous versions of notebooks/scripts
│
├── forecast_data_preparation_combined.ipynb   # End-to-end data preprocessing and feature computation
├── model_essential_nonessential_combined.ipynb # Random Forest model for essential/nonessential genes
├── forecast_utils.py           # Shared utility functions (feature merging, sampling, scaling, etc.)
├── mycobrowser.py              # Module for fetching sequence and gene metadata from Mycobrowser
│
├── requirements.txt
├── README.md
└── .gitignore

````



##  Key Workflows

### 1. **Data Preparation**
File: [`forecast_data_preparation_combined.ipynb`](forecast_data_preparation_combined.ipynb)

- Parses WHO 2021 and WHO 2023 mutation catalog  
- Converts three-letter to one-letter amino acid codes  
- Generates mutated protein sequences  
- Computes:
  - Proximity to R-conferring sites (from PDB)
  - LLR (ESM2 log-likelihood)
  - AAIndex distances  
  - Rosetta-based thermostability  
- Outputs integrated feature tables for modeling.

### 2. **Modeling and Prediction**
File: [`model_essential_nonessential_combined.ipynb`](model_essential_nonessential_combined.ipynb)

- Trains **classifier models** separately for essential and nonessential genes and on combined genes.   
- Evaluates via AUC, sensitivity, specificity
- Predicts confidence categories for Category 3 (“Uncertain”) mutations  
- Supports **weighted sampling** to handle class imbalance.

### 3. **Utility Scripts**
- `forecast_utils.py`: helper functions for feature merging, resampling, normalization, and evaluation  
- `mycobrowser.py`: retrieves reference gene/protein info for sequence-based computations


##  Data Access

| Type | Description | Link |
|------|--------------|------|
|  **Training Data** | Preprocessed feature tables (ΔZ, LLR, Rosetta, AAIndex, etc.) used for model training | [Download ZIP](https://umass.sharepoint.com/:u:/r/sites/AntibioticResistanceForcast/Shared%20Documents/resistance_forecast_training_data.zip?csf=1&web=1&e=hA8E8i) |
|  **Prediction Results** | Model outputs and predicted confidence categories for Category 3 mutations | [Download ZIP](https://umass.sharepoint.com/:u:/r/sites/AntibioticResistanceForcast/Shared%20Documents/resistance_forecast_predictions_csv.zip?csf=1&web=1&e=xLaya8) |

Each ZIP file contains per-gene CSVs and combined summary tables compatible with the notebooks in this repository.



## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/resistance_forecast.git
cd resistance_forecast
````

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Prepare data
Run the following notebook to generate all training and evaluation datasets:

```
forecast_data_preparation_combined.ipynb
```

### 4. Train and evaluate models
Then run:
```
model_essential_nonessential_combined.ipynb
```
This notebook trains classifier models for feature sets and performs evaluation and prediction.


##  Evaluation Metrics

* **AUC** and **ROC Curves**
* **Sensitivity / Specificity**
* **Confusion Matrix**
* **Feature Importances**
* **Violin plots of feature dist**
* **SHAP-based Interpretability** [in progress]




##  Requirements

* Python 3.8+
* Key dependencies:

  ```bash
  pandas
  numpy
  torch
  esm
  scikit-learn
  biopython
  matplotlib
  seaborn
  ```

Install all at once:

```bash
pip install -r requirements.txt
```

---

##  Contact

For questions or collaborations, please reach out to:
 **[mtasmin@umass.edu](mailto:mtasmin@umass.edu)**

---

##  Citation

If you use this repository or dataset, please cite:

> Tasmin, M. et al. *Resistance Forecast: Structure-Informed Prediction of Antibiotic Resistance in Mycobacterium tuberculosis* (2025, in prep.)





