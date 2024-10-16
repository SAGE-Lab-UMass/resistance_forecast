

# Resistance Forecast in *Mycobacterium tuberculosis* Antibiotics

This project focuses on predicting the confidence category of unlabeled mutations in *Mycobacterium tuberculosis* (M. tuberculosis), which are associated with antibiotic resistance. By leveraging features derived from biological and structural data, such as proximity to known resistance-conferring mutations, changes in thermostability, and log-likelihood scores from protein language models, we aim to classify mutations into categories that indicate their likelihood of being associated with drug resistance.

## Objective

The primary task is to predict the confidence category for mutations in *M. tuberculosis* genes that are currently labeled as having "uncertain significance." These categories are based on evidence regarding the mutations’ association with antibiotic resistance and are classified into five levels:

- **1) Assoc w R**: Strongly associated with resistance.
- **2) Assoc w R - Interim**: Interim evidence suggests association with resistance.
- **3) Uncertain significance**: This category is used for test data with no clear evidence.
- **4) Not assoc w R - Interim**: Interim evidence suggests no association with resistance.
- **5) Not assoc w R**: Strong evidence suggests no association with resistance.

## Dataset Overview

This project includes mutations from 17 essential proteins involved in *M. tuberculosis* antibiotic resistance. These proteins are critical to understanding the mechanisms of drug resistance. The total number of labeled data points is **6172**, and the target phenotype distribution is as follows:

- **5699** mutations in category 3 (Uncertain significance)
- **225** mutations in category 2 (Assoc w R - Interim)
- **186** mutations in category 1 (Assoc w R)
- **41** mutations in category 5 (Not assoc w R)
- **21** mutations in category 4 (Not assoc w R - Interim)

## Features

To predict the confidence category of mutations, several important features were extracted from biological data. These features are key in distinguishing resistant mutations from non-resistant ones:

1. **Frequency in Population**: The proportion of resistant to susceptible variants calculated as:
   \[
   \text{Frequency} = \frac{\text{Present\_R}}{\text{Present\_R} + \text{Present\_S}}
   \]
   
2. **Proximity to Nearest Non-Self R-Conferring Mutation**: This feature represents the structural proximity of the mutation to known resistance-conferring mutations, derived from protein structure data (PDB). The distance is expressed as a floating-point value.

3. **Delta-Z in Embedding Model (ESM2)**: The Euclidean distance between the embedding of the wildtype and mutated sequences, using a pretrained protein language model (ESM2). Larger values indicate significant changes in the protein's functional or structural properties.

4. **Change in Thermostability**: Changes in the protein's stability due to mutation, calculated using Rosetta modeling. This floating-point value represents how a mutation impacts the protein’s structural stability.

5. **Log Likelihood Score (ESM2)**: A measure of how likely the mutation is to be conserved based on evolutionary alignment, calculated from the log-likelihood score derived from ESM2 models. Higher scores suggest a greater likelihood that the mutation is evolutionarily favored.

6. **Amino Acid Index (AA PCA)**: The Euclidean distance between the wildtype and mutant amino acids based on physicochemical properties (PCA of amino acid indices). This feature captures the biochemical impact of the mutation.

These features provide a comprehensive view of how each mutation may affect the protein's structure and function, which in turn informs the model's prediction of the mutation’s role in antibiotic resistance.

# Protein Mutation Data Preparation

This repository contains scripts for preparing mutation-related data, calculating feature values, and generating derived datasets for analyzing antibiotic resistance in *Mycobacterium tuberculosis*. The primary focus is to convert mutation information, process protein sequences, calculate multiple biological features, and output relevant data for further analysis. The data preparation steps for this project can be found in the [forecast_data_preparation.ipynb](forecast_data_preparation.ipynb) file. This notebook walks through the steps of generating the mutation data, processing protein sequences, and computing the features used in the prediction model.

<details>
  <summary>Click to expand for detailed data preparation steps</summary>

## Folder Structure

- **data/**: Contains various input datasets, including protein sequences, WHO catalog data, and derived features.
- **mutated_sequences/**: Generated protein sequences with mutations applied.

## Prerequisites

- Python 3.8+
- Libraries: `pandas`, `numpy`, `re`, `csv`, `os`, `torch`, `esm`, `Bio`, `evcouplings`, `transformers`, `matplotlib`.

## Steps and Functions

### 1. **Three-to-One Amino Acid Conversion**
   This function converts three-letter amino acid mutation codes into one-letter codes using regular expressions.

   ```python
   def convert_three_to_one(mutation)
   ```

   Example:
   ```python
   convert_three_to_one("atpE_p.Asp28Gly")  # Output: p.D28G
   ```

### 2. **WHO Catalog Parsing**
   - Loads the WHO catalog (`HO-UCN-TB-2023.7-eng.xlsx`) and extracts relevant columns like `gene`, `mutation`, `drug`, and `confidence level`.
   - Filters target genes and applies the mutation conversion.
   - Saves the processed data to `mutations_with_one_letter_all_confidence.csv`.

### 3. **Frequency Calculation**
   - Computes the frequency of each variant in the population as:
     \[
     \text{Frequency} = \frac{\text{Present\_R}}{\text{Present\_R} + \text{Present\_S}}
     \]
   - Output stored in `all_proteins_frequency_catalog.csv`.

### 4. **Mutated Sequences Generation**
   - For each protein, applies mutations based on the WHO catalog and generates the mutated protein sequences in FASTA format.
   - Output directory: `mutated_sequences/`.

### 5. **Delta-Z Value Calculation**
   - Uses embeddings from the pretrained ESM-2 model (`esm2_t6_8M_UR50D`) to compute the delta-z value, which measures the Euclidean distance between wildtype and mutated sequence embeddings.
   - Resulting values are saved to `all_delta_z_values.csv`.

### 6. **Proximity to Nearest Resistance-Conferring Mutations**
   - For each gene, calculates the proximity to the nearest resistance-conferring mutations using protein structure distance maps.
   - Output: `all_proteins_freq_details_proximity.csv`.

### 7. **Amino Acid Index Distance**
   - Computes Euclidean distances between wildtype and mutant amino acids based on their physicochemical properties (amino acid indices).
   - Output stored in `all_proteins_freq_details_proximity_aaindex.csv`.

### 8. **Log-Likelihood Ratio (LLR) Calculation**
   - Uses the ESM-2 model to calculate the log-likelihood ratio for each mutation by comparing the wildtype and mutated amino acid probabilities.
   - Results stored in `all_proteins_freq_details_proximity_aaindex_llr.csv`.

### 9. **Thermostability Score Calculation**
   - Extracts thermostability scores for specific genes from Rosetta modeling and merges them with the mutation dataset.
   - Output stored as `*_thermostability.csv` for each gene.

## Output Files

- **mutations_with_one_letter_all_confidence.csv**: Processed mutation data with 3-to-1 amino acid conversion.
- **all_proteins_frequency_catalog.csv**: Data with frequency calculations for variants.
- **all_delta_z_values.csv**: Delta-Z values for mutation embeddings.
- **all_proteins_freq_details_proximity.csv**: Data with proximity to nearest resistance-conferring mutations.
- **all_proteins_freq_details_proximity_aaindex.csv**: Mutation data enriched with amino acid index distances.
- **all_proteins_freq_details_proximity_aaindex_llr.csv**: Mutation data with LLR scores.
- **{gene_name}_thermostability.csv**: Gene-specific data with thermostability scores.
</details>



# Random Forest Classification on Mutation Data

This repository contains a Random Forest classifier implementation to predict antibiotic resistance phenotypes in *Mycobacterium tuberculosis* based on various mutation-related features. The script processes previously prepared datasets, trains a model, and evaluates its performance using several metrics, including accuracy, feature importance, confusion matrix, and ROC curve. The implemented code can be found in the [random_forest.ipynb](random_forest.ipynb) file. 
<details>
  <summary>Click to expand for detailed training and classification steps. </summary>

## Prerequisites

- Python 3.8+
- Libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`, `torch`, `transformers`.

## Data Preparation

The input data for this Random Forest classifier is prepared from previous steps, where various biological features were extracted. These features include:
1. **Delta-Z**: Euclidean distance between wildtype and mutated protein sequence embeddings.
2. **Frequency**: Proportion of resistant vs susceptible variants.
3. **Proximity to Resistance-Conferring Mutations**: Structural distance to known resistance mutations.
4. **Amino Acid Index Distance**: Physicochemical property distance between wildtype and mutated amino acids.
5. **Log-Likelihood Ratio (LLR)**: Change in log-probability between wildtype and mutated residues.
6. **Thermostability**: Rosetta-based predicted protein stability changes due to mutations.

## Random Forest Model

### Steps

1. **Data Loading**: Load the phenotype data from `all_proteins_freq_details_proximity_aaindex_llr.csv` or similar files.
   
   ```python
   catalog_data = pd.read_csv('/path/to/your/data.csv')
   ```

2. **Data Cleaning**: Remove duplicates and handle missing values, especially in the numerical columns. The target variable (`confidence` level) is mapped to numerical labels.

3. **Feature Selection**: The following features are used for training the model:
   - `delta_z`
   - `frequency`
   - `Proximity_to_R_Conferring`
   - `aa_index_dist`
   - `llr_score`
   - `thermostability`

4. **Data Scaling**: Features are scaled using `StandardScaler`.

5. **Train-Test Split**: The dataset is split into training and testing sets using an 70-30% split.

6. **Random Forest Training**: A `RandomForestClassifier` with 100 trees is trained on the dataset.

   ```python
   rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
   rf_classifier.fit(X_train_split, y_train_split)
   ```

7. **Evaluation Metrics**:
   - **Confusion Matrix**: To visualize the true vs predicted classifications.
   - **Classification Report**: Provides precision, recall, and F1-score for each class.
   - **Accuracy Score**: Overall accuracy of the classifier.
   - **ROC Curve**: Receiver Operating Characteristic curves are plotted for each class to assess the model's ability to distinguish between classes.

8. **Feature Importance**: The model's feature importance is visualized using a bar chart.

### Output Files and Visualizations

- **Confusion Matrix**: Displays the accuracy of predictions for each confidence category.
- **Classification Report**: Shows precision, recall, F1-score, and support for each category.
- **ROC Curve**: Plots the trade-off between true positive rate and false positive rate for each class.
- **Feature Importance**: Bar plot showing the relative importance of each feature in the Random Forest model.

### Code Structure

```python
# Data loading and cleaning
catalog_data = pd.read_csv('/path/to/your/data.csv')

# Feature selection and scaling
X_train = train_data[numeric_columns]
y_train = train_data['confidence']

# Train-test split
X_train_split, X_test_split, y_train_split, y_test_split = train_test_split(X_scaled, y_train, test_size=0.3, random_state=42)

# Train the Random Forest classifier
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
rf_classifier.fit(X_train_split, y_train_split)

# Model evaluation
print(confusion_matrix(y_test_split, y_pred_split))
print(classification_report(y_test_split, y_pred_split))
```

## Category 3 Prediction

The model also predicts entries in the dataset labeled as "3) Uncertain significance" by applying the trained Random Forest model to this subset.

```python
# Predict the confidence categories for the entries labeled as category 3
category_3_data['predicted_confidence'] = rf_classifier.predict(X_category_3)
```

## Weighted Sampling

To address class imbalance, weighted sampling is performed, ensuring that under-represented classes have a proportionally larger impact on the model's learning process.

```python
# Perform weighted sampling using pandas
resampled_df = filtered_data.sample(n=len(filtered_data),
                                    replace=True,
                                    weights=filtered_data['weights'],
                                    random_state=42)
```

## Additional Features

- **Thermostability Feature Integration**: The model can include a sixth feature, `thermostability`, derived from Rosetta-based protein stability scores.
</details>


## Contact

For questions or feedback, please reach out to the mtasmin@umass.edu

