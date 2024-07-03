
# Protein Sequence Mutation and Delta-Z Calculation

This repository contains scripts and data for generating all possible mutations from wildtype protein sequences and calculating the delta-Z values using ESM-2 embeddings. The data and scripts are based on a study of Mycobacterium tuberculosis genes.

## Table of Contents
- [Background](#background)
- [Data](#data)
- [Scripts](#scripts)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

## Background
This project aims to identify and analyze the effects of mutations on the protein sequences of specific genes in Mycobacterium tuberculosis. Using the ESM-2 model, embeddings are generated for both wildtype and mutated sequences, and the delta-Z values are calculated to quantify the changes.

## Data
### Input Files
1. `protein_sequences.csv`: Contains wildtype protein sequences for genes of interest.
2. `mutations_with_one_letter.csv`: WHO catalog with mutation information.
3. `delta_z_values.csv`: Delta-Z values calculated for each mutation.

### Output Files
1. `mutated_sequences/`: Directory containing FASTA files with wildtype and mutated sequences.
2. `final_delta_z_with_catalog_info.csv`: CSV file with delta-Z values and additional information from the WHO catalog.

## Scripts
1. `generate_mutations.py`: Generates mutated sequences from wildtype protein sequences using WHO catalog information.
2. `delta_z_calculation.py`: Calculates delta-Z values for wildtype and mutated sequences using ESM-2 embeddings.
3. `merge_catalog_data.py`: Merges delta-Z values with WHO catalog information.

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/Protein-Sequence-Mutation-Delta-Z.git
   cd Protein-Sequence-Mutation-Delta-Z
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Convert amino acid representation (3-to-1):
   ```bash
   python aa_change_3_1.py
   ```

4. Generate mutated sequences:
   ```bash
   python generate_mutations.py
   ```

4. Calculate delta-Z values:
   ```bash
   python delta_z_calculation.py
   ```

5. Merge delta-Z values with catalog information:
   ```bash
   python merge_catalog_data.py
   ```

## Results
The `/data/final_delta_z_with_catalog_info.csv` file contains the delta-Z values along with the corresponding genome position, nucleotide change, and amino acid change for each mutation.

## Contributing
Contributions are welcome! Please fork this repository and submit pull requests for any enhancements or bug fixes.

