import csv
import pandas as pd
import os

# Read the CSV files
protein_sequences_df = pd.read_csv('/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/protein_sequences.csv')
mutations_df = pd.read_csv('/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/mutations_with_one_letter_all.csv')

# Output directory
output_dir = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/mutated_sequences'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to apply a mutation
def apply_mutation(wildtype_seq, mutation):
    mutation_info = mutation.split('.')
    if len(mutation_info) != 2:
        return None
    mutation_type, change = mutation_info
    if mutation_type != 'p':
        return None
    original_aa = change[0]
    position = int(''.join(filter(str.isdigit, change)))
    new_aa = change[-1]
    if position <= 0 or position > len(wildtype_seq):
        return None
    mutated_seq = wildtype_seq[:position - 1] + new_aa + wildtype_seq[position:]
    return original_aa, new_aa, position, mutated_seq

# Process each gene
for index, row in protein_sequences_df.iterrows():
    gene = row['Gene']
    rv_id = row['RV']
    wildtype_seq = row['Protein Sequence']
    
    # Output FASTA file for the gene
    output_file = os.path.join(output_dir, f"{rv_id}_{gene}.fasta")
    
    with open(output_file, mode='w') as file:
        # Write wildtype sequence
        file.write(f">{rv_id}|{gene}|Wildtype\n")
        file.write(f"{wildtype_seq}\n")
        
        # Get mutations for the specific gene
        gene_mutations = mutations_df[mutations_df['gene'] == gene]
        
        for _, mut_row in gene_mutations.iterrows():
            mutation = mut_row['one_letter_mutation']
            result = apply_mutation(wildtype_seq, mutation)
            if result:
                original_aa, new_aa, position, mutated_seq = result
                mutation_label = f"{gene}_p.{original_aa}{position}{new_aa}"
                file.write(f">{rv_id}|{gene}|{mutation_label}\n")
                file.write(f"{mutated_seq}\n")
                print(f"Mutation {mutation} applied for gene {gene}.")

print(f"Mutated sequences have been saved to {output_dir}.")
