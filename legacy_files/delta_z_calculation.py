import os
import torch
import numpy as np
import pandas as pd
import esm
from Bio import SeqIO

# Load the pretrained ESM-2 model
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval().to('cuda')

# Function to calculate embeddings
def get_embedding(sequence):
    sequence = sequence.upper()  # Ensure sequence is in uppercase
    sequence = sequence.replace('*', '')  # Remove any stop codons
    batch_labels, batch_strs, batch_tokens = batch_converter([("sequence", sequence)])
    with torch.no_grad():
        results = model(batch_tokens.to('cuda'), repr_layers=[6])
    token_embeddings = results["representations"][6].cpu().numpy()
    return token_embeddings.mean(axis=1)

# Directory containing the FASTA files
fasta_dir = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/mutated_sequences'

# Output CSV file
output_file = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/delta_z_values.csv'

# Initialize the results list
results = []

# Process each FASTA file
for fasta_file in os.listdir(fasta_dir):
    if fasta_file.endswith('.fasta'):
        file_path = os.path.join(fasta_dir, fasta_file)
        records = list(SeqIO.parse(file_path, "fasta"))
        print("loaded file: ", file_path)
        
        # Extract RV and gene names from the file name
        rv_id, gene_name = fasta_file.replace('.fasta', '').split('_')
        print("continue")
        
        # Find the wildtype sequence
        wildtype_record = None
        for record in records:
            if "Wildtype" in record.description:
                wildtype_record = record
                break
        
        if wildtype_record is None:
            print(f"No wildtype sequence found in {fasta_file}")
            continue
        
        wildtype_seq = str(wildtype_record.seq).upper()
        wildtype_embedding = get_embedding(wildtype_seq)
        
        for record in records:
            if record.description == wildtype_record.description:
                continue
            mutation = record.description
            mutated_seq = str(record.seq).upper().replace('*', '')
            try:
                mutated_embedding = get_embedding(mutated_seq)
                delta_z = np.linalg.norm(mutated_embedding - wildtype_embedding)
                results.append([fasta_file, mutation, delta_z, rv_id, gene_name])
            except KeyError as e:
                print(f"Error processing {mutation} in {fasta_file}: {e}")

# Save the results to a CSV file
results_df = pd.DataFrame(results, columns=['Filename', 'Mutation', 'Delta Z', 'RV', 'Gene'])
results_df.to_csv(output_file, index=False)

print(f"Delta Z values have been saved to {output_file}.")
