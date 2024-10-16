import pandas as pd
import torch
import numpy as np
from transformers import AutoTokenizer, EsmForMaskedLM
import re

# Function to load protein sequence based on the gene name
def load_protein_sequence(gene_name, protein_seq_file):
    protein_df = pd.read_csv(protein_seq_file)
    protein_sequence = protein_df.loc[protein_df['gene'] == gene_name, 'protein_sequence'].values[0]
    return protein_sequence

# Function to load phenotype data based on the gene name
def load_phenotype_data(gene_name, base_path):
    # phenotype_data_path = f'{base_path}/{gene_name}_phenotype_data.csv'
    phenotype_data_path = f'{base_path}/{gene_name}_aaindex_data.csv'
    phenotype_data = pd.read_csv(phenotype_data_path)
    return phenotype_data, phenotype_data_path

# Function to extract mutation information
def extract_mutation_info(mutation):
    match = re.match(r"([A-Z])(\d+)([A-Z])", mutation)
    if match:
        wt_residue = match.group(1)
        position = int(match.group(2)) - 1  # Convert to 0-based index
        mt_residue = match.group(3)
        return wt_residue, position, mt_residue
    else:
        return None, None, None

# Function to calculate LLR for a single mutation
def calculate_llr_for_mutation(wt_residue, mt_residue, position, input_ids, model, tokenizer):
    # Mask the target position
    masked_input_ids = input_ids.clone()
    masked_input_ids[0, position + 1] = tokenizer.mask_token_id  # +1 for special token offset

    # Get logits for the masked token
    with torch.no_grad():
        logits = model(masked_input_ids).logits

    # Calculate log probabilities
    probabilities = torch.nn.functional.softmax(logits[0, position + 1], dim=0)  # +1 for special token offset
    log_probabilities = torch.log(probabilities)

    # Get the log probability of the wild-type residue
    log_prob_wt = log_probabilities[tokenizer.convert_tokens_to_ids(wt_residue)].item()

    # Calculate LLR for the mutant residue
    log_prob_mt = log_probabilities[tokenizer.convert_tokens_to_ids(mt_residue)].item()
    llr_score = log_prob_mt - log_prob_wt

    return llr_score

# Function to calculate LLR scores for unique mutations
def calculate_llr_scores(phenotype_data, protein_sequence, model, tokenizer):
    input_ids = tokenizer.encode(protein_sequence, return_tensors="pt")
    unique_mutations = phenotype_data[['wt_residue', 'mutation_position', 'mutant_residue']].drop_duplicates()

    llr_scores = {}
    for _, row in unique_mutations.iterrows():
        wt_residue = row['wt_residue']
        position = row['mutation_position']
        mt_residue = row['mutant_residue']
        llr_score = calculate_llr_for_mutation(wt_residue, mt_residue, position, input_ids, model, tokenizer)
        llr_scores[(wt_residue, position, mt_residue)] = llr_score

    # Map the LLR scores back to the phenotype data
    phenotype_data['LLR'] = phenotype_data.apply(
        lambda x: llr_scores[(x['wt_residue'], x['mutation_position'], x['mutant_residue'])], axis=1
    )

    return phenotype_data

# Function to save the updated phenotype data
def save_phenotype_data(phenotype_data, output_path):
    phenotype_data.to_csv(output_path, index=False)

# Main function to run the analysis
def main(gene_name, protein_seq_file, phenotype_base_path):
    # Load protein sequence
    protein_sequence = load_protein_sequence(gene_name, protein_seq_file)

    # Load phenotype data
    phenotype_data, phenotype_data_path = load_phenotype_data(gene_name, phenotype_base_path)

    # Add mutation info columns
    phenotype_data['wt_residue'], phenotype_data['mutation_position'], phenotype_data['mutant_residue'] = zip(*phenotype_data['one_letter_mutation'].apply(extract_mutation_info))

    # Load the ESM-2 model and tokenizer
    model_name = "facebook/esm2_t30_150M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name)

    # Calculate LLR scores for unique mutations
    phenotype_data = calculate_llr_scores(phenotype_data, protein_sequence, model, tokenizer)
    
    output_path = f'/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/{gene_name}_phenotype_data.csv'
    # Save the updated phenotype data
    save_phenotype_data(phenotype_data, phenotype_data_path)

# Run the analysis
if __name__ == "__main__":
    gene_name = 'pncA'  # Change this as needed
    protein_seq_file = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/protein_sequences.csv'  # Update path
    phenotype_base_path = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data'

    main(gene_name, protein_seq_file, phenotype_base_path)
