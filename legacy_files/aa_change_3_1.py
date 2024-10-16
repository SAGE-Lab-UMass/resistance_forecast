import pandas as pd
import re

# Mapping from three-letter to one-letter amino acid codes
three_to_one = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H',
    'Ile': 'I', 'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q',
    'Arg': 'R', 'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
}

def convert_three_to_one(mutation):
    match = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", mutation)
    if not match:
        return None
    original_aa, position, new_aa = match.groups()
    try:
        return f"p.{three_to_one[original_aa]}{position}{three_to_one[new_aa]}"
    except KeyError:
        return None



# Load the mutations CSV
mutations_df = pd.read_csv('/work/pi_annagreen_umass_edu/mahbuba/phenotype_prediction/data/updated_high_confidence_catalog.csv')

# Apply the conversion
mutations_df['one_letter_mutation'] = mutations_df['mutation'].apply(convert_three_to_one)

# Filter out rows with invalid mutations
mutations_df = mutations_df.dropna(subset=['one_letter_mutation'])

# Save the updated DataFrame
mutations_df.to_csv('mutations_with_one_letter_all.csv', index=False)

print(mutations_df[['mutation', 'one_letter_mutation']])
