import pandas as pd
import re, sys

input_file = sys.argv[1] # "variant.csv"
output_file = sys.argv[2] # "new_file.csv"
offset = int(sys.argv[3]) # The offset between the Rosetta and Original residue index
offset_residue_num = int(sys.argv[4]) # The residue number from which the offset will start

# input the variant file obtained from ddG_calculation.py
data = pd.read_csv(input_file)
variant = data["variant"]

# Funtion to change the numbering from Rosetta index to the origianl PDB index
def add_to_number_in_variant(variant, offset, offset_start):
    """
    Convert residue numbering from Rosetta to the original PDB indexing.

    Input Parameters:
        variants (str): currnet variant in the list
        offset (int): The offset between the Rosetta and Original residue index    
        offset_start (int): The residue number from which the offset will start

    Output Parameters:
        new_variant (str): Original residue index in the current varaint
    """

    # It divides the variant into wild-type residue alphabet, Rosetta residue index, mutated residue alphabet
    match = re.match(r"([A-Z]+)(\d+)([A-Z]+)", variant)
    prefix, number, suffix = match.groups()

    # The offset between the Rosetta and Original residue index
    if int(number) >= offset_start:
        new_number = int(number) + offset

    return f"{prefix}{new_number}{suffix}"

new_variants = [add_to_number_in_variant(var, offset, offset_residue_num) for var in variant]
data["variant"] = new_variants
print(data["variant"]) 

# new datafram with correct residue index
data.to_csv(output_file, index=False)
