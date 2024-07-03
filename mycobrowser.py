import requests
from bs4 import BeautifulSoup
import csv

def get_protein_sequence(gene_id):
    url = f"https://mycobrowser.epfl.ch/genes/{gene_id}"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Failed to retrieve data for {gene_id}")
        return None
    
    soup = BeautifulSoup(response.text, 'html.parser')
    # Locate the specific div containing the protein sequence
    panel_divs = soup.find_all('div', class_='panel panel-default')
    
    for panel_div in panel_divs:
        heading = panel_div.find('div', class_='panel-heading')
        if heading and 'Protein sequence' in heading.text:
            pre_tag = panel_div.find('pre')
            if pre_tag:
                sequence = pre_tag.text.strip().split('\n', 1)[1]  # Remove the first line
                return sequence
    
    print(f"Protein sequence not found for {gene_id}")
    return None

# List of gene IDs
gene_ids = [
    'Rv1305',
    'Rv3919c',
    'Rv0682',
    'Rv0701',
    'Rv0005',
    'Rv0667',
    'Rv2043c',
    'Rv1694',
    'Rv1484',
    'Rv0407',
    'Rv1908c',
    'Rv3795',
    'Rv0678',
    'Rv2535c',
    'Rv3547',
    'Rv0006',
    'Rv3854c'
]

# Retrieve protein sequences and save to a CSV file
with open('protein_sequences.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Gene ID', 'Protein Sequence'])

    for gene_id in gene_ids:
        sequence = get_protein_sequence(gene_id)
        if sequence:
            writer.writerow([gene_id, sequence])
            print(f"Protein sequence for {gene_id} saved.")
        else:
            writer.writerow([gene_id, 'Not found'])
            print(f"Protein sequence for {gene_id} not found.")
