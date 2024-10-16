import pandas as pd

mutations_file = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/mutations_with_one_letter_all.csv'
delta_z_file = '/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/delta_z_values.csv'


mutations_df = pd.read_csv(mutations_file)
delta_z_df = pd.read_csv(delta_z_file)


# Ensure the 'Mutation' column is a string and handle missing values
delta_z_df['Mutation'] = delta_z_df['Mutation'].astype(str).fillna('')

# Remove the gene name from the beginning of the 'Mutation' column
delta_z_df['Mutation'] = delta_z_df['Mutation'].apply(lambda x: x.split('_')[1] if '_' in x else x)

# Filter out entries with NaN mutations before merging
filtered_delta_z_df = delta_z_df[~delta_z_df['Mutation'].isna()]


# Merge the filtered Delta-Z DataFrame with the filtered mutations DataFrame to get the additional columns
merged_df = pd.merge(filtered_delta_z_df, mutations_df, left_on=['Gene', 'Mutation'], right_on=['gene', 'one_letter_mutation'], how='left')
# Save the final merged DataFrame to a CSV file
merged_df.to_csv('/work/pi_annagreen_umass_edu/mahbuba/resistance_forecast/data/final_delta_z_with_catalog_info.csv', index=False)

print(f"Updated Delta-Z values with additional columns have been saved.")
