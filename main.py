import pandas as pd

# Load the metadata file (for example: /yourpath/metadata_wargo_melanoma.txt)
df = pd.read_csv(METADATA_FILE_PATH, sep='\t', header=None)

# 1. Select rows where the second column contains the word 'age'
age_filter = df[df[2].str.contains('age')]

# 2. Split the values in the second column and create a new DataFrame with each value in separate columns
split_values = [x.split(";") for x in age_filter[2]]
split_df = pd.DataFrame(split_values)

# 3. Filter out rows where the 10th column contains null values
filtered_df = split_df[split_df[10].notna()]

# 4. Further filter the DataFrame to exclude rows where the 10th column contains the word 'Wargo'
filtered_df_no_wargo = filtered_df[~filtered_df[10].str.contains('Wargo')]

# 5. Create a new 'ID' column by extracting IDs from the 10th column using a list comprehension
filtered_df_no_wargo['ID'] = [x.split('=')[1] for x in filtered_df_no_wargo[10]]

# 6. Load the report file and extract 'sample_alias' values
report_df = pd.read_csv(REPORT_FILE_PATH, sep='\t')
sample_aliases = [str(alias) for alias in report_df['sample_alias'].values]

# 7. Select rows in the first DataFrame where 'ID' matches any of the 'sample_alias' values
matched_rows = filtered_df_no_wargo[filtered_df_no_wargo['ID'].isin(sample_aliases)]

# 8. Rename the columns of the matched DataFrame for clarity
new_column_names = ['age', 'idh', 'phenotype', 'type', 'pi', 'pt', 'prog', 'race', 'sex', 'stage', 'id', 'time', 'treat', 'treat_sub', 'cd8', 'ena', '-', 'ID']
matched_rows.columns = new_column_names

# 9. Extract 'R' and 'NR' values from the 'phenotype' column and create a new column 'R_NR'
matched_rows['R_NR'] = [x.split('=')[1] for x in matched_rows['phenotype']]

# 10. Select columns 'run_accession', 'sample_accession', and 'sample_alias' from the report DataFrame and create a new DataFrame
new_report_df = pd.DataFrame(report_df[['run_accession', 'sample_accession', 'sample_alias']].values)

# 11. Rename the columns of the new DataFrame to ['variable', 'sample_accession', 'ID']
new_column_names_report_df = ['variable', 'sample_accession', 'ID']
new_report_df.columns = new_column_names_report_df

# 12. Merge the metadata DataFrame with the new report DataFrame based on the 'ID' column
new_report_df["ID"] = new_report_df["ID"].astype(str)
merged_df = new_report_df.merge(matched_rows, on='ID')

# 13. Rename columns of the merged DataFrame for clarity
final_column_names = ['variable', 'sample_accession', 'ID', 'age', '1', 'phenotype', '3', '4', '5', '6', '7', 'sex', '9', '10', '11', '12', '13', '14', '15', '16']
merged_df.columns = final_column_names

# 14. Extract 'R' and 'NR' values from the 'phenotype' column, 'sex', and 'age' columns
merged_df['R_NR'] = [x.split('=')[1] for x in merged_df['phenotype']]
merged_df['sex'] = [x.split('=')[1] for x in merged_df['sex']]
merged_df['age'] = [x.split('=')[1] for x in merged_df['age']]

# 15. Select and keep only the columns 'variable', 'age', 'R_NR', and 'sex'
final_df = merged_df[['variable', 'age', 'R_NR', 'sex']]

# 16. Export the filtered DataFrame to a TSV file
final_df.to_csv(OUTPUT_FILE_PATH, sep='\t', index=False)

# 17. Generate a download link for the TSV file (if using Google Colab or similar environments)
from IPython.display import FileLink
FileLink(OUTPUT_FILE_PATH)

# Alternatively, for Google Colab, you can use:
# from google.colab import files
# files.download(OUTPUT_FILE_PATH)
