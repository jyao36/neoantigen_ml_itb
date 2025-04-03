## Test 

## Merge immuno pipeline output files for a given patient to prepare for prediction with the ML model 

## Description: 
# The script merges the following files (using 10146-CO070-0054 as an example):
# 1. class1 10146-0054-TumorDNA.all_epitopes.aggregated.tsv
# 2. class1 10146-0054-TumorRNA.all_epitopes.tsv
# 3. class2 10146-0054-TumorRNA.all_epitopes.aggregated.tsv

# %% Import dependencies
from pathlib import Path
import pandas as pd

# %% Set up directories
patient_id = "10146-0052" # the patient ID used for the folder name 
pvacseq_results_dir = f"/Volumes/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/{patient_id}/gcp_immuno/final_results/pVACseq"
print(pvacseq_results_dir)


# %% Find and read files
# Define a reusable function to find and read files
def find_and_read_file(directory, pattern, file_description):
    """
    Finds the first file matching the given pattern in the specified directory
    and reads it into a pandas DataFrame.

    Args:
        directory (Path): The directory to search in.
        pattern (str): The glob pattern to match files.
        file_description (str): A description of the file for logging purposes.

    Returns:
        pd.DataFrame: The loaded DataFrame, or None if no file is found.
    """
    files = list(directory.glob(pattern))
    if files:
        file_path = files[0]  # Take the first matching file
        print(f"Found {file_description}: {file_path}")
        df = pd.read_csv(file_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False) # keep “None" values as "None" which we will change to 0 later. Others will be treated as missing values. 
        print(f"Loaded {file_description} into DataFrame with shape: {df.shape}")
        return df
    else:
        print(f"No {file_description} found matching pattern '{pattern}' in {directory}")
        return None


# Class 1 aggregated file
mhc1_dir = Path(f"{pvacseq_results_dir}/mhc_i/")
mhc1_agg_df = find_and_read_file(mhc1_dir, "*.all_epitopes.aggregated.tsv", "class 1 aggregated file")

# Class 1 all epitopes file
mhc1_allepi_df = find_and_read_file(mhc1_dir, "*.all_epitopes.tsv", "class 1 all epitopes file")

# Class 2 aggregated file
mhc2_dir = Path(f"{pvacseq_results_dir}/mhc_ii/")
mhc2_agg_df = find_and_read_file(mhc2_dir, "*.all_epitopes.aggregated.tsv", "class 2 aggregated file")



# %% Reformat column names 

# Rename columns in mhc1_agg_df
mhc1_agg_df.rename(columns={
    "Best Peptide": "Best Peptide class1", 
    "Best Transcript": "Best Transcript class1", 
    "IC50 MT": "IC50 MT class1", 
    "IC50 WT": "IC50 WT class1", 
    "%ile MT": "%ile MT class1", 
    "%ile WT": "%ile WT class1"
    }, inplace=True)


# Rename columns in mhc2_agg_df
mhc2_agg_df.rename(columns={
    "Best Peptide": "Best Peptide class2", 
    "Best Transcript": "Best Transcript class2", 
    "IC50 MT": "IC50 MT class2", 
    "IC50 WT": "IC50 WT class2", 
    "%ile MT": "%ile MT class2", 
    "%ile WT": "%ile WT class2"
    }, inplace=True)

columns_to_keep = [
    "ID",  # Always include "ID"
    *[col for col in mhc2_agg_df.columns if col.startswith("D") and "DNA VAF" not in col],  # Columns starting with "D" but without "DNA VAF"
    *[col for col in mhc2_agg_df.columns if col.endswith("class2")]  # Renamed columns ending with "class2"
]
# Select the subset of columns
mhc2_agg_subset = mhc2_agg_df[columns_to_keep]

# Compare the number of rows
if mhc1_agg_df.shape[0] != mhc2_agg_subset.shape[0]:
    print("Warning: Class 1 aggregated file DOES NOT have the same number of rows as Class 2 aggregated file.")




# %% Merge dataframes
# Merge class 1 and class 2 aggregated dataframes
merged_df = pd.merge(mhc1_agg_df, mhc2_agg_subset, on="ID", how="left")
# Merge with class 1 all epitopes dataframe


# Split the "ID" column in merged_df into separate columns for matching to all epitopes file
merged_df_split = merged_df.copy()
merged_df_split[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']] = merged_df_split['ID'].str.split('-', expand=True)

# Convert "Start" and "Stop" columns to integers
merged_df_split['Start'] = merged_df_split['Start'].astype(int)
merged_df_split['Stop'] = merged_df_split['Stop'].astype(int)

# Perform an inner join with mhc1_allepi_df on the specified columns
merged_all = pd.merge(
    merged_df_split,
    mhc1_allepi_df,
    how='inner',
    left_on=[
        'Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 
        'Best Transcript class1', 'Best Peptide class1', 'Allele'
    ],
    right_on=[
        'Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 
        'Transcript', 'MT Epitope Seq', 'HLA Allele'
    ]
)

# Drop the redundant columns created from splitting "ID"
merged_all.drop(columns=['Chromosome', 'Start', 'Stop', 'Reference', 'Variant'], inplace=True)

# Print the resulting DataFrame
print("Merged DataFrame:")
print(merged_all.head())
# %%
