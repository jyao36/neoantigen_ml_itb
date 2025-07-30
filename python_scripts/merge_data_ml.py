# --------------------------------------------------------------------------------------------------
# Script: merge_data_ml.py
# Description: Merge immuno pipeline output files for a given patient and prepare for ML prediction.
# --------------------------------------------------------------------------------------------------

## Test 

## Merge immuno pipeline output files for a given patient to prepare for prediction with the ML model 

## Description: 
# The script merges the following files (using 10146-CO070-0054 as an example):
# 1. class1 10146-0054-TumorDNA.all_epitopes.aggregated.tsv
# 2. class1 10146-0054-TumorRNA.all_epitopes.tsv
# 3. class2 10146-0054-TumorRNA.all_epitopes.aggregated.tsv

## For testing purposes: 
# remember to mount to the gillandersw server before running the script.

# %% Import dependencies
from pathlib import Path
import pandas as pd
import joblib
import pickle
import numpy as np
import os

# --------------------------------------------------------------------------------------------------
# Section 1: Set up directories and patient information
# --------------------------------------------------------------------------------------------------

patient_id = "10146-0062"  # The patient ID used for the folder name
patient_dir = Path(f"/Volumes/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/{patient_id}")
pvacseq_results_dir = patient_dir / "gcp_immuno" / "final_results" / "pVACseq"
print(f"pVACseq results directory: {pvacseq_results_dir}")


# --------------------------------------------------------------------------------------------------
# Section 2: Define utility function to find and read files
# --------------------------------------------------------------------------------------------------

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
        df = pd.read_csv(file_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False) # keep "None" values as "None" which we will change to 0 later. Others will be treated as missing values. 
        print(f"Loaded {file_description} into DataFrame with shape: {df.shape}")
        return df
    else:
        print(f"No {file_description} found matching pattern '{pattern}' in {directory}")
        return None


# --------------------------------------------------------------------------------------------------
# Section 3: Read input files (class 1/class 2 aggregated and all epitopes)
# --------------------------------------------------------------------------------------------------

mhc1_dir = pvacseq_results_dir / "mhc_i"
mhc1_agg_df = find_and_read_file(mhc1_dir, "*.all_epitopes.aggregated.tsv", "class 1 aggregated file")
mhc1_allepi_df = find_and_read_file(mhc1_dir, "*.all_epitopes.tsv", "class 1 all epitopes file")
# FOR TESTING:------------------------------------------------------------------------------------
# Add a column 'Gene of Interest' filled with NA for testing
if mhc1_allepi_df is not None:
    mhc1_allepi_df['Gene of Interest'] = "False"
#-------------------------------------------------------------------------------------------------

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

# Columns to keep in mhc2_agg_df
columns_to_keep = [
    "ID",  # Always include "ID"
    *[col for col in mhc2_agg_df.columns if col.startswith("D") and "DNA VAF" not in col],  # Columns starting with "D" but without "DNA VAF"
    *[col for col in mhc2_agg_df.columns if col.endswith("class2")]  # Renamed columns ending with "class2"
]
#columns_to_keep = ["ID", "Best Peptide class2", "Best Transcript class2", "IC50 MT class2", "IC50 WT class2", "%ile MT class2", "%ile WT class2"]
#columns_to_keep += [col for col in c2_file.columns if col.startswith("D") and col != "DNA VAF"]
# Select the subset of columns
mhc2_agg_subset = mhc2_agg_df[columns_to_keep]

# Compare the number of rows
if mhc1_agg_df.shape[0] != mhc2_agg_subset.shape[0]:
    print("Warning: Class 1 aggregated file DOES NOT have the same number of rows as Class 2 aggregated file.\n Extra rows in Class 1 aggregated file will have Evaluation as Pending.")




# %% Merge dataframes
# Merge class 1 and class 2 aggregated dataframes
merged_df = pd.merge(mhc1_agg_df, mhc2_agg_subset, on="ID")
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
merged_all.drop(columns=["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "MT Epitope Seq", "HLA Allele"])

# Print the resulting DataFrame
#print("Merged DataFrame:")
#print(merged_all.head())


# %% Reduce the number of columns in the merged_all DataFrame


## ---------------------------------------------------------------------------------------------------------------
# this section use to be for getting driver gene status, but it will be replaced with Susanna's updated edits to the all.epitopes.tsv file
# the new all.epitopes.tsv file will have a column called "Gene of Interest"
# this should contain "True" or "False" for the status of driver_gene


## ---------------------------------------------------------------------------------------------------------------

# %% Transformations (corresponds to 04_reduce_columns_for_training.py)
"""
Based on 04_reduce_columns_for_training.Rmd: 
New columns added: 
* `Pos` for values with the form integer-integer (e.g. 2-10), only keep the first position
* `Prob.Pos` is the position of mutation, changed to 0 if the value is "None", 1 otherwise
    * Change "None" in the `Prob.Pos` column to 0.
    * Create a new column `Prob.match` based on whether the integer(s) in `Prob.Pos` match the value in the `Pos` column.
    * Modify the `Prob.Pos` column to keep only the first integer if there are multiple.
* `Ref.Match` = factor(Ref.Match, levels = 0:1, labels = c("FALSE", "TRUE"))
* `TSL = ifelse(is.na(TSL), 6, TSL)`, change NAs into 6 category
* `Evaluation` changed to 1 if "Accept", 0 otherwise (also changed to a factor variable)
* `driver_gene`: now changed to "Gene of Interest", 1 if "Gene" is a driver gene (in the metadata), 0 otherwise
"""


# Transformations
# NOTE: Pos may take form "#-#", so we need to extract the first integer
if merged_all["Pos"].dtype == "object":
    merged_all["Pos"] = merged_all["Pos"].str.extract(r"^(\d+)").astype("Int64")  # Extract the first integer from Pos
# Convert Pos to float64
merged_all['Pos'] = merged_all['Pos'].astype(float)



# NOTE: Prob Pos may take form "#,#", so we need to extract the first integer
merged_all["Prob Pos"] = (
    merged_all["Prob Pos"]
    .fillna("0")  # Replace NaN with "0"
    .astype(str)  # Ensure all values are strings
    .replace("None", "0")  # Replace "None" with "0"
    .str.split(",")  # Split by commas
    .apply(lambda x: int(float(x[0])) if x[0].replace('.', '', 1).isdigit() else 0)  # Handle floats and integers
)
# Check for values in the format of #-# and keep only the first integer
#merged_all['Prob Pos'] = merged_all['Prob Pos'].apply(
    #lambda x: int(x.split('-')[0]) if isinstance(x, str) and '-' in x else x
#)

# Create a new column `Prob.match` based on whether the integer(s) in `Prob.Pos` match the value in the `Pos` column
#merged_all["Prob match"] = merged_all.apply(
    #lambda row: 1 if pd.notna(row["Pos"]) and row["Pos"] in [int(x) for x in str(row["Prob Pos"]).split(",") if x.isdigit()] else 0,
    #axis=1
#) # modified to handel Na in Pos column
merged_all["Prob match"] = merged_all.apply(
    lambda row: "True" if pd.notna(row["Pos"]) and row["Pos"] in [int(x) for x in str(row["Prob Pos"]).split(",") if x.isdigit()] else "False",
    axis=1
)
# If they are currently strings "True"/"False", convert them:
merged_all['Prob match'] = merged_all['Prob match'].map({'True': True, 'False': False}).astype(bool)
merged_all['Gene of Interest'] = merged_all['Gene of Interest'].map({'True': True, 'False': False}).astype(bool)
# Replace NA in TSL with 6, change into int
merged_all["TSL"] = merged_all["TSL"].fillna(6).astype(int)





# Convert Ref Match to lowercase and map to 0/1
# NOTE: In Python, True is equivalent to 1 and False is equivalent to 0. 
# Using .astype(int) converts the boolean values directly to integers.

#merged_all['Ref Match'] = merged_all['Ref Match'].str.upper().map({'FALSE': 0, 'TRUE': 1})


# Create a new column `driver_gene` based on whether the Gene column matches any value in the driver list
#driver_genes = set(driver.split(','))  # Assuming `driver` is a comma-separated string of driver genes
#merged_all['driver_gene'] = merged_all['Gene'].apply(lambda gene: 1 if gene in driver_genes else 0)


# select columns to keep from merged_all
columns_keep = [
    "ID",
    "Evaluation",
    "TSL",
    "Pos",
    "Prob Pos",
    "Num Passing Peptides",
    "RNA Expr",
    "RNA VAF",
    "Allele Expr",
    "RNA Depth",
    "DNA VAF",
    "Ref Match",
    "Biotype",
    "Variant Type",
    "Peptide Length",
    "Best MT IC50 Score",
    "Corresponding WT IC50 Score",
    "Corresponding Fold Change",
    "Best MT Percentile",
    "Corresponding WT Percentile",
    "Median Fold Change",
    "cysteine_count"
]

# Select the final columns
final_columns = columns_keep + [
    col for col in merged_all.columns if col.startswith(("IC50", "%ile", "MHCflurry", "MHCnuggetsI", "NetMHC", "PickPocket", "SMM"))
] + ['Prob match', 'Gene of Interest']

merged_all = merged_all[final_columns]

# Unify the column names
merged_all.rename(columns={
    "NetMHCpanEL WT Score": "NetMHCpanEL WT IC50 Score",
    "NetMHCpanEL MT Score": "NetMHCpanEL MT IC50 Score",
    "MHCflurry WT Score": "MHCflurry WT IC50 Score",
    "MHCflurry MT Score": "MHCflurry MT IC50 Score",
    "MHCnuggetsI WT Score": "MHCnuggetsI WT IC50 Score",
    "MHCnuggetsI MT Score": "MHCnuggetsI MT IC50 Score",
    "NetMHC WT Score": "NetMHC WT IC50 Score",
    "NetMHC MT Score": "NetMHC MT IC50 Score",
    "NetMHCcons WT Score": "NetMHCcons WT IC50 Score",
    "NetMHCcons MT Score": "NetMHCcons MT IC50 Score",
    "NetMHCpan WT Score": "NetMHCpan WT IC50 Score",
    "NetMHCpan MT Score": "NetMHCpan MT IC50 Score",
    "PickPocket WT Score": "PickPocket WT IC50 Score",
    "PickPocket MT Score": "PickPocket MT IC50 Score",
    "SMM WT Score": "SMM WT IC50 Score",
    "SMM MT Score": "SMM MT IC50 Score",
    "SMMPMBEC WT Score": "SMMPMBEC WT IC50 Score",
    "SMMPMBEC MT Score": "SMMPMBEC MT IC50 Score"
}, inplace=True)



# ---------------------------------------------------------------------------------------------------------------

# %% Cleaning up the merged_all DataFrame to prepare for prediction with the ML model (corresponds to 05_cleaning_imputation.py)

# --- Load trained imputer and label encoders (adjust paths as needed) ---
imputer_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/model_artifacts/trained_imputer.joblib")
label_encoders_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/model_artifacts/label_encoders.pkl")

imputer = joblib.load(imputer_path)
with open(label_encoders_path, "rb") as f:
    label_encoders = pickle.load(f)

# --- Prepare columns for imputation ---

exclude_columns = ["ID", "Evaluation"]  # Add any other columns you want to exclude
columns_to_impute = merged_all.columns.difference(exclude_columns)

excluded_data = merged_all[exclude_columns].copy()
data_to_impute = merged_all[columns_to_impute].copy()

# --- Apply label encoding to categorical columns ---
categorical_columns = data_to_impute.select_dtypes(include=['category', 'object']).columns
for col in categorical_columns:
    if col in label_encoders:
        le = label_encoders[col]
        # Handle unseen categories by using the most common category
        data_to_impute.loc[:, col] = data_to_impute[col].map(
            lambda x: le.transform([x])[0] if x in le.classes_ else le.transform([le.classes_[0]])[0]
        )

# --- Impute missing values ---
imputed_data = imputer.transform(data_to_impute)
imputed_data = pd.DataFrame(imputed_data, columns=columns_to_impute)

# --- Combine imputed and excluded columns ---
post_imputed_data = pd.concat([excluded_data.reset_index(drop=True), imputed_data.reset_index(drop=True)], axis=1)


# %% Predictions
# --- Load the trained Random Forest model ---
model_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/model_artifacts/rf_model_final.pkl")
rf_model = joblib.load(model_path)


# --- Prepare data for prediction ---
# If 'ID' or 'Evaluation' are in post_imputed_data, drop them for prediction
predict_cols = [col for col in post_imputed_data.columns if col not in ['ID', 'Evaluation', 'patient_id']]

# Make prediction with rf_model
rf_pred = rf_model.predict_proba(post_imputed_data[predict_cols])[:, 1]

# Add predictions to DataFrame
post_imputed_data['Accept_pred_prob'] = rf_pred
post_imputed_data['Evaluation_pred'] = np.where(
    post_imputed_data['Accept_pred_prob'].isna(), "Pending",
    np.where(
        post_imputed_data['Accept_pred_prob'] >= 0.55, "Accept",
        np.where(
            post_imputed_data['Accept_pred_prob'] > 0.30, "Review", "Reject"
        )
    )
)

# --- Merge with original ITB review file ---
mhc1_dir = Path(f"{pvacseq_results_dir}/mhc_i/")
mhc1_agg_df_itb = find_and_read_file(mhc1_dir, "*.all_epitopes.aggregated.tsv", "class 1 aggregated file")
#itb_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/data/itb_review")
#itb_file_path = list(itb_path.glob(f"*{patient_id}*.tsv"))
#if not itb_file_path:
    #raise FileNotFoundError(f"No itb_review file found for patient: {patient_id}")
#itb_file = pd.read_csv(itb_file_path[0], sep="\t", dtype=str) # force all columns as string

# Join the predicted Evaluation back to the original itb_review file
final_df = mhc1_agg_df_itb.drop(columns=['Evaluation']).merge(
    post_imputed_data[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
).rename(columns={'Evaluation_pred': 'Evaluation'})
final_df['Comments'] = "Probability of Accept: " + final_df['Accept_pred_prob'].round(3).astype(str)
final_df = final_df.drop(columns=['Accept_pred_prob'])

# Set Evaluation to "Pending" if missing
final_df.loc[(final_df['Evaluation'].isna()), 'Evaluation'] = "Pending"
final_df.loc[final_df['Evaluation'] == "Pending", 'Comments'] = "Unable to make prediction with ML model"




# %%
# --- Export the result ---
outdir = patient_dir / "ml_predicted_evaluations"
os.makedirs(outdir, exist_ok=True)
out_name = f"{patient_id}_predict_newThreshold_pvacview.tsv"
final_df.to_csv(outdir / out_name, sep="\t", index=False, na_rep="NA")
