# %% Summary of the script: 
# NOTE: This file is for training samples only! Not used during prediction (for prediction use script 05-2_cleaning_prediciton.py)

# Goal:  

# * Merge all patient's data into one  
# * Impute missing NA values
# * Merge samples together for later training and testing

# %% Import dependencies
from pathlib import Path
import pandas as pd
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.experimental import enable_iterative_imputer  # This will enable the experimental feature and allow you to use IterativeImputer
from sklearn.impute import IterativeImputer
import joblib  # Add this import
import pickle

# %% Set working and output directory ---------------------------------------------------------------------------------------------------------
# Set the working directory
project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
os.chdir(project_dir)

# Define output directory
outdir = project_dir / "output_python" / "05_cleaning_imputation.py"
os.makedirs(outdir, exist_ok=True)

# %% Load metadata ---------------------------------------------------------------------------------------------------------
meta = pd.read_csv(project_dir / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv")

# Replace "." with " " in metadata column names
meta.columns = meta.columns.str.replace('.', ' ', regex=False)

# Get a list of patient IDs for ML study
patient_id_list = meta.loc[meta["include_ML"] == "Yes", "patient_id"].tolist()
#print(f"Number of patients included in analysis: {len(patient_id_list)}")


# %% Define function to process each file and unify column names -------------------------------------------------------------------------
def merge_files_ml2(patient_id_test, in_dir):
    # Locate the reduced file for the patient
    reduced_file_path = list(Path(in_dir).glob(f"*{patient_id_test}*.tsv"))
    if not reduced_file_path:
        print(f"File for Patient ID {patient_id_test} not found. Skipping...")
        return None

    # Load the reduced file
    reduced_data = pd.read_csv(reduced_file_path[0], sep="\t")
    reduced_data["patient_id"] = patient_id_test

    # Unify the column names
    reduced_data.rename(columns={
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

    return reduced_data

# %% Process all files and combine into a single dataframe -------------------------------------------------------------------------
in_dir = project_dir / "output_python" / "04_reduce_columns_for_training.py"
all_files2 = [merge_files_ml2(patient_id, in_dir) for patient_id in patient_id_list]
merged_df2 = pd.concat([df for df in all_files2 if df is not None], ignore_index=True)


# %% Replace NA with existing values
replace_map = [
    ("NetMHC WT IC50 Score", "IC50 WT class1"),
    ("SMM WT IC50 Score", "IC50 WT class1"),
    ("SMMPMBEC WT IC50 Score", "IC50 WT class1"),
    ("NetMHC MT IC50 Score", "IC50 MT class1"),
    ("SMM MT IC50 Score", "IC50 MT class1"),
    ("SMMPMBEC MT IC50 Score", "IC50 MT class1"),
    ("NetMHC WT Percentile", "percentile WT class1"),
    ("SMM WT Percentile", "percentile WT class1"),
    ("SMMPMBEC WT Percentile", "percentile WT class1"),
    ("NetMHC MT Percentile", "percentile MT class1"),
    ("SMM MT Percentile", "percentile MT class1"),
    ("SMMPMBEC MT Percentile", "percentile MT class1"),
]

for target_col, backup_col in replace_map:
    if target_col in merged_df2.columns and backup_col in merged_df2.columns:
        merged_df2[target_col] = merged_df2[target_col].fillna(merged_df2[backup_col])


# %% The following chunk is used to compare to the R version, comment out the following code block
"""
reduced_data_to_compare_path = project_dir / "output_python" / "transform_training_data4python.ipynb" / 'pre_imputed_data.csv'
reduced_data_to_compare = pd.read_csv(reduced_data_to_compare_path)

# Compare column names between merged_df2 and reduced_data_to_compare
merged_columns = set(merged_df2.columns)
reduced_columns = set(reduced_data_to_compare.columns)

# Find columns present in merged_df2 but not in reduced_data_to_compare
only_in_merged = merged_columns - reduced_columns
print(f"Columns only in merged_df2: {only_in_merged}")

# Find columns present in reduced_data_to_compare but not in merged_df2
only_in_reduced = reduced_columns - merged_columns
print(f"Columns only in reduced_data_to_compare: {only_in_reduced}")

# Find common columns
common_columns = merged_columns & reduced_columns
print(f"Common columns: {common_columns}")

# Compare column types between merged_df2 and reduced_data_to_compare
merged_dtypes = merged_df2.dtypes
reduced_dtypes = reduced_data_to_compare.dtypes

# Find columns with differing types
differing_types = {
    col: (merged_dtypes[col], reduced_dtypes[col])
    for col in merged_columns & reduced_columns
    if merged_dtypes[col] != reduced_dtypes[col]
}

if differing_types:
    print("Columns with differing types:")
    for col, types in differing_types.items():
        print(f"{col}: merged_df2 type = {types[0]}, reduced_data_to_compare type = {types[1]}")
else:
    print("All common columns have the same types.")

# Set pandas options to display all rows and columns
pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.max_colwidth', None)  # Show full column content
"""


# %% Encode labels for categorical columns 
categorical_columns = merged_df2.select_dtypes(include=['category', 'object']).columns.drop(['Evaluation', 'ID', 'patient_id'])
print("Categorical columns:", categorical_columns)

# Label encode categorical columns
label_encoders = {}
for col in categorical_columns:
    le = LabelEncoder()
    merged_df2[col] = le.fit_transform(merged_df2[col])
    label_encoders[col] = le  # Store the encoder for future use

with open(outdir / "label_encoders.pkl", "wb") as f:
    pickle.dump(label_encoders, f)
# dump another copy of it in output_python/model_artifacts
with open(project_dir / "output_python" / "model_artifacts" / "label_encoders.pkl", "wb") as f:
    pickle.dump(label_encoders, f)

# %% Impute missing values ----------------------------------------------------------------------------------
# Exclude columns from imputation
exclude_columns = ["ID", "patient_id", "Evaluation"]
columns_to_impute = merged_df2.columns.difference(exclude_columns)

# Separate the columns to exclude
excluded_data = merged_df2[exclude_columns]
data_to_impute = merged_df2[columns_to_impute]

# Initialize the IterativeImputer
imputer = IterativeImputer(random_state=918)

# Impute missing values only in the selected columns
imputed_data = imputer.fit_transform(data_to_impute)

# Save the trained imputer for 05-2_cleaning_prediction.py
joblib.dump(imputer, outdir / "trained_imputer.joblib")
# dump another copy of it in output_python/model_artifacts
joblib.dump(imputer, project_dir / "output_python" / "model_artifacts" / "trained_imputer.joblib")

# Convert the imputed data back to a DataFrame
imputed_data = pd.DataFrame(imputed_data, columns=columns_to_impute)

# Combine the excluded columns back with the imputed data
post_imputed_data = pd.concat([excluded_data.reset_index(drop=True), imputed_data.reset_index(drop=True)], axis=1)

# save data
post_imputed_data.to_csv(outdir / "df_imputed.csv", index=False)
# %%


