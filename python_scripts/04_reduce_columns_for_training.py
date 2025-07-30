# %%
# Summary of the script:
# - Sets the working directory and defines the output directory for reduced data.
# - Loads metadata from a CSV file and retrieves patient IDs for ML study and external validation.
# - Defines a function `reduce_columns` to:
#   - Locate and load patient-specific merged files.
#   - Extract driver genes for each patient from metadata.
#   - Perform data transformations, including:
#     - Extracting and converting specific columns (e.g., "Pos", "Prob Pos").
#     - Adding new columns like "Prob match" and "Gene of Interest".
#   - Select and retain relevant columns for training.
#   - Save the reduced data to the output directory.
# - Processes all patient IDs using the `reduce_columns` function.

# %% Import libraries ---------------------------------------------------------------------------------------------------------
from pathlib import Path
import pandas as pd
import os

# %% Set working and output directory ---------------------------------------------------------------------------------------------------------
# Set the working directory
project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
os.chdir(project_dir)

# Define output directory
outdir = project_dir / "output_python" / "04_reduce_columns_for_training.py"

os.makedirs(outdir, exist_ok=True)

# %% Load metadata ---------------------------------------------------------------------------------------------------------
# Load metadata
meta = pd.read_csv(project_dir / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv")

# %% Get a list of patient IDs for ML study and external validation ---------------------------------------------------------------------------------------------------------
# Get the list of patient IDs to process
patient_id_list = meta.loc[
    (meta["include_ML"] == "Yes") | (meta["external_validation"] == "Yes"),
    "patient_id"
].tolist()


# %%
# Define the reduce_columns function
def reduce_columns(patient_id_test, in_dir, out_dir):
    # Find the patient's file
    matching_files = [f for f in os.listdir(in_dir) if f"{patient_id_test}" in f and f.endswith(".tsv")]
    
    if not matching_files:
        print(f"Patient ID {patient_id_test} not found in the directory. Skipping...")
        return None
    else:
        print(f"Patient ID {patient_id_test} reducing files")
    
    # Read the merged file for the patient
    merged_full_file_path = next((Path(in_dir) / f for f in matching_files), None)
    merged_data = pd.read_csv(merged_full_file_path, sep="\t")
        
    # Extract driver genes for the patient from metadata
    driver = meta.loc[meta["patient_id"] == patient_id_test, "driver_genes"].values
    driver_genes = driver[0].split(",") if len(driver) > 0 else []
    
    # Columns to keep
    columns_keep = [
        "ID", "Evaluation", "TSL", "Pos", "Prob Pos", "Num Passing Peptides", "RNA Expr", "RNA VAF",
        "Allele Expr", "RNA Depth", "DNA VAF", "Ref Match", "Biotype", "Variant Type", "Peptide Length",
        "Best MT IC50 Score", "Corresponding WT IC50 Score", "Corresponding Fold Change", "Best MT Percentile",
        "Corresponding WT Percentile", "Median Fold Change", "cysteine_count"
    ]
    
    # Process the data
    # NOTE: Pos may take form "#-#", so we need to extract the first integer
    if merged_data["Pos"].dtype == "object":
        merged_data["Pos"] = merged_data["Pos"].str.extract(r"^(\d+)").astype("Int64")
    
    # NOTE: Prob Pos may take form "#,#", so we need to extract the first integer
    merged_data["Prob Pos"] = (
        merged_data["Prob Pos"]
        .fillna("0")  # Replace NaN with "0"
        .astype(str)  # Ensure all values are strings
        .replace("None", "0")  # Replace "None" with "0"
        .str.split(",")  # Split by commas
        .apply(lambda x: int(float(x[0])) if x[0].replace('.', '', 1).isdigit() else 0)  # Handle floats and integers
    )
    
    merged_data["Prob match"] = merged_data.apply(
        lambda row: "True" if pd.notna(row["Pos"]) and row["Pos"] in [int(x) for x in str(row["Prob Pos"]).split(",") if x.isdigit()] else "False",
        axis=1
    )
    
    merged_data["TSL"] = merged_data["TSL"].fillna(6).astype(int)
    
    #merged_data["Ref Match"] = merged_data["Ref Match"].str.upper().map({"FALSE": 0, "TRUE": 1})
    
    merged_data["driver_gene"] = merged_data["Gene"].apply(lambda gene: "True" if gene in driver_genes else "False")
    
    # Rename "driver_gene" to "Gene of Interest"
    merged_data.rename(columns={"driver_gene": "Gene of Interest"}, inplace=True)
    
    # Select columns
    new_data = merged_data[
        columns_keep +
        [col for col in merged_data.columns if col.startswith(("IC50", "%ile", "MHCflurry", "MHCnuggetsI", "NetMHC", "PickPocket", "SMM"))] + 
        ["Prob match", "Gene of Interest"]
    ]
    
    # Generate file names
    out_name = f"{patient_id_test}_reduced.tsv"
    
    # Save the output
    new_data.to_csv(out_dir / out_name, sep="\t", index=False)
    print(f"Reduced file saved for Patient ID {patient_id_test}")


in_dir = project_dir / "output_python" / "03_merge_data.py"
#reduce_columns("TWJF-10146-MO011-0024", in_dir, outdir)  # Test the function with one patient ID
# %% Process all patient IDs using multiprocessing ---------------------------------------------------------------------------------------------------------

results = [reduce_columns(patient_id, in_dir, outdir) for patient_id in patient_id_list]
# %%
