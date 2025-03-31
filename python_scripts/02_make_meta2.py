# %% import libraries
# import libraries 
import pandas as pd
import numpy as np
from pathlib import Path
import os
from typing import Optional

#%% Set directories 
# Set the working directory
# Ensure project_dir is correctly defined
project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
os.chdir(project_dir)

# Define output directory
outdir = Path.cwd() / "output_python" / "02_make_meta2.ipynb"
os.makedirs(outdir, exist_ok=True)

print(f"Output directory created at: {outdir}")


# Define file path
file_path = Path("data/meta_manual.csv")

# Load the CSV file into a pandas DataFrame
meta_manual = pd.read_csv(file_path)

# Display the DataFrame
#meta_manual.head()  # Shows the first few rows


# Filter rows where include_ML is "Yes" or external_validation is "Yes"
patient_id_ml = meta_manual.loc[
    (meta_manual["include_ML"] == "Yes") | (meta_manual["external_validation"] == "Yes"),
    "patient_id"
].tolist()

# Print or check the output
print(patient_id_ml)

# %%
# Update counts ----------------------------------------------------------------------------------------

def extract_num_peptides(patient_id, project_dir):
    # Path to the peptide review files
    itb_path = project_dir / "data" / "itb_review"

    # Find files that match the patient ID
    itb_files = list(itb_path.glob(f"*{patient_id}*.tsv"))

    if not itb_files:
        print(f"[SKIP] No file found for patient {patient_id}")
        return None

    print(f"[PROCESS] Counting peptides for patient {patient_id}")

    try:
        # Attempt something that might raise an error
        # Read the first matching file
        itb_file = pd.read_csv(itb_files[0], sep="\t")

        # Compute counts
        n_peptides = len(itb_file)
        n_pass = (itb_file['Tier'] == "Pass").sum() if 'Tier' in itb_file.columns else None
        n_accept = (itb_file['Evaluation'] == "Accept").sum() if 'Evaluation' in itb_file.columns else None

        return pd.DataFrame({
            'patient_id': [patient_id],
            'n_peptides': [n_peptides],
            'n_pass': [n_pass],
            'n_accept': [n_accept]
        })

    except Exception as e: # Exception is the general type of error. as e gives you the actual error object (e), so you can print or log the message like: file not found, column missing, etc.
        print(f"[ERROR] Failed to process {patient_id}: {e}")
        return None
    
#extract_num_peptides('TWJF-5120-06', project_dir)

# Apply to each patient and concatenate results
count_df_list = [extract_num_peptides(pid, project_dir) for pid in patient_id_ml]
count_df = pd.concat([df for df in count_df_list if df is not None], ignore_index=True)

# Merge back into the original metadata
meta_count = meta_manual.merge(count_df, on="patient_id", how="left")

# Display the final DataFrame
#meta_count.head()





# %% get purity
def extract_purity_from_driver(patient_id, project_dir, meta_count):
    # Step 1: Get driver genes for this patient
    match = meta_count[meta_count['patient_id'] == patient_id]
    if match.empty:
        print(f"[SKIP] No metadata found for patient {patient_id}")
        return None

    driver_genes_str = match.iloc[0].get('driver_genes', None)
    if pd.isna(driver_genes_str):
        print(f"[SKIP] No driver genes for patient {patient_id}")
        return None

    # Step 2: Locate the file
    itb_path = project_dir / "data" / "itb_review"
    itb_files = list(itb_path.glob(f"*{patient_id}*.tsv"))
    if not itb_files:
        print(f"[SKIP] No file found for patient {patient_id}")
        return None

    print(f"[PROCESS] Calculating purity for patient {patient_id}")
    try:
        df = pd.read_csv(itb_files[0], sep="\t")

        # Ensure columns exist
        if 'DNA VAF' not in df.columns or 'Gene' not in df.columns:
            print(f"[WARN] Missing DNA VAF or Gene column for {patient_id}")
            return None

        # Convert VAF safely
        df['DNA VAF'] = pd.to_numeric(df['DNA VAF'], errors='coerce')

        # Parse driver genes
        driver_genes = [g.strip() for g in str(driver_genes_str).split(',')]

        # Subset to only driver genes
        driver_vaf_df = df[df['Gene'].isin(driver_genes)][['ID', 'Gene', 'DNA VAF']].dropna()

        vaf_for_purity = pd.DataFrame()

        # === Begin Decision Logic ===
        if len(driver_vaf_df) == 1 and driver_vaf_df['DNA VAF'].iloc[0] < 0.5:
            vaf_for_purity = driver_vaf_df.copy()
            vaf_for_purity['purity_from_driver'] = "YES"

        elif len(driver_vaf_df) == 1 and driver_vaf_df['DNA VAF'].iloc[0] >= 0.5:
            below_05 = df[df['DNA VAF'] < 0.5].sort_values('DNA VAF', ascending=False)
            if not below_05.empty:
                highest_vaf_row = below_05[['ID', 'Gene', 'DNA VAF']].iloc[0:1].copy()
                highest_vaf_row['purity_from_driver'] = "NO"
                vaf_for_purity = highest_vaf_row

        elif len(driver_vaf_df) > 1 and any(driver_vaf_df['DNA VAF'] < 0.5):
            below_05 = driver_vaf_df[driver_vaf_df['DNA VAF'] < 0.5].sort_values('DNA VAF', ascending=False)
            if not below_05.empty:
                highest = below_05.iloc[0:1].copy()
                highest['purity_from_driver'] = "YES"
                vaf_for_purity = highest

        elif len(driver_vaf_df) > 1 and not any(driver_vaf_df['DNA VAF'] < 0.5):
            below_05 = df[df['DNA VAF'] < 0.5].sort_values('DNA VAF', ascending=False)
            if not below_05.empty:
                highest = below_05[['ID', 'Gene', 'DNA VAF']].iloc[0:1].copy()
                highest['purity_from_driver'] = "NO"
                vaf_for_purity = highest

        elif len(driver_vaf_df) == 0:
            below_05 = df[df['DNA VAF'] < 0.5].sort_values('DNA VAF', ascending=False)
            if not below_05.empty:
                highest = below_05[['ID', 'Gene', 'DNA VAF']].iloc[0:1].copy()
                highest['purity_from_driver'] = "NO"
                vaf_for_purity = highest

        # Finalize output
        if not vaf_for_purity.empty:
            vaf_for_purity['purity'] = vaf_for_purity['DNA VAF'] * 2
            vaf_for_purity['patient_id'] = patient_id
            return vaf_for_purity[['patient_id', 'ID', 'Gene', 'DNA VAF', 'purity_from_driver', 'purity']]

        else:
            print(f"[WARN] No valid purity row for patient {patient_id}")
            return None

    except Exception as e:
        print(f"[ERROR] {patient_id}: {e}")
        return None
# %% test the function
#extract_purity_from_driver("TWJF-5120-06", project_dir, meta_count)
# %%
purity_df_list = [
    extract_purity_from_driver(pid, project_dir, meta_count)
    for pid in meta_count['patient_id']
]

purity_df = pd.concat([df for df in purity_df_list if df is not None], ignore_index=True)


# Merge the purity results back into the main metadata table
# Merge full purity information into the metadata
meta_full = meta_count.merge(purity_df, on='patient_id', how='left')


# %% Save output 

# Define output file path (if not already defined)
output_path = outdir / "metadata_count_purity.xlsx"

# Save to Excel
meta_full.to_excel(output_path, index=False)

print(f"Saved to: {output_path}")

# Define output file path for CSV
output_path_csv = outdir / "metadata_count_purity.csv"

# Save to CSV
meta_full.to_csv(output_path_csv, index=False)

print(f"Saved to: {output_path_csv}")