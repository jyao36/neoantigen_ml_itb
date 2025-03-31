# %% 
import pandas as pd
import numpy as np
from pathlib import Path
import os
from multiprocessing import Pool
import re
import multiprocessing
# %%
# Set the working directory
project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
os.chdir(project_dir)

# Define output directory
outdir = project_dir / "output_python" / "03_merge_data.py"

os.makedirs(outdir, exist_ok=True)

print(f"Output directory created at: {outdir}")

# %%
# Load metadata
meta = pd.read_csv(project_dir / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv")
#meta.head()

# Get a list of patient IDs for ML study and external validation
patient_id_list = meta.loc[
    (meta["include_ML"] == "Yes") | (meta["external_validation"] == "Yes"),
    "patient_id"
].tolist()
print(patient_id_list)

# %%
def merge_patient_df(patient_id, project_dir, out_dir):
    try:
        # Find the patient's file in itb_review folder using regex
        itb_path = project_dir / "data" / "itb_review"
        itb_file_path = list(itb_path.glob(f"*{patient_id}*.tsv"))
        
        if not itb_file_path:
            print(f"Patient ID {patient_id} not found in the directory. Skipping...")
            return None
        else:
            print(f"Patient ID {patient_id} merging files")

        # Load ITB review file
        itb_file = pd.read_csv(itb_file_path[0], sep="\t").rename(columns={
            "Best Peptide": "Best Peptide class1",
            "Best Transcript": "Best Transcript class1",
            "IC50 MT": "IC50 MT class1",
            "IC50 WT": "IC50 WT class1",
            "%ile MT": "percentile MT class1",
            "%ile WT": "percentile WT class1"
        })

        # Load all_epitopes file
        epi_path = project_dir / "data" / "all_epitopes"
        epi_file_path = list(epi_path.glob(f"*{patient_id}*.tsv"))
        epi_file = pd.read_csv(epi_file_path[0], sep="\t")

        # Load class2 file and rename columns to class 2
        c2_path = project_dir / "data" / "class2"
        c2_file_path = list(c2_path.glob(f"*{patient_id}*.tsv"))
        c2_file = pd.read_csv(c2_file_path[0], sep="\t").rename(columns={
            "Best Peptide": "Best Peptide class2",
            "Best Transcript": "Best Transcript class2",
            "IC50 MT": "IC50 MT class2",
            "IC50 WT": "IC50 WT class2",
            "%ile MT": "percentile MT class2",
            "%ile WT": "percentile WT class2"
        })
        # Filter columns to include ID, columns starting with "D" except "DNA VAF", and renamed columns
        columns_to_keep = ["ID", "Best Peptide class2", "Best Transcript class2", "IC50 MT class2", "IC50 WT class2", "percentile MT class2", "percentile WT class2"]
        columns_to_keep += [col for col in c2_file.columns if col.startswith("D") and col != "DNA VAF"]
        c2_file = c2_file.filter(items=columns_to_keep)

        # Check if itb_file and c2_file have the same number of rows
        if itb_file.shape[0] != c2_file.shape[0]:
            print("itb_review file and class2 file rows DO NOT match, but continue processing...")

        # Merge class2 file to itb file
        itb_c2 = pd.merge(itb_file, c2_file, on="ID")

        # Join epi file columns
        # Separate values in the "ID" column into separate columns
        itb_c2[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']] = itb_c2['ID'].str.split('-', expand=True)
        
        # Convert "Start" and "Stop" columns to integers
        itb_c2['Start'] = itb_c2['Start'].astype(int)
        itb_c2['Stop'] = itb_c2['Stop'].astype(int)
        
        # Merge class2 file to itb file
        itb_c2_epi = pd.merge(
            itb_c2, epi_file,
            left_on=["Chromosome", "Start", "Stop", "Reference", "Variant", "Best Transcript class1", "Best Peptide class1", "Allele"],
            right_on=["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "MT Epitope Seq", "HLA Allele"],
        how='inner'
    )

        
        # Remove the redundant columns that were made from splitting ID
        # NOTE: This step is necessary to avoid duplicate columns when saving the output (R removes the redundant columns "Transcript", "MT Epitope Seq", "HLA Allele" automatically after joining)
        itb_c2_epi = itb_c2_epi.drop(columns=["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "MT Epitope Seq", "HLA Allele"])

        
        
        # Check that the row numbers match between the merged file and original itb_file
        if itb_c2_epi.shape[0] != itb_file.shape[0]:
            print("itb_review file and merged file rows DO NOT match, but continue saving...")

        # Generate file name
        out_name = f"{patient_id}_merged.tsv"
        # Save the output into the outdir
        itb_c2_epi.to_csv(out_dir / out_name, sep="\t", index=False)
    except Exception as e:
        print(f"Error processing patient ID {patient_id}: {e}")

# %%
# Test the function with one patient
#merge_patient_df(patient_id_list[0], project_dir, outdir)
#merge_patient_df(patient_id_list[1], project_dir, outdir)
#merge_patient_df("CTEP-10146-MD017-0052", project_dir, outdir)

#%%
# Get the number of available CPU cores
#num_cores = multiprocessing.cpu_count()

# Run the function on all patients
#with Pool(num_cores-4) as p:
#    p.starmap(merge_patient_df, [(pid, project_dir, outdir) for pid in patient_id_list])
# %%
results = [merge_patient_df(pid, project_dir, outdir) for pid in patient_id_list]
# %%
