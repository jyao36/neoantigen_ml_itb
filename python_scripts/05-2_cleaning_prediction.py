# Note: This file is for external validation/prediction samples only! Not used during training (for training use script 05_cleaning_imputation.py)

# %%
from pathlib import Path
import pandas as pd
import joblib
import pickle
import os
from sklearn.preprocessing import LabelEncoder

# %%
def process_prediction_samples(patient_id, project_dir, in_dir, outdir):
    """
    Process and impute missing values for prediction samples using the trained imputer.
    
    Parameters:
    -----------
    patient_id : str
        Patient ID to process
    project_dir : Path
        Project root directory
    in_dir : Path
        Directory containing input files
    outdir : Path
        Directory to save output files
    """
    try:
        # Load the trained imputer and label encoders
        imputer_path = project_dir / "output_python" / "05_cleaning_imputation.py" / "trained_imputer.joblib"
        label_encoders_path = project_dir / "output_python" / "05_cleaning_imputation.py" / "label_encoders.pkl"
        
        imputer = joblib.load(imputer_path)
        with open(label_encoders_path, "rb") as f:
            label_encoders = pickle.load(f)
        
        # Find and load the patient's file
        file_path = list(in_dir.glob(f"*{patient_id}*.tsv"))
        if not file_path:
            print(f"File for Patient ID {patient_id} not found. Skipping...")
            return None
        
        print(f"Processing Patient ID {patient_id}")
        df = pd.read_csv(file_path[0], sep="\t", low_memory=False)
        
        # Add patient_id column
        df["patient_id"] = patient_id
        
        # Unify the column names
        df.rename(columns={
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
            if target_col in df.columns and backup_col in df.columns:
                df[target_col] = df[target_col].fillna(df[backup_col])

        
        # Define columns to exclude from imputation
        exclude_columns = ["ID", "patient_id", "Evaluation"]
        columns_to_impute = df.columns.difference(exclude_columns)
        
        # Create explicit copies of the data
        excluded_data = df[exclude_columns].copy()
        data_to_impute = df[columns_to_impute].copy()
        
        # Apply label encoding to categorical columns
        categorical_columns = data_to_impute.select_dtypes(include=['category', 'object']).columns
        for col in categorical_columns:
            if col in label_encoders:
                le = label_encoders[col]
                # Handle unseen categories by using the most common category
                data_to_impute.loc[:, col] = data_to_impute[col].map(
                    lambda x: le.transform([x])[0] if x in le.classes_ else le.transform([le.classes_[0]])[0]
                )
        
        # Transform the data using the trained imputer
        imputed_data = imputer.transform(data_to_impute)
        
        # Convert back to DataFrame
        imputed_data = pd.DataFrame(imputed_data, columns=columns_to_impute)
        
        # Combine with excluded columns
        post_imputed_data = pd.concat([excluded_data.reset_index(drop=True), 
                                     imputed_data.reset_index(drop=True)], axis=1)
        
        # Save the imputed data
        out_name = f"{patient_id}_imputed.csv"
        post_imputed_data.to_csv(outdir / out_name, index=False)
        
        return post_imputed_data
        
    except Exception as e:
        print(f"Error processing patient ID {patient_id}: {e}")
        return None

# %%
def main():
    # Set up directories
    project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
    in_dir = project_dir / "output_python" / "04_reduce_columns_for_training.py"
    outdir = project_dir / "output_python" / "05-2_cleaning_prediction.py"
    os.makedirs(outdir, exist_ok=True)
    
    # Load metadata to get patient IDs for external validation
    meta = pd.read_csv(project_dir / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv")
    patient_id_list = meta.loc[meta["external_validation"] == "Yes", "patient_id"].tolist()
    
    # Process all patient IDs
    results = [process_prediction_samples(pid, project_dir, in_dir, outdir) for pid in patient_id_list]
    
    # Print summary
    successful = sum(1 for r in results if r is not None)
    print(f"\nProcessed {successful} out of {len(patient_id_list)} patients successfully")

if __name__ == "__main__":
    main()
# %%
