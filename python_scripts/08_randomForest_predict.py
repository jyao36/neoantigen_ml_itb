#%% Load libraries ----------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import joblib
from pathlib import Path
import os
from sklearn.preprocessing import LabelEncoder

#%% Load the trained Random Forest model and metadata -----------------------------------------------------
# Load the trained Random Forest model
model_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/07_ml_randomForest.ipynb/rf_model.pkl")
rf_model = joblib.load(model_path)

# Define directories
project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
in_dir = project_dir / "output" / "05-2_cleaning_prediction.Rmd"
outdir = project_dir / "output_python" / "08_randomForest_predict.py"
os.makedirs(outdir, exist_ok=True)

# Read metadata
meta = pd.read_csv(project_dir / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv")

# Get a list of patient IDs for external validation
patient_id_list = meta.loc[meta["external_validation"] == "Yes", "patient_id"].tolist()

# Define the label encoders and mappings used during training
label_encoders = {}

# Actual mappings used during training
class_mappings = {
    'Biotype': ['IG_V_gene', 'nonsense_mediated_decay', 'protein_coding'],
    'Variant.Type': ['FS', 'inframe_del', 'inframe_ins', 'missense'],
    'Prob.match': ['NO', 'YES'],
    'driver_gene': ['NO', 'YES']
}

for col, classes in class_mappings.items():
    le = LabelEncoder()
    le.classes_ = np.array(classes)
    label_encoders[col] = le

#%% Prediction (Accept VS Reject model) -------------------------------------------------------------------
def rf_predict_on_new(patient_id, project_dir, in_dir, outdir, label_encoders):
    # Read cleaned data that is used for prediction
    cleaned_file_path = list(in_dir.glob(f"*{patient_id}*.xlsx"))
    if not cleaned_file_path:
        raise FileNotFoundError(f"No cleaned data file found for patient: {patient_id}")
    cleaned_data = pd.read_excel(cleaned_file_path[0])
    
    # Identify categorical columns
    categorical_cols = cleaned_data.select_dtypes(include=['object', 'category']).columns
    
    # Apply label encoding to the categorical columns
    for col in categorical_cols:
        if col in label_encoders and isinstance(label_encoders[col], LabelEncoder):
            le = label_encoders[col]
            cleaned_data[col] = le.transform(cleaned_data[col])
    
    # Make prediction with rf_model
    rf_pred = rf_model.predict_proba(cleaned_data.drop(columns=['ID', 'patient_id', 'Evaluation']))[:, 1]
    
    # Change predicted probabilities into "Accept", "Reject", "Review" based on set thresholds
    final_pred = cleaned_data.copy()
    final_pred['Accept_pred_prob'] = rf_pred
    final_pred['Evaluation_pred'] = np.where(
        final_pred['Accept_pred_prob'].isna(), "Pending",
        np.where(
            final_pred['Accept_pred_prob'] >= 0.60, "Accept",
            np.where(
                final_pred['Accept_pred_prob'] > 0.20, "Review", "Reject"
            )
        )
    )
    
    # Read the original itb_review file
    itb_path = project_dir / "data" / "itb_review"
    itb_file_path = list(itb_path.glob(f"*{patient_id}*.tsv"))
    if not itb_file_path:
        raise FileNotFoundError(f"No itb_review file found for patient: {patient_id}")
    itb_file = pd.read_csv(itb_file_path[0], sep="\t")
    
    # Join the predicted Evaluation back to the original itb_review file
    final_df = itb_file.drop(columns=['Evaluation']).merge(
        final_pred[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
    ).rename(columns={'Evaluation_pred': 'Evaluation'})
    final_df['Comments'] = "Probability of Accept: " + final_df['Accept_pred_prob'].round(3).astype(str)
    final_df = final_df.drop(columns=['Accept_pred_prob'])
    
    # Export the itb_review file with new predictions
    out_name = f"{patient_id}_predict_newThreshold.tsv"
    final_df.to_csv(outdir / out_name, sep="\t", index=False)
    
    # Make another version that keeps the Reject and Accept predicted probabilities
    final_df2 = itb_file.drop(columns=['Evaluation']).merge(
        final_pred[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
    )
    out_name2 = f"{patient_id}_predict_newThreshold2.tsv"
    final_df2.to_csv(outdir / out_name2, sep="\t", index=False)

# Apply the function to all patient IDs
for patient_id in patient_id_list:
    rf_predict_on_new(patient_id, project_dir, in_dir, outdir, label_encoders)
