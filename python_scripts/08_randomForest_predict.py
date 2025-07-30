#%% Load libraries ----------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import joblib
from pathlib import Path
import os
from sklearn.preprocessing import LabelEncoder
import pickle

#%% Load the trained Random Forest model and metadata -----------------------------------------------------
# Load the trained Random Forest model
#model_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/model_artifacts/rf_model_final.pkl")
model_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/transform_training_data_model.ipynb/rf_model.pkl")
rf_model = joblib.load(model_path)

# Load the trained imputer and label encoders
imputer_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/model_artifacts/trained_imputer.joblib")
label_encoders_path = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project/output_python/model_artifacts/label_encoders.pkl")

imputer = joblib.load(imputer_path)
with open(label_encoders_path, "rb") as f:
    trained_label_encoders = pickle.load(f)

# Define directories
project_dir = Path("/Users/Jennie/Desktop/WashU/Rotation_labs/Griffith Lab/Neoantigen ML project")
in_dir = project_dir / "output_python" / "05-2_cleaning_prediction.py"
outdir = project_dir / "output_python" / "08_randomForest_predict.py"
os.makedirs(outdir, exist_ok=True)

# Read metadata
meta = pd.read_csv(project_dir / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv")
# If running with R, use the following line instead:
#meta = pd.read_csv(project_dir / "output" / "02_make_meta2.Rmd" / "metadata_count_purity.csv")

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
    cleaned_file_path = list(in_dir.glob(f"*{patient_id}*.csv"))
    if not cleaned_file_path:
        raise FileNotFoundError(f"No cleaned data file found for patient: {patient_id}")
    cleaned_data = pd.read_csv(cleaned_file_path[0])
    
    # --- Imputation section ---
    # Prepare columns for imputation
    exclude_columns = ["ID", "Evaluation", "patient_id"]  # Add any other columns you want to exclude
    columns_to_impute = cleaned_data.columns.difference(exclude_columns)
    
    excluded_data = cleaned_data[exclude_columns].copy()
    data_to_impute = cleaned_data[columns_to_impute].copy()
    
    # Apply label encoding to categorical columns for imputation
    categorical_columns = data_to_impute.select_dtypes(include=['category', 'object']).columns
    for col in categorical_columns:
        if col in trained_label_encoders:
            le = trained_label_encoders[col]
            # Handle unseen categories by using the most common category
            data_to_impute.loc[:, col] = data_to_impute[col].map(
                lambda x: le.transform([x])[0] if x in le.classes_ else le.transform([le.classes_[0]])[0]
            )
    
    # Impute missing values
    imputed_data = imputer.transform(data_to_impute)
    imputed_data = pd.DataFrame(imputed_data, columns=columns_to_impute)
    
    # Combine imputed and excluded columns
    cleaned_data = pd.concat([excluded_data.reset_index(drop=True), imputed_data.reset_index(drop=True)], axis=1)
    
    # Identify categorical columns for final encoding
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
            final_pred['Accept_pred_prob'] >= 0.55, "Accept",
            np.where(
                final_pred['Accept_pred_prob'] > 0.30, "Review", "Reject"
            )
        )
    )
    
    # Read the original itb_review file
    itb_path = project_dir / "data" / "itb_review"
    itb_file_path = list(itb_path.glob(f"*{patient_id}*.tsv"))
    if not itb_file_path:
        raise FileNotFoundError(f"No itb_review file found for patient: {patient_id}")
    itb_file = pd.read_csv(itb_file_path[0], sep="\t", dtype=str) # force all columns as string to avoid dtype issues
    
    # Join the predicted Evaluation back to the original itb_review file
    final_df = itb_file.drop(columns=['Evaluation']).merge(
        final_pred[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
    ).rename(columns={'Evaluation_pred': 'Evaluation'})
    final_df['Comments'] = "Probability of Accept: " + final_df['Accept_pred_prob'].round(3).astype(str)
    final_df = final_df.drop(columns=['Accept_pred_prob'])
    
    # There are cases where the Evaluation is NA and Accept_pred_prob is NaN. 
    # This may happen when class1 and class2 files have different number of rows. 
    # (usually when class1 aggregated file has more rows than class2, but in 03_merge_data.py the rows are removed since the model needs all info for prediction)
    # In this case, we will set the Evaluation to "Pending"
    final_df.loc[(final_df['Evaluation'].isna()), 'Evaluation'] = "Pending"
    # For rows where Evaluation is "Pending", set Comments accordingly
    final_df.loc[final_df['Evaluation'] == "Pending", 'Comments'] = "Unable to make prediction with ML model"
    
    # Export the itb_review file with new predictions
    out_name = f"{patient_id}_predict_newThreshold_pvacview.tsv"
    final_df.to_csv(outdir / out_name, sep="\t", index=False, na_rep="NA") # Export as TSV with NA rather than empty cells
    
    # Make another version that keeps the Reject and Accept predicted probabilities
    final_df2 = itb_file.drop(columns=['Evaluation']).merge(
        final_pred[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
    )
    out_name2 = f"{patient_id}_predict_newThreshold_prob.tsv"
    final_df2.to_csv(outdir / out_name2, sep="\t", index=False, na_rep="NA")
    

# For testing
#rf_predict_on_new("10146-0061", project_dir, in_dir, outdir, label_encoders)
# Apply the function to all patient IDs
for patient_id in patient_id_list:
    rf_predict_on_new(patient_id, project_dir, in_dir, outdir, label_encoders)
