from pathlib import Path
import pandas as pd

def create_output_dir():
    """Create output directory"""
    outdir = Path.cwd() / "output_python" / "summary.py"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

def read_data():
    """Read metadata and imputed data"""
    meta_path = Path.cwd() / "output_python" / "02_make_meta2.ipynb" / "metadata_count_purity.csv"
    imputed_path = Path.cwd() / "output_python" / "05_cleaning_imputation.py" / "df_imputed.csv"
    # external validation data (imputated and actual Evaluation)
    external_validation_dir = Path.cwd() / "output_python" / "05-2_cleaning_prediction.py"
    itb_review_dir = Path.cwd() / "data" / "itb_review_new"
    
    # Read metadata CSV
    meta = pd.read_csv(meta_path)
    
    # Read imputed training and testing data
    df_imputed = pd.read_csv(imputed_path.with_suffix('.csv'))
    # Process training/testing data: change Pending to Reject and remove Review rows
    df_imputed['Evaluation'] = df_imputed['Evaluation'].replace('Pending', 'Reject')
    df_imputed = df_imputed[df_imputed['Evaluation'] != 'Review']
    df_imputed['data_source'] = "training and testing"
    
    # Get list of external validation patient IDs from metadata
    external_validation_patients = meta[meta['external_validation'] == "Yes"]['patient_id'].tolist()
    
    # Read external validation data
    external_validation_dfs = []
    
    for patient_id in external_validation_patients:
        # Find the imputed file for this patient
        imputed_file_pattern = f"{patient_id}_imputed.csv"
        imputed_files = list(external_validation_dir.glob(imputed_file_pattern))
        
        if not imputed_files:
            print(f"Warning: No imputed file found for patient {patient_id}")
            continue
            
        imputed_file = imputed_files[0]
        
        # Find the corresponding itb_review file
        itb_review_pattern = f"{patient_id}.Annotated.Neoantigen_Candidates.tsv"
        itb_review_files = list(itb_review_dir.glob(itb_review_pattern))
        
        if not itb_review_files:
            print(f"Warning: No itb_review file found for patient {patient_id}")
            continue
            
        itb_review_file = itb_review_files[0]
        
        try:
            # Read imputed data
            df_external_imputed = pd.read_csv(imputed_file)
            
            # Read itb_review data (TSV format)
            df_itb_review = pd.read_csv(itb_review_file, sep='\t')
            
            # Join Evaluation column from itb_review to imputed data using ID
            df_external_imputed = df_external_imputed.merge(
                df_itb_review[['ID', 'Evaluation']], 
                on='ID', 
                how='left', 
                suffixes=('', '_new')
            )
            
            # Replace the old Evaluation column with the new one
            df_external_imputed['Evaluation'] = df_external_imputed['Evaluation_new']
            df_external_imputed = df_external_imputed.drop('Evaluation_new', axis=1)
            
            # Change Pending to Reject in external validation data
            df_external_imputed['Evaluation'] = df_external_imputed['Evaluation'].replace('Pending', 'Reject')
            
            # Add data source column
            df_external_imputed['data_source'] = "external validation"
            
            external_validation_dfs.append(df_external_imputed)
            print(f"Loaded external validation data for patient: {patient_id}")
            
        except Exception as e:
            print(f"Error processing patient {patient_id}: {e}")
    
    # Merge all external validation data
    if external_validation_dfs:
        df_external_merged = pd.concat(external_validation_dfs, ignore_index=True)
        print(f"Loaded {len(external_validation_dfs)} external validation files")
        print(f"External validation data shape: {df_external_merged.shape}")
    else:
        df_external_merged = pd.DataFrame()
        print("No external validation data found")
    
    # Merge training/testing data with external validation data
    if not df_external_merged.empty:
        # Ensure both dataframes have the same columns
        common_columns = list(set(df_imputed.columns) & set(df_external_merged.columns))
        df_imputed_common = df_imputed[common_columns]
        df_external_common = df_external_merged[common_columns]
        
        # Concatenate the dataframes
        df_combined = pd.concat([df_imputed_common, df_external_common], ignore_index=True)
        print(f"Combined data shape: {df_combined.shape}")
        print(f"Training/testing data: {len(df_imputed_common)} rows")
        print(f"External validation data: {len(df_external_common)} rows")
    else:
        df_combined = df_imputed
        print("No external validation data to merge")
    
    return meta, df_combined

def get_trial_summary(meta, df_combined, include_external=True, data_source_filter=None):
    """Calculate summary statistics for trials"""
    # Filter patient IDs based on inclusion criteria
    if include_external:
        condition = (meta['include_ML'] == "Yes") | (meta['external_validation'] == "Yes")
    else:
        condition = meta['include_ML'] == "Yes"
    
    patient_id_all = meta[condition][['patient_id', 'trial', 'tumor_type']]
    
    # Process combined data
    df_trial = df_combined[df_combined['Evaluation'] != "Review"].copy()
    df_trial['Evaluation'] = (df_trial['Evaluation'] == "Accept").astype(int)
    
    # Filter by data source if specified
    if data_source_filter:
        df_trial = df_trial[df_trial['data_source'] == data_source_filter]
    
    # Merge with patient information
    df_trial = df_trial.merge(patient_id_all, on='patient_id', how='inner')
    
    # Calculate summary statistics
    summary = (df_trial.groupby(['trial', 'tumor_type'])
              .agg({
                  'patient_id': 'nunique',
                  'Evaluation': ['count', lambda x: sum(x == 1), lambda x: sum(x == 0)]
              })
              .reset_index())
    
    # Rename columns
    summary.columns = ['trial', 'tumor_type', 'total_patients', 'total_peptides', 
                      'count_accept', 'count_reject']
    
    # Add count_review column
    # For external validation data, count Review evaluations
    if data_source_filter == "external validation" or include_external:
        # Get the original data with Review evaluations
        df_with_review = df_combined.copy()
        if data_source_filter:
            df_with_review = df_with_review[df_with_review['data_source'] == data_source_filter]
        
        df_with_review = df_with_review.merge(patient_id_all, on='patient_id', how='inner')
        
        # Count Review evaluations
        review_counts = (df_with_review[df_with_review['Evaluation'] == 'Review']
                        .groupby(['trial', 'tumor_type'])
                        .size()
                        .reset_index(name='count_review'))
        
        # Merge with summary
        summary = summary.merge(review_counts, on=['trial', 'tumor_type'], how='left')
        summary['count_review'] = summary['count_review'].fillna(0).astype(int)
    else:
        # For training/testing data, set count_review to NA
        summary['count_review'] = 'NA'
    
    # Add totals row
    totals_row = pd.DataFrame({
        'trial': ['TOTAL'],
        'tumor_type': ['ALL'],
        'total_patients': [summary['total_patients'].sum()],
        'total_peptides': [summary['total_peptides'].sum()],
        'count_accept': [summary['count_accept'].sum()],
        'count_reject': [summary['count_reject'].sum()]
    })
    
    # Handle count_review in totals row
    if data_source_filter == "external validation" or include_external:
        if 'count_review' in summary.columns and summary['count_review'].dtype != 'object':
            totals_row['count_review'] = [summary['count_review'].sum()]
        else:
            totals_row['count_review'] = ['NA']
    else:
        totals_row['count_review'] = ['NA']
    
    # Concatenate summary with totals row
    summary_with_totals = pd.concat([summary, totals_row], ignore_index=True)
    
    return summary_with_totals

def main():
    # Create output directory
    outdir = create_output_dir()
    
    # Read data (now includes external validation)
    meta, df_combined = read_data()
    
    # Get summary for all data (including external validation)
    overall_summary = get_trial_summary(meta, df_combined, include_external=True)
    print("\nOverall Summary (Including External Validation):")
    print(overall_summary)
    
    # Get summary for training set only
    training_summary = get_trial_summary(meta, df_combined, include_external=False)
    print("\nTraining Set Summary:")
    print(training_summary)
    
    # Get summary for external validation only
    external_summary = get_trial_summary(meta, df_combined, include_external=True, data_source_filter="external validation")
    print("\nExternal Validation Summary:")
    print(external_summary)
    
    # Save combined dataset
    df_combined.to_csv(outdir / "combined_data_with_external_validation.csv", index=False)
    
    # Export feature list (column names) as a single column
    features_df = pd.DataFrame({'features': df_combined.columns.tolist()})
    features_df.to_csv(outdir / "features_list.csv", index=False)
    print(f"Exported {len(df_combined.columns)} features to features_list.csv")
    
    # Optionally save results
    overall_summary.to_csv(outdir / "overall_summary.csv", index=False)
    training_summary.to_csv(outdir / "training_summary.csv", index=False)
    external_summary.to_csv(outdir / "external_validation_summary.csv", index=False)

if __name__ == "__main__":
    main()