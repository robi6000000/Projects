import pandas as pd
import subprocess
import sys
import os

# setup
manifest_path = "./data/manifest/Cristiano_manifest.csv"
samples_folder = "./data/source/cristiano"
rows_folder = "./data/rows_sample_temp"
matrix_folder = "./data/matrix"
sample_prep_script = "./scripts/sample_prep_script.sh"

os.makedirs(rows_folder, exist_ok=True)
os.makedirs(matrix_folder, exist_ok=True)

# Load manifest
manifest = pd.read_csv(manifest_path)
print(f"number of samples in manifest: {len(manifest)}")
print()

processed = 0
failed = 0

# processing each sample
for idx, row in manifest.iterrows():
    seqrun_id = row['seqrun_id']
    sample_name = row['sample_name']
    stage = row['stage']
    
    sample_id = f"EE{seqrun_id}"
    frag_file = os.path.join(samples_folder, f"{sample_id}.hg19.frag.tsv.bgz")
    feature_file = os.path.join(rows_folder, f"{sample_id}_features.csv")
    
    # checks
    if os.path.exists(feature_file):
        print(f"{sample_id} already processed")
        continue
    if not os.path.exists(frag_file):
        print(f"{sample_id} missing")
        continue
    
    # create features foer each sample
    print(f"Processing {sample_id} ({processed + 1}/{len(manifest)}, stage: {stage})...")
    
    try:
        result = subprocess.run(["bash", sample_prep_script, frag_file, sample_id])
        processed += 1
    except subprocess.CalledProcessError as e:
        failed += 1
        print(f"stderr: {e.stderr}")
    
    print()
    
    # # test
    # if processed >= 5:
    #     break


print(f"Total in manifest: {len(manifest)}")
print(f"Processed: {processed}")
print(f"Failed: {failed}")
print()

# matrix construciton
# load rows into a list of dataframes
rows = []
for file in os.listdir(rows_folder):
    if file.endswith("_features.csv"):
        df_row = pd.read_csv(os.path.join(rows_folder, file))
        rows.append(df_row)
matrix = pd.concat(rows, ignore_index=True) 

# add stage, disease and tissue metadata
matrix['seqrun_id'] = matrix['sample_id'].str.replace("EE", "")
matrix = matrix.merge(manifest[['seqrun_id','sample_name', 'stage', 'disease', 'tissue']], on='seqrun_id', how='left')

metadata_cols = ['sample_id', 'seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']
feature_cols = [col for col in matrix.columns if col not in metadata_cols]
matrix = matrix[metadata_cols + feature_cols]   #just to reorder for testing

matrix.to_csv(os.path.join(matrix_folder, "feature_matrix.csv"), index=False)
print(f"Matrix shape: {matrix.shape}")