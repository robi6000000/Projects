import pandas as pd
import numpy as np
import sys
import os
import glob

features_folder = "./data/cristiano_features"
matrix_folder = "./data/matrix/by_feature"
manifest_path = "./data/manifest/Cristiano_manifest.csv"

os.makedirs(matrix_folder, exist_ok=True)

def build_feature_matrix(feature_name):
    """
    Build feature matrix for a single fragmentation pattern.
    """
    
    feature_dir = f"{features_folder}/{feature_name}"
    
    print(f"Building matrix for: {feature_name}")
    
    # Get all sample files
    feature_files = sorted(glob.glob(f"{feature_dir}/*.csv"))
    print(f"Found {len(feature_files)} files")
    
    if len(feature_files) == 0:
        print(f"No files - {feature_name}")
        return
    
    # Load first file to get structure
    first_file = feature_files[0]
    
    df = pd.read_csv(first_file, index_col=0)
    print(f"Feature shape: {df.shape}")
    if len(df.columns) == 1:
        column_names = [f"{feature_name}_region_{i}" for i in df.index]
    else:
        column_names = []
        for row in df.index:
            for col_name in df.columns:
                clean_col = str(col_name).replace(' ', '')
                column_names.append(f"{feature_name}_{row}_bin_{clean_col}")
    
    print(f"Total features: {len(column_names)}")
    print("Load samples")
    rows = []
    sample_ids = []
    
    for i, file_path in enumerate(feature_files):
        if i % 10 == 0:
            print(f"  {i+1}/{len(feature_files)}")
        
        sid = os.path.basename(file_path).replace(f"_{feature_name}.csv", "")
        df = pd.read_csv(file_path, index_col=0)
        
        # Flatten
        values = df.values.flatten()
        if i % 50 == 0:
            print(f"    Sample {sid} - values: {values[:30]}")
        rows.append(values)
        sample_ids.append(sid)
    
    matrix = pd.DataFrame(rows, columns=column_names)
    matrix.insert(0, 'sample_id', sample_ids)
    print(f"Matrix done: {matrix.shape}")
    
    # Add metadata
    print("Metadata")
    manifest = pd.read_csv(manifest_path)
    matrix['seqrun_id'] = matrix['sample_id'].str.replace('EE', '').astype(int)
    matrix = matrix.merge(
        manifest[['seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']],
        on='seqrun_id',
        how='left'
    )
    
    # Reorder to metadata first
    metadata_cols = ['sample_id', 'seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']
    feature_cols = [c for c in matrix.columns if c not in metadata_cols]
    matrix = matrix[metadata_cols + feature_cols]
    
    # Save
    output_path = f"{matrix_folder}/matrix_{feature_name}.csv"
    print(f"Saving to {output_path}")
    matrix.to_csv(output_path, index=False)
    print(f"Saved: {matrix.shape}")
    print(f"File size: {os.path.getsize(output_path) / (1024**2):.1f} MB")

if __name__ == "__main__":    
    if len(sys.argv) > 1:
        # Run specific feature from command line
        feature = sys.argv[1]
        build_feature_matrix(feature)
    else:
        features = ['length', 
                    # 'pfe', 'fsr', 'fsd', 'coverage', 
                    'ends', 
                    # 'ocf', 
                    # 'ifs', 'wps', 'edm'
                    ]
        
        for feature in features:
            try:
                build_feature_matrix(feature)
            except Exception as e:
                print(f"Error building {feature}: {e}")
                continue
        
        print("Finished all features")