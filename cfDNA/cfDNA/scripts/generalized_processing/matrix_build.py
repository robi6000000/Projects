import pandas as pd
import numpy as np
import sys
import os
import glob

features_folder = "./data/internal_features"
matrix_folder = "./data/matrix/internal_by_feature"
metadata_path = "./data/manifest/internal_metadata.csv"


def build_feature_matrix(feature_name, mapq_filter=None, features_folder=features_folder, matrix_folder=matrix_folder, metadata_path=metadata_path):
    """
    Build feature matrix for a single fragmentation pattern.
    """

    if mapq_filter is not None:
        feature_dir = f"{features_folder}/{feature_name}_mapq{mapq_filter}"
        output_name = f"{feature_name}_mapq{mapq_filter}"
    else:
        feature_dir = f"{features_folder}/{feature_name}"
        output_name = feature_name

    print(f"Building matrix for: {output_name}")
    
    # get all sample files in the given folder
    feature_files = sorted(glob.glob(f"{feature_dir}/*.csv"))
    print(f"Found {len(feature_files)} files")
    
    if len(feature_files) == 0:
        print(f"No files - {output_name}")
        return
    
    metadata = pd.read_csv(metadata_path)
    valid_ids = set(metadata['sample_id'].astype(str))

    feature_files = [
        f for f in feature_files
        if os.path.basename(f).replace(f"_{output_name}.csv", "") in valid_ids
    ]
    print(f"After metadata filter: {len(feature_files)} files")

    if len(feature_files) == 0:
        print(f"No files after filtering - {output_name}")
        return
    # load first file to get structure (matrix vs vector feature)
    first_file = feature_files[0]
    
    df = pd.read_csv(first_file, index_col=0)
    print(f"Feature shape: {df.shape}")
    if len(df.columns) == 1:
        column_names = [f"{output_name}_region_{i}" for i in df.index]
    else:
        column_names = []
        for row in df.index:
            for col_name in df.columns:
                # fixing inconsistent colnames having spaces
                clean_col = str(col_name).replace(' ', '')
                column_names.append(f"{output_name}_{row}_bin_{clean_col}")
    
    print(f"Total features: {len(column_names)}")
    print("Load samples")
    rows = []
    sample_ids = []
    
    for i, file_path in enumerate(feature_files):
        if i % 10 == 0:
            print(f"  {i+1}/{len(feature_files)}")
        
        sid = os.path.basename(file_path).replace(f"_{output_name}.csv", "")
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
    matrix = matrix.merge(
        metadata[['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true']],
        on='sample_id',
        how='left'
    )

    # Reorder to metadata first
    metadata_cols = ['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true']
    feature_cols = [c for c in matrix.columns if c not in metadata_cols]
    matrix = matrix[metadata_cols + feature_cols]
    
    # Save
    os.makedirs(matrix_folder, exist_ok=True)
    output_path = f"{matrix_folder}/matrix_{output_name}.csv"
    print(f"Saving to {output_path}")
    matrix.to_csv(output_path, index=False)
    print(f"Saved: {matrix.shape}")
    print(f"File size: {os.path.getsize(output_path) / (1024**2):.1f} MB")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Args: feature mapq_filter features_folder matrix_folder
        feature = sys.argv[1]
        mapq_filter = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2].lower() not in ["none", "0"] else None
        features_folder = sys.argv[3] if len(sys.argv) > 3 else features_folder
        matrix_folder = sys.argv[4] if len(sys.argv) > 4 else matrix_folder
        metadata_path = sys.argv[5] if len(sys.argv) > 5 else metadata_path

        build_feature_matrix(feature, mapq_filter, features_folder, matrix_folder, metadata_path)
    else:
        features = ['length',
                    'pfe', 'fsr', 'fsd', 'coverage',
                    'ends',
                    'ocf',
                    'ifs', 'wps', 'edm',
                    'iedm', 'eedm'
                    ]
        mapq_filters = [30, 20, 10, None]

        for feature in features:
            for mapq_filter in mapq_filters:
                try:
                    build_feature_matrix(feature, mapq_filter, features_folder, matrix_folder, metadata_path)
                except Exception as e:
                    print(f"Error building {feature} with MAPQ filter {mapq_filter}: {e}")
                    continue

        print("Finished all features")