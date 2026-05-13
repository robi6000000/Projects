import pandas as pd
import numpy as np
import sys
import os
import glob

EXCLUDE_FEATURES = {'fsr'}

# set current path in os
os.chdir("/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA")
probs_folder = "./data/matrix/svm_by_feature/svm_linear_pca150.0_gc_cv1_mapq30/cv/"
output_path = "./data/matrix/svm_by_feature/svm_linear_pca150.0_gc_cv1_mapq30/meta_matrix/meta_matrix_nofsr.csv"
metadata_path = "./data/manifest/Cristiano_metadata.csv"


def build_meta_matrix(probs_folder, output_path, metadata_path):
    prob_files = sorted(glob.glob(f"{probs_folder}/*.csv"))
    print(f"Found {len(prob_files)} probability files")
    probs_df = pd.DataFrame()
    for file in prob_files:
        feature_name = os.path.basename(file).split("_")[0]
        if feature_name in EXCLUDE_FEATURES:
            print(f"Skipping excluded feature: {feature_name}")
            continue
        df = pd.read_csv(file, index_col=0)
        probs_df[feature_name] = df["probability"]
    print("Combined probabilities shape:", probs_df.shape)
    metadata = pd.read_csv(metadata_path)
    metadata.set_index("sample_id", inplace=True)
    metadata_cols = ["disease", "dataset", "material", "stage", "cancer_true"]
    meta_matrix = metadata[metadata_cols].join(probs_df, how="inner")
    print("Meta-matrix shape after join:", meta_matrix.shape)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    meta_matrix.to_csv(output_path)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Args: probs_folder output_path metadata_path
        probs_folder = sys.argv[1]
        output_path = sys.argv[2]
        metadata_path = sys.argv[3]

        build_meta_matrix(probs_folder, output_path, metadata_path)
    else:

        build_meta_matrix(probs_folder, output_path, metadata_path)

        print("Finished all features")
