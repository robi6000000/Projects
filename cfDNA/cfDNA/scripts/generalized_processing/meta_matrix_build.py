import pandas as pd
import numpy as np
import sys
import os
import glob

# set current path in os
os.chdir("/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA")
probs_folder = "./data/matrix/svm_by_feature/svm_linear_pca150.0_gc_cv1_mapq30/cv/"
output_folder = "./data/matrix/svm_by_feature/svm_linear_pca150.0_gc_cv1_mapq30/meta_matrix/"
os.makedirs(output_folder, exist_ok=True)
metadata_path = "./data/manifest/Cristiano_metadata.csv"


def build_meta_matrix(probs_folder, output_folder, metadata_path):
    prob_files = sorted(glob.glob(f"{probs_folder}/*.csv"))
    print(f"Found {len(prob_files)} probability files")
    probs_df = pd.DataFrame()
    for file in prob_files:
        df = pd.read_csv(file, index_col=0)
        feature_name = os.path.basename(file).split("_")[-3]
        probs_df[feature_name] = df["probability"]
    print("Combined probabilities shape:", probs_df.shape)
    metadata = pd.read_csv(metadata_path)
    metadata.set_index("sample_id", inplace=True)
    metadata_cols = ["disease", "dataset", "material", "stage", "cancer_true"]
    meta_matrix = metadata[metadata_cols].join(probs_df, how="inner")
    print("Meta-matrix shape after join:", meta_matrix.shape)
    meta_matrix.to_csv(f"{output_folder}/meta_matrix.csv")
if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Args: probs_folder output_folder metadata_path
        probs_folder = sys.argv[1]
        output_folder = sys.argv[2]
        metadata_path = sys.argv[3]
        os.makedirs(output_folder, exist_ok=True)

        build_meta_matrix(probs_folder, output_folder, metadata_path)
    else:
        
        build_meta_matrix(probs_folder, output_folder, metadata_path)

        print("Finished all features")