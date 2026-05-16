import pandas as pd
import numpy as np
import sys
import os
import glob

os.chdir("/data/projects/liquid_biopsy/Projects/cfDNA_v2")
probs_folder = "./data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_cv1_mapq30/cv/"
output_path = "./data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_cv1_mapq30/meta_matrix/meta_matrix.csv"
metadata_path = "./data/metadata/Cristiano_metadata.csv"


def build_meta_matrix(probs_folder, output_path, metadata_path):
    prob_files = sorted(glob.glob(f"{probs_folder}/*.csv"))
    print(f"Found {len(prob_files)} probability files")
    probs_df = pd.DataFrame()
    for file in prob_files:
        df = pd.read_csv(file, index_col=0)
        feature_name = os.path.basename(file).split("_")[0]
        probs_df[feature_name] = df["probability"]
    metadata = pd.read_csv(metadata_path)
    metadata.set_index("sample_id", inplace=True)
    metadata_cols = ["disease", "dataset", "material", "stage", "cancer_true"]
    meta_matrix = metadata[metadata_cols].join(probs_df, how="inner")
    print("Meta-matrix shape after join:", meta_matrix.shape)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    meta_matrix.to_csv(output_path)
if __name__ == "__main__":
    if len(sys.argv) > 1:
        # example: python meta_matrix_build.py probs_folder output_path metadata_path
        probs_folder = sys.argv[1]
        output_path = sys.argv[2]
        metadata_path = sys.argv[3]

        build_meta_matrix(probs_folder, output_path, metadata_path)
    else:
        
        build_meta_matrix(probs_folder, output_path, metadata_path)
        print("Finished all features")