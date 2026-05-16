import os
import sys
import gc
import time
import pickle
import pandas as pd
import numpy as np
import pyarrow.csv as pv
from modules.gen_cfdna_model import CFDNAModel

'''
Script for applying a trained SVM model to a test feature matrix using the CFDNAModel class

Parameters:
- test_matrix_path: path to test feature matrix (csv)
- model_path: path to trained model pickle (.pkl)
- output_probs_path: path to output prediction csv
- gc_path: path to gc content csv (required if the trained model uses gc_correction)

Example test:
python scripts/generalized_modelling/gen_svm_test.py data/internal_ly_test_runs/feature_matrices/matrix_length_mapq30.csv data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/models/svm_linear_pca150.0_gc_mapq30_length.pkl data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/test/prevelynch_test/length_probs.csv data/processing/gc_content_per_region.csv
'''
if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        print("Usage: python scripts/generalized_modelling/gen_svm_test.py test_matrix_path model_path output_probs_path [gc_path?]")
        sys.exit(1)

    matrix_path = sys.argv[1]
    model_path = sys.argv[2]
    output_path = sys.argv[3]
    gc_path = sys.argv[4] if len(sys.argv) == 5 else None

    print(f"Predict config (v2): matrix_path={matrix_path}, model_path={model_path}, output_path={output_path}, gc_path={gc_path}")

    with open(model_path, 'rb') as f:
        trained_model = pickle.load(f)

    gc_content = None
    if trained_model.get('gc_correction'):
        if not gc_path:
            raise ValueError(
                "Trained model uses gc_correction=True but no gc_path was provided. "
                "Pass the GC content CSV as the 4th argument."
            )
        gc_content = pd.read_csv(gc_path)

    print("Loading test matrix with pyarrow.csv (multi-threaded)")
    t0 = time.time()
    table = pv.read_csv(
        matrix_path,
        # using larger block size because of 1.6m col matrix
        read_options=pv.ReadOptions(use_threads=True, block_size=256 * 1024 * 1024),
    )
    print(f"pyarrow.read_csv done in {time.time() - t0:.1f}s, table shape({table.num_rows},{table.num_columns})")

    mx = table.to_pandas(self_destruct=True)
    del table
    gc.collect()

    mx = mx.set_index(mx.columns[0])

    metadata_cols = ['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true']
    model_feature_cols = trained_model.get('feature_columns')
    if model_feature_cols:
        # need to test the runtime because with
        model_feature_cols = list(model_feature_cols)
        col_set = set(mx.columns)
        missing_features = [c for c in model_feature_cols if c not in col_set]
        if missing_features:
            raise ValueError(f"Test matrix missing {len(missing_features)} required feature columns.")

        present_metadata = [c for c in metadata_cols if c in col_set]
        mx = mx[present_metadata + model_feature_cols]

    print(f"Test matrix shape: {mx.shape}")

    model = CFDNAModel(
        mx,
        gc_content=gc_content,
        feature=trained_model.get('feature'),
        kernel=trained_model.get('kernel'),
        gc_correction=bool(trained_model.get('gc_correction', False)),
        pca=bool(trained_model.get('pca', False)),
        pca_components=trained_model.get('pca_components', 0.95),
        cv_repeats=1,
    )

    del mx
    gc.collect()

    print("Running predict (scaler + GC + PCA + SVC)")
    t3 = time.time()
    output = model.predict(trained_model)

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    output.to_csv(output_path, index=False)
    print(f"Test predictions saved {output_path}")
