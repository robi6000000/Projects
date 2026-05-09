import os
import sys
import gc
import pickle
import pandas as pd
import numpy as np
from modules.gen_cfdna_model import CFDNAModel

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python scripts/generalized_modelling/gen_svm_test.py "
            "[test_matrix_path] [model_path] [output_probs_path]"
        )
        sys.exit(1)

    matrix_path = sys.argv[1]
    model_path = sys.argv[2]
    output_path = sys.argv[3]

    print(f"Predict config: matrix_path={matrix_path}, model_path={model_path}, output_path={output_path}")

    with open(model_path, 'rb') as f:
        trained_model = pickle.load(f)

    print("Loading test matrix")
    mx = pd.read_csv(matrix_path, index_col=0, low_memory=False)

    metadata_cols = ['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true']
    numeric_cols = [c for c in mx.columns if c not in metadata_cols]
    mx[numeric_cols] = mx[numeric_cols].apply(pd.to_numeric, errors='coerce').astype(np.float32)

    model_feature_cols = trained_model.get('feature_columns')
    if model_feature_cols:
        model_feature_cols = list(model_feature_cols)
        missing_features = [c for c in model_feature_cols if c not in mx.columns]
        if missing_features:
            raise ValueError(
                f"Test matrix is missing {len(missing_features)} required feature columns. "
                f"First missing columns: {missing_features[:10]}"
            )

        present_metadata = [c for c in metadata_cols if c in mx.columns]
        mx = mx[present_metadata + model_feature_cols]

    print(f"Test matrix loaded: shape={mx.shape}")

    model = CFDNAModel(
        mx,
        gc_content=None,
        feature=trained_model.get('feature'),
        kernel=trained_model.get('kernel'),
        gc_correction=bool(trained_model.get('gc_correction', False)),
        pca=bool(trained_model.get('pca', False)),
        pca_components=trained_model.get('pca_components', 0.95),
        cv_repeats=1,
    )

    del mx
    gc.collect()

    output = model.predict(trained_model)

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    output.to_csv(output_path, index=False)
    print(f"Test predictions saved {output_path}")
