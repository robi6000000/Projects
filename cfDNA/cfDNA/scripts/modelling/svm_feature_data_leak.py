import sys
import pandas as pd
import numpy as np
from modules.cfdna_model_data_leak import CFDNAModel

GC_CORRECT_FEATURES = ['pfe', 'fsr', 'coverage', 'ends', 'ocf', 'ifs', 'wps']

if __name__ == "__main__":
    if len(sys.argv) != 7:
        # format of input arguments: matrix_path, gc_content_path, feature_name, kernel, gc_correction (True/False), pca (True/False)
        print("Incorrect number of arguments.")
        sys.exit(1)

    matrix_path = sys.argv[1]
    gc_path = sys.argv[2]
    feature = sys.argv[3]
    kernel = sys.argv[4]
    gc_correction = sys.argv[5].lower() == 'true'
    pca = sys.argv[6].lower() == 'true'


    print(f"Loading matrix for {feature}")
    mx = pd.read_csv(matrix_path, index_col=0, low_memory=False)
    gc_content = pd.read_csv(gc_path)

    gc_correction = feature in GC_CORRECT_FEATURES
    model = CFDNAModel(mx, gc_content, feature=feature,
                       kernel=kernel, gc_correction=gc_correction, pca=pca)

    print(f"Training SVM for {feature}")
    results = model.train_svm()

    output = pd.DataFrame({
        'sample_id': results['sample_ids'],
        'probability': results['probabilities'],
        'label': results['labels']
    })

    metadata_cols = ['sample_name', 'stage', 'disease', 'tissue']
    for col in metadata_cols:
        if col in mx.columns:
            output[col] = mx[col].values
    output = output[['sample_id'] + metadata_cols + ['label', 'probability']]

    # out_path = f"./data/matrix/svm_by_feature/svm_{feature}_probs.csv"
    out_path = f"./data/matrix/svm_rbf_data_leak_by_feature/svm_rbf_gc_{feature}_probs.csv"
    output.to_csv(out_path, index=False)
    print(f"probabilities saved {out_path}")