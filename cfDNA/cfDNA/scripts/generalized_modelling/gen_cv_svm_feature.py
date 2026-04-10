import os
import sys
import gc
import pickle
import pandas as pd
import numpy as np
from modules.gen_cfdna_model import CFDNAModel

GC_CORRECT_FEATURES = ['pfe', 'fsr', 'coverage', 'ends', 'ocf', 'ifs', 'wps']

if __name__ == "__main__":
    if len(sys.argv) < 9:
        print(
            "Usage: python scripts/generalized_modelling/gen_cv_svm_feature.py "
            "[option:cv/train/test] [matrix_path] [parent_matrix_path] [gc_path] [feature] [kernel] "
            "[gc_correction: bool] [pca: bool] [pca_components] [cv_repeats] [mapq_filter] [model_path] [test_output_path]"
        )
        sys.exit(1)

    option = sys.argv[1].lower() # cv, train, or test
    matrix_path = sys.argv[2]
    # if matrix is in data/matrix_internal_gs/by_feature/, parent_path is data/matrix_internal_gs
    parent_matrix_path = sys.argv[3]
    gc_path = sys.argv[4]
    feature = sys.argv[5]
    kernel = sys.argv[6]
    gc_correction = sys.argv[7].lower() == 'true'
    pca = sys.argv[8].lower() == 'true'
    pca_components = float(sys.argv[9]) if len(sys.argv) >= 10 else 0.95
    cv_repeats = int(sys.argv[10]) if len(sys.argv) >= 11 else 10
    mapq_filter = sys.argv[11] if len(sys.argv) >= 12 else 'None'

    name = f"svm_{kernel}"
    if pca:
        name += f"_pca{pca_components}"
    if gc_correction:
        name += "_gc"
    if cv_repeats != 10:
        name += f"_cv{cv_repeats}"
    if str(mapq_filter).lower() != 'none':
        name += f"_mapq{mapq_filter}"

    base_output_dir = os.path.join(parent_matrix_path, "svm_by_feature", name)
    cv_dir = os.path.join(base_output_dir, "cv")
    model_dir = os.path.join(base_output_dir, "models")
    test_dir = os.path.join(base_output_dir, "test")
    os.makedirs(cv_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(test_dir, exist_ok=True)

    model_path = sys.argv[12] if len(sys.argv) >= 13 else os.path.join(model_dir, f"{name}_{feature}.pkl")
    test_output_path = sys.argv[13] if len(sys.argv) >= 14 else os.path.join(test_dir, f"{name}_{feature}_probs.csv")


    print(f"Config: option={option}, feature={feature}, kernel={kernel}, gc_correction={gc_correction}, pca={pca}, pca_components={pca_components}, cv_repeats={cv_repeats}, mapq_filter={mapq_filter}")

    print(f"Loading matrix for {feature}")
    mx = pd.read_csv(matrix_path, index_col=0, low_memory=False)
    metadata_cols = ['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true']
    numeric_cols = [c for c in mx.columns if c not in metadata_cols]
    mx[numeric_cols] = mx[numeric_cols].apply(pd.to_numeric, errors='coerce').astype(np.float32)

    mem_gb = mx.memory_usage(deep=True).sum() / (1024 ** 3)
    print(f"Matrix loaded: shape={mx.shape}, approx_memory={mem_gb:.2f} GB")

    gc_content = pd.read_csv(gc_path)

    model = CFDNAModel(mx, gc_content, feature=feature,
                       kernel=kernel, gc_correction=gc_correction, 
                       pca=pca, pca_components=pca_components, cv_repeats=cv_repeats)

    metadata = mx[[c for c in metadata_cols if c in mx.columns]].copy()
    # delete large dataframe to free memory
    del mx
    gc.collect()

    if option == 'cv':
        print(f"Cross-validating SVM for {feature}")
        results = model.cv_svm()
        output = pd.DataFrame({
            'sample_id': results['sample_ids'],
            'probability': results['probabilities'],
            'label': results['labels']
        })
        out_path = os.path.join(cv_dir, f"{name}_{feature}_cv_probs.csv")
        output.to_csv(out_path, index=False)
        print(f"CV results saved {out_path}")
    elif option == 'train':
        print(f"Training model for {feature}")
        trained_model = model.train_model()
        with open(model_path, 'wb') as f:
            pickle.dump(trained_model, f)
        print(f"Trained model saved {model_path}")
    elif option == 'test':
        print(f"Testing model for {feature}")
        with open(model_path, 'rb') as f:
            trained_model = pickle.load(f)
        try:
            output = model.test_model(trained_model)
            if isinstance(output, pd.DataFrame):
                output.to_csv(test_output_path, index=False)
                print(f"Test predictions saved {test_output_path}")
            else:
                print("Test method placeholder returned no DataFrame yet.")
        except NotImplementedError as e:
            print(str(e))
            print("No test output written yet.")
    else:
        raise ValueError(f"Unknown option: {option}. Use cv, train, or test.")