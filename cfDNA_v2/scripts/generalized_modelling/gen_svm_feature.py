import os
import sys
import gc
import time
import pickle
import pandas as pd
import numpy as np
import pyarrow.csv as pv
from modules.gen_cfdna_model import CFDNAModel

GC_CORRECT_FEATURES = ['pfe', 'coverage', 'ends', 'ocf', 'ifs', 'wps']

''' 
Script for performing Crossvalidation or training an SVM model using the CFDNAModel class

Parameters:
- option: cv/train
- matrix_path: path to input feature matrix (csv)
- parent_matrix_path: path to parent directory of matrix, used for constructing output directories
- gc_path: path to gc content csv 
- feature: name of feature being used, used for output file naming
- kernel: kernel for SVM, (linear/rbf)
- gc_correction: whether to apply gc correction (true/false)
- pca: whether apply PCA (true/false)
- pca_components: how many pca components to use (can be either integer or float from [0,1])
- cv_repeats: number of repeats for cross-validation (integer)
- mapq_filter: mapq score threshold 
- run_tag: optional suffix for output directory and files

Example train:
python scripts/generalized_modelling/gen_svm_feature.py train data/cristiano_runs/feature_matrices/matrix_length_mapq30.csv data/cristiano_runs data/processing/gc_content_per_region.csv length linear true true 150 10 30 none
Example cv:
python scripts/generalized_modelling/gen_svm_feature.py cv data/cristiano_runs/feature_matrices/matrix_length_mapq30.csv data/cristiano_runs data/processing/gc_content_per_region.csv length linear true true 150 10 30 none
'''
if __name__ == "__main__":
    if len(sys.argv) < 12:
        print("Incorrect number of arguments")
        sys.exit(1)

    option = sys.argv[1].lower()
    matrix_path = sys.argv[2]
    parent_matrix_path = sys.argv[3]
    gc_path = sys.argv[4]
    feature = sys.argv[5]
    kernel = sys.argv[6]
    gc_correction = sys.argv[7].lower() == 'true'
    pca = sys.argv[8].lower() == 'true'
    pca_components = float(sys.argv[9])
    cv_repeats = int(sys.argv[10])
    mapq_filter = sys.argv[11]
    if str(mapq_filter).lower() in ['0', 'none']:
        mapq_filter = None
    run_tag = sys.argv[12] if len(sys.argv) >= 13 else None
    if run_tag and run_tag.lower() in ('none', ''):
        run_tag = None

    name = f"svm_{kernel}"
    if pca:
        name += f"_pca{pca_components}"
    if gc_correction:
        name += "_gc"
    if cv_repeats != 10:
        name += f"_cv{cv_repeats}"
    if str(mapq_filter).lower() != 'none':
        name += f"_mapq{mapq_filter}"
    if run_tag:
        name += f"_{run_tag}"

    if gc_correction and feature not in GC_CORRECT_FEATURES:
        print(f"gc_correction=True but '{feature}' is not in GC_CORRECT_FEATURES {GC_CORRECT_FEATURES}. Disabling GC correction.")
        gc_correction = False

    base_output_dir = os.path.join(parent_matrix_path, "model_runs", name)
    cv_dir = os.path.join(base_output_dir, "cv")
    model_dir = os.path.join(base_output_dir, "models")
    os.makedirs(cv_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)

    print(f"Config: option={option}, feature={feature}, kernel={kernel}, gc_correction={gc_correction}, pca={pca}, pca_components={pca_components}, cv_repeats={cv_repeats}, mapq_filter={mapq_filter}")

    print(f"Loading matrix for {feature} (pyarrow.csv, multi-threaded)")
    t0 = time.time()
    table = pv.read_csv(
        matrix_path,
        read_options=pv.ReadOptions(use_threads=True, block_size=256 * 1024 * 1024),
    )
    print(f"pyarrow.read_csv finished {time.time() - t0:.1f}s, table shape={table.num_rows}x{table.num_columns}")

    t1 = time.time()
    mx = table.to_pandas(self_destruct=True)
    del table
    gc.collect()
    print(f"to_pandas done in {time.time() - t1:.1f}s")

    mx = mx.set_index(mx.columns[0])

    mem_gb = mx.memory_usage(deep=True).sum() / (1024 ** 3)
    print(f"Matrix loaded: shape={mx.shape}, approx_memory={mem_gb:.2f} GB")

    gc_content = pd.read_csv(gc_path)

    model = CFDNAModel(mx, gc_content, feature=feature,
                       kernel=kernel, gc_correction=gc_correction,
                       pca=pca, pca_components=pca_components, cv_repeats=cv_repeats)
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
        out_path = os.path.join(cv_dir, f"{feature}_cv_probs.csv")
        output.to_csv(out_path, index=False)
        print(f"CV results saved {out_path}")
    elif option == 'train':
        print(f"Training model for {feature}")
        trained_model = model.train_model()
        model_path = os.path.join(model_dir, f"{name}_{feature}.pkl")
        with open(model_path, 'wb') as f:
            pickle.dump(trained_model, f)
        print(f"Trained model saved {model_path}")
    else:
        raise ValueError(f"Unknown option - Use cv/train")
