import os
import sys
import gc
import pandas as pd
import numpy as np
from modules.cfdna_model import CFDNAModel

GC_CORRECT_FEATURES = ['pfe', 'fsr', 'coverage', 'ends', 'ocf', 'ifs', 'wps']

if __name__ == "__main__":
    if len(sys.argv) < 8:
        print(
            "Usage: python scripts/modelling/svm_feature.py "
            "<matrix_path> <gc_path> <feature> <kernel> <gc_correction:true|false> "
            "<pca:true|false> [pca_components] [cv_repeats] [mapq_filter]"
        )
        sys.exit(1)

    matrix_path = sys.argv[1]
    gc_path = sys.argv[2]
    feature = sys.argv[3]
    kernel = sys.argv[4]
    gc_correction = sys.argv[5].lower() == 'true'
    pca = sys.argv[6].lower() == 'true'
    pca_components = float(sys.argv[7]) if len(sys.argv) >= 8 else 0.95
    cv_repeats = int(sys.argv[8]) if len(sys.argv) >= 9 else 10
    mapq_filter = sys.argv[9] if len(sys.argv) == 10 else 'None'

    print(f"Config: feature={feature}, kernel={kernel}, gc_correction={gc_correction}, pca={pca}, pca_components={pca_components}, cv_repeats={cv_repeats}, mapq_filter={mapq_filter}")

    print(f"Loading matrix for {feature}")
    mx = pd.read_csv(matrix_path, index_col=0, low_memory=False)
    metadata_cols = ['seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']
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

    print(f"Cross-validating SVM for {feature}")
    results = model.cv_svm()

    output = pd.DataFrame({
        'sample_id': results['sample_ids'],
        'probability': results['probabilities'],
        'label': results['labels']
    })

    metadata_cols = ['sample_name', 'stage', 'disease', 'tissue']
    for col in metadata_cols:
        if col in metadata.columns:
            output[col] = metadata[col].values
    output = output[['sample_id'] + metadata_cols + ['label', 'probability']]

    # adjust specific folder name
    subfolder = f"svm_{kernel}"
    if pca:
        subfolder += f"_pca{pca_components}"
    if gc_correction:
        subfolder += "_gc"
    if cv_repeats != 10:
        subfolder += f"_cv{cv_repeats}"
    if str(mapq_filter).lower() != 'none':
        subfolder += f"_mapq{mapq_filter}"

    output_folder = f"./data/matrix/svm_by_feature/{subfolder}"
    os.makedirs(output_folder, exist_ok=True)

    out_path = f"{output_folder}/{subfolder}_{feature}_probs.csv"
    
    output.to_csv(out_path, index=False)
    print(f"Probabilities saved {out_path}")