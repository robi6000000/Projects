
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
import sys
import re


zhou_scores = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}


def extract_feature_and_mapq(file_name):
    """
    Parse filenames like:
      svm_linear_mapq45_edm_cv_probs.csv
      svm_linear_edm_cv_probs.csv
    Returns: (feature, mapq_value_or_None)
    """
    m = re.match(r"^svm_[^_]+(?:_mapq(\d+))?_(.+)_cv_probs\.csv$", file_name)
    if not m:
        return None, None
    mapq = m.group(1)
    feature = m.group(2)
    return feature, mapq


def should_use_file(file_name, mapq_filter=None):
    if not file_name.endswith("_probs.csv"):
        return False

    feature, mapq = extract_feature_and_mapq(file_name)
    if feature is None:
        return False

    if mapq_filter is None:
        return True

    return str(mapq_filter) == str(mapq)

def evaluate_results(file):
    df = pd.read_csv(file)
    # print(df.columns)
    y_true = df['label'].values
    y_pred = df['probability'].values
    auc = roc_auc_score(y_true, y_pred)
    return auc


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("incorrect number of arguments")
        sys.exit(1)
    svm_dir = sys.argv[1].rstrip("/")
    mapq_filter = sys.argv[2] if len(sys.argv) > 2 else None
    parts = svm_dir.split("/")
    result_name = "_".join(parts[-2:]) if len(parts) >= 2 else parts[-1]
    results = []
    for file in sorted(os.listdir(svm_dir)):
        if not should_use_file(file, mapq_filter):
            continue

        feature, file_mapq = extract_feature_and_mapq(file)
        if feature is None:
            continue

        if file_mapq is not None:
            feature_name = f"mapq{file_mapq}_{feature}"
            zhou_key = feature
        else:
            feature_name = feature
            zhou_key = feature

        auc = evaluate_results(os.path.join(svm_dir, file))
        zhou_auc = zhou_scores.get(zhou_key, None)
        results.append({'feature': feature_name, 'auc': auc, 'zhou_auc': zhou_auc})
    if results == []:
        print(f"No probs files found")
        sys.exit(1)
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by='auc', ascending=False).reset_index(drop=True)
    if not os.path.exists("data/model_results"):
        os.makedirs("data/model_results")
    output_path = f"data/model_results/{result_name}_results.csv"
    results_df.to_csv(output_path, index=False)
    print(results_df)
