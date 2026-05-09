
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
import sys


zhou_scores = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}


def parse_probs_filename(file_name):
    suffixes = ["_cv_probs.csv", "_probs.csv"]
    suffix = next((s for s in suffixes if file_name.endswith(s)), None)
    if suffix is None:
        return None, None

    stem = file_name[: -len(suffix)]
    parts = stem.split("_")
    if len(parts) < 2:
        return None, None
    feature = parts[-1]

    mapq = None
    for token in parts:
        if token.startswith("mapq") and token[4:].isdigit():
            mapq = token[4:]

    return feature, mapq


def should_use_file(file_name, mapq_filter=None):
    if not (file_name.endswith("_probs.csv") or file_name.endswith("_cv_probs.csv")):
        return False

    feature, mapq = parse_probs_filename(file_name)
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

        feature, file_mapq = parse_probs_filename(file)
        if feature is None:
            continue

        feature_name = feature
        zhou_key = feature

        auc = evaluate_results(os.path.join(svm_dir, file))
        zhou_auc = zhou_scores.get(zhou_key, None)
        results.append({'feature': feature_name, 'auc': auc, 'zhou_auc': zhou_auc})
    if results == []:
        print(f"No probs files found")
        sys.exit(1)
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by='feature', ascending=False).reset_index(drop=True)
    if not os.path.exists("data/model_results"):
        os.makedirs("data/model_results")
    output_path = f"data/model_results/{result_name}_results.csv"
    results_df.to_csv(output_path, index=False)
    print(results_df)
