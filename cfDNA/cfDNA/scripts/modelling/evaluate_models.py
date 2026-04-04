
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
import sys


zhou_scores = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}


def should_use_file(file_name, mapq_filter=None):
    """
    check if the file does or does not contain mapq filtering
    """
    if not file_name.endswith(".csv"):
        return False

    stem = file_name[:-4]
    if mapq_filter is None:
        return "_mapq" not in stem

    mapq_token = f"_mapq{mapq_filter}"
    return mapq_token in stem

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

        # specific feature naming containing _
        if "ext_poem_prem" in file:
            feature = "ext_poem_prem" if mapq_filter is None else f"mapq{mapq_filter}_ext_poem_prem"
        elif "poem_prem" in file:
            feature = "poem_prem" if mapq_filter is None else f"mapq{mapq_filter}_poem_prem"
        elif "wps_compute" in file:
            feature = "wps_compute" if mapq_filter is None else f"mapq{mapq_filter}_wps_compute"
        else:
            feature = file.split("_")[-2] if mapq_filter is None else f"mapq{mapq_filter}_{file.split('_')[-2]}"
        auc = evaluate_results(os.path.join(svm_dir, file))
        zhou_auc = zhou_scores.get("_".join(feature.split("_")[1:]), None)
        results.append({'feature': feature, 'auc': auc, 'zhou_auc': zhou_auc})
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
