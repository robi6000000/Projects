
from sklearn.metrics import roc_auc_score
import pandas as pd
import os

svm_dir = "./data/matrix/svm_rbf_data_leak_by_feature/"
zhou_auc = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}

def evaluate_results(file):
    df = pd.read_csv(file)
    print(df.columns)
    y_true = df['label'].values
    y_pred = df['probability'].values
    auc = roc_auc_score(y_true, y_pred)
    return auc

if __name__ == "__main__":
    zhou_auc = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
    }
    results = []
    for file in os.listdir(svm_dir):
        if file.endswith(".csv"):
            feature = file.split("_")[3]
            auc = evaluate_results(os.path.join(svm_dir, file))
            zhou = zhou_auc.get(feature, None)
            results.append({'feature': feature, 'auc': auc, 'zhou_auc': zhou})
    results_df = pd.DataFrame(results)
    print(results_df)
