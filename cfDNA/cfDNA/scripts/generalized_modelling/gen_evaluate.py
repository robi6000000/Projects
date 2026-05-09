import os
import sys
import pandas as pd
from sklearn.metrics import roc_auc_score

# Zhou et al. (2022) reference AUCs — remove before final thesis figures.
ZHOU_SCORES = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}

def compute_auc(filepath):
    df = pd.read_csv(filepath)
    return roc_auc_score(df['label'].values, df['probability'].values)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gen_evaluate.py [directory] [output_csv]")
        sys.exit(1)

    input_dir = sys.argv[1]
    out_filepath = sys.argv[2]

    results = []
    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith('.csv'):
            continue
        filepath = os.path.join(input_dir, filename)
        feature = os.path.basename(filepath).split('_')[0]
        auc = compute_auc(filepath)
        results.append({'feature': feature, 'auc': auc, 'zhou_auc': ZHOU_SCORES.get(feature)})

    df = pd.DataFrame(results)
    if os.path.dirname(out_filepath):
        os.makedirs(os.path.dirname(out_filepath), exist_ok=True)
    df.to_csv(out_filepath, index=False)
    print(df.to_string(index=False))
    print(f"\nSaved to {out_filepath}")
