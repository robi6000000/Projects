import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score

FEATURE_NAMES = ['length', 'fsd', 'fsr', 'pfe', 'coverage', 'ends', 'ocf',
                 'ifs', 'wps', 'edm', 'iedm', 'eedm', 'eoedm', 'cposedm']

RENAME     = {'cposedm': 'ocedm'}
RENAME_INV = {v: k for k, v in RENAME.items()}

ZHOU_SCORES = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}

NEW_FEATURES = {'iedm', 'eedm', 'eoedm', 'ocedm'}


def load_roc_data(directory, feature_names=None, suffix='_probs.csv'):
    """Load per-feature ROC data. Returns list of dicts with name/fpr/tpr/auc, sorted by AUC desc."""
    if feature_names is None:
        feature_names = FEATURE_NAMES
    roc_data, missing = [], []
    for feature in feature_names:
        path = os.path.join(directory, f"{feature}{suffix}")
        if not os.path.exists(path):
            missing.append(feature)
            continue
        df = pd.read_csv(path)[['probability', 'label']]
        fpr, tpr, _ = roc_curve(df['label'], df['probability'])
        auc = roc_auc_score(df['label'], df['probability'])
        roc_data.append({'name': RENAME.get(feature, feature), 'fpr': fpr, 'tpr': tpr, 'auc': auc})
    roc_data.sort(key=lambda x: x['auc'], reverse=True)
    if missing:
        print(f"Missing: {missing}")
    return roc_data


def load_probs(directory, feature_names=None, suffix='_cv_probs.csv'):
    """Load raw probability DataFrames per feature. Returns dict: feature -> DataFrame."""
    if feature_names is None:
        feature_names = FEATURE_NAMES
    probs = {}
    for feature in feature_names:
        path = os.path.join(directory, f"{feature}{suffix}")
        if not os.path.exists(path):
            continue
        probs[feature] = pd.read_csv(path)
    return probs


def bootstrap_auc_ci(y_true, y_prob, n=1000, seed=42):
    """Bootstrap 95% CI for AUC. Returns (lower, upper)."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_prob[idx]))
    return tuple(np.percentile(aucs, [2.5, 97.5]))


def sens_at_spec(fpr, tpr, target_spec):
    """Sensitivity (TPR) at a given specificity. target_spec e.g. 0.95."""
    mask = fpr <= (1.0 - target_spec)
    return float(tpr[mask].max()) if mask.any() else 0.0


def build_auc_table(roc_dict):
    """
    Build AUC summary DataFrame from {condition_label: roc_data_list}.
    Returns DataFrame: rows=features, columns=conditions.
    """
    records = {cond: {d['name']: d['auc'] for d in roc_data}
               for cond, roc_data in roc_dict.items()}
    df = pd.DataFrame(records)
    df.index.name = 'feature'
    return df


def build_sensitivity_table(roc_dict, specificities=(0.95, 0.85)):
    """
    Build sensitivity-at-fixed-specificity table from {condition_label: roc_data_list}.
    Returns DataFrame: rows=features, columns=condition_sens{spec_pct}.
    """
    records = {}
    for cond, roc_data in roc_dict.items():
        for spec in specificities:
            col = f"{cond}_sens{int(spec * 100)}"
            records[col] = {d['name']: sens_at_spec(d['fpr'], d['tpr'], spec) for d in roc_data}
    df = pd.DataFrame(records)
    df.index.name = 'feature'
    return df
