
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
import sys

# start by renaming the file to evaluate.py

# the zhou scores will later be removed, only leaving out own auc scores after we get all the results(the jobs finish), and after we compare them. we dont want the final 
# results which we put into the thesis to have a table with zhou results plastered everywhere.
zhou_scores = {
    'length': 0.8741, 'edm': 0.9736, 'fsd': 0.9271,
    'coverage': 0.9638, 'ends': 0.9639, 'ocf': 0.9467,
    'ifs': 0.9653, 'wps': 0.9658, 'pfe': 0.9579, 'fsr': 0.9441
}

# the logic of this script should be that it will evaluate auc for all files in a given directory, which have _cv.csv at the end
# files in a given directory will be named like this 'feature_cv.csv' (if we are evaluating per-feature cv results)
# or like this 'ensemble_cv.csv' (if we are evaluating ensemble cv results) 
# an argument determining which of these it is will be passed: ensemble = true/false

def parse_filename(file, ensemble=False):
    filename = os.path.basename(file)
    if ensemble:
        feature = "ensemble"
    else:
        feature = filename.split("_")[0]
    return feature

def evaluate_cv_feature(file):
    df = pd.read_csv(file)
    # print(df.columns)
    y_true = df['label'].values
    y_pred = df['probability'].values
    auc = roc_auc_score(y_true, y_pred)
    return auc

def evaluate_cv_ensemble(file):
    # ensemble will only return a single probability per sample, meaning the result will just be 1 auc score, 
    # need to figure out what to do with this
    return 

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("incorrect number of arguments")
        sys.exit(1)
    ensemble = sys.argv[1].lower() == "true"
    cv_directory = sys.argv[2]
    out_filepath = sys.argv[3]
    results = []
    for file in os.listdir(cv_directory):
        if ensemble and file.endswith("_cv.csv"):
            feature = parse_filename(file, ensemble=True)
            auc = evaluate_cv_ensemble(os.path.join(cv_directory, file))
            results.append((feature, auc))
        elif not ensemble and file.endswith("_cv.csv"):
            feature = parse_filename(file, ensemble=False)
            auc = evaluate_cv_feature(os.path.join(cv_directory, file))
            results.append((feature, auc))
    results_df = pd.DataFrame(results, columns=["feature", "auc"])
    results_df.to_csv(out_filepath, index=False)


