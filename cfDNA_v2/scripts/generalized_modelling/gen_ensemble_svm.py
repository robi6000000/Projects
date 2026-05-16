import os
import sys
import gc
import pickle
import pandas as pd
from modules.gen_cfdna_model import CFDNAModel

if __name__ == "__main__":
    if len(sys.argv) < 9:
        print("Incorrect arguments")
        sys.exit(1)

    option = sys.argv[1].lower() # cv, train, test
    meta_matrix_path = sys.argv[2] # path to metamatrix csv
    kernel = sys.argv[3] # kernel type for svm, linear/rbf
    cv_repeats = int(sys.argv[4]) # number of repeats for CV, only used for cv option
    dataset_tag = None if sys.argv[5].lower() == 'none' else sys.argv[5] # for specifying test datset
    model_path_arg = None if sys.argv[6].lower() == 'none' else sys.argv[6] # for specifying trained model path
    output_path_arg = None if sys.argv[7].lower() == 'none' else sys.argv[7] # for specifying test output path
    run_tag = None if sys.argv[8].lower() == 'none' else sys.argv[8] # optional suffix for directory/files
    # originally ensemble SVM had standardization on by default, currently comparing results without std. 
    # ensemble input matrix is already a matrix of probabilities
    standardize = sys.argv[9].lower() == 'true' if len(sys.argv) >= 10 else True


    parent_output_dir = os.path.dirname(os.path.dirname(meta_matrix_path))
    ensemble_dir_name = f"ensemble_{run_tag}" if run_tag else "ensemble"

    base_output_dir = os.path.join(parent_output_dir, ensemble_dir_name)
    cv_dir = os.path.join(base_output_dir, "cv")
    model_dir = os.path.join(base_output_dir, "models")
    test_dir = os.path.join(base_output_dir, "test")
    os.makedirs(cv_dir,    exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(test_dir,  exist_ok=True)

    model_path = model_path_arg  or os.path.join(model_dir, "ensemble.pkl")
    _test_stem = f"{dataset_tag}_preds" if dataset_tag else "preds"
    test_output_path = output_path_arg or os.path.join(test_dir, f"{_test_stem}.csv")

    print(f"Config: option-{option}, kernel={kernel}, cv_repeats={cv_repeats}, "
          f"dataset_tag={dataset_tag}, run_tag={run_tag}, standardize={standardize}")
    print(f"Output dir: {base_output_dir}")

    mx = pd.read_csv(meta_matrix_path, index_col=0, low_memory=False)
    print(f"Meta-matrix loaded: shape={mx.shape}")

    model = CFDNAModel(mx, gc_content=None, kernel=kernel,
                       gc_correction=False, pca=False, cv_repeats=cv_repeats,
                       standardize=standardize)
    del mx
    gc.collect()

    if option == "cv":
        print("Cross-validating ensemble SVM")
        results = model.cv_svm()
        output = pd.DataFrame({
            "probability": results["probabilities"],
            "label":       results["labels"]
        }, index=results["sample_ids"])
        out_path = os.path.join(cv_dir, "cv_probs.csv")
        output.to_csv(out_path)
        print(f"CV results saved: {out_path}")

    elif option == "train":
        print("Training ensemble SVM on full meta-matrix")
        trained_model = model.train_model()
        with open(model_path, "wb") as f:
            pickle.dump(trained_model, f)
        print(f"Trained model saved: {model_path}")

    elif option == "test":
        print("Predicting with ensemble SVM")
        with open(model_path, "rb") as f:
            trained_model = pickle.load(f)
        output = model.predict(trained_model)
        output.to_csv(test_output_path)
        print(f"Predictions saved: {test_output_path}")

    else:
        raise ValueError(f"Unknown option: {option}. Use cv, train, or test.")
