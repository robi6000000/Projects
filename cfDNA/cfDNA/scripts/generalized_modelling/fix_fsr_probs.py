import pandas as pd

manifest_path = "data/manifest/Cristiano_metadata.csv"
fsr_path = "data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/cv/fsr_cv_probs.csv"

fsr = pd.read_csv(fsr_path)
fsr[['sample_id', 'probability', 'label']].to_csv(fsr_path, index=False)
print(f"Reverted to {len(fsr)} rows with columns: sample_id, probability, label")
print(fsr[['sample_id', 'probability', 'label']].head(3).to_string(index=False))
