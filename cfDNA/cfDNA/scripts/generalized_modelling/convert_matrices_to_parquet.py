"""
One-time conversion of all feature matrix CSVs to Parquet (float32).
Originals are not deleted. Existing parquet files are skipped.

Usage:
    python scripts/generalized_modelling/convert_matrices_to_parquet.py [matrix_dir]

Default matrix_dir: data/matrix/by_feature
"""

import os
import sys
import glob
import pandas as pd
import numpy as np

matrix_dir = sys.argv[1] if len(sys.argv) > 1 else "data/matrix/by_feature"
metadata_cols = ['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true',
                 'seqrun_id', 'sample_name', 'tissue']

csv_files = sorted(glob.glob(os.path.join(matrix_dir, "*.csv")))
print(f"Found {len(csv_files)} CSV files in {matrix_dir}")

for csv_path in csv_files:
    parquet_path = csv_path.replace(".csv", ".parquet")
    if os.path.exists(parquet_path):
        print(f"  SKIP (exists): {os.path.basename(parquet_path)}")
        continue

    print(f"  Converting: {os.path.basename(csv_path)} ... ", end="", flush=True)
    df = pd.read_csv(csv_path, index_col=0, low_memory=False)

    # Cast numeric columns to float32, leave metadata as-is
    numeric_cols = [c for c in df.columns if c not in metadata_cols]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').astype(np.float32)

    df.to_parquet(parquet_path)
    size_mb = os.path.getsize(parquet_path) / (1024 ** 2)
    print(f"done ({size_mb:.0f} MB)")

print("Conversion complete.")
