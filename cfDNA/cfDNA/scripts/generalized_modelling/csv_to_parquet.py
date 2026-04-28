"""
Convert a single feature matrix CSV to Parquet (float32).
Usage: python csv_to_parquet.py <csv_path>
"""
import sys
import pandas as pd
import numpy as np
import os
import time

if len(sys.argv) < 2:
    print("Usage: python csv_to_parquet.py <csv_path>")
    sys.exit(1)

csv_path = sys.argv[1]
parquet_path = csv_path.replace(".csv", ".parquet")
metadata_cols = ['sample_id', 'disease', 'dataset', 'material', 'stage', 'cancer_true',
                 'seqrun_id', 'sample_name', 'tissue']

if os.path.exists(parquet_path):
    print(f"Parquet already exists, skipping: {parquet_path}")
    sys.exit(0)

csv_mb = os.path.getsize(csv_path) / 1024**2
print(f"Reading {csv_path}  ({csv_mb:.0f} MB) ...")
t0 = time.perf_counter()
df = pd.read_csv(csv_path, index_col=0, low_memory=False)
numeric_cols = [c for c in df.columns if c not in metadata_cols]
df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').astype(np.float32)
read_time = time.perf_counter() - t0
print(f"Read in {read_time:.1f}s  shape={df.shape}")

print(f"Writing {parquet_path} ...")
t1 = time.perf_counter()
df.to_parquet(parquet_path)
write_time = time.perf_counter() - t1
parquet_mb = os.path.getsize(parquet_path) / 1024**2
print(f"Done in {write_time:.1f}s  →  {parquet_mb:.0f} MB  ({parquet_mb/csv_mb*100:.1f}% of CSV)")
