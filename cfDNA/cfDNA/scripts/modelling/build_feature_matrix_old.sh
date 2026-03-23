#!/bin/bash
set -euo pipefail

FEATURE=$1  #length, pfe, fsr, fsd, coverage, ends, ocf, ifs, wps, edm
FEATURE_DIR="./data/cristiano_features/${FEATURE}"
OUTPUT_DIR="./data/matrix/by_feature"
MANIFEST="./data/manifest/Cristiano_manifest.csv"

mkdir -p "$OUTPUT_DIR"

echo "Building HDF5 matrix for: $FEATURE"

# get feature files
FEATURE_FILES=$(ls "$FEATURE_DIR"/*.csv | sort)
FILE_COUNT=$(echo "$FEATURE_FILES" | wc -l)
echo "Found $FILE_COUNT files"

# get structure from first file
FIRST_FILE=$(echo "$FEATURE_FILES" | head -1)
SAMPLE_ID=$(basename "$FIRST_FILE" | sed "s/_${FEATURE}.csv//")

echo "Processing first file: $SAMPLE_ID"

# export 
export FEATURE
export FIRST_FILE
export SAMPLE_ID

# Read and flatten first sample
python3 << 'PYTHON_FLATTEN_EOF'
import pandas as pd
import numpy as np
import os

feature = os.environ['FEATURE']
first_file = os.environ['FIRST_FILE']
sample_id = os.environ['SAMPLE_ID']

# Load feature matrix
df = pd.read_csv(first_file, index_col=0)

# fix column names
column_names = []
for row_idx in df.index:
    for col_name in df.columns:
        # Clean up column name
        clean_col = str(col_name).replace('[', '').replace(']', '').replace('(', '').replace(')', '').replace(', ', '_').replace(' ', '_')
        column_names.append(f"{feature}_{row_idx}_bin_{clean_col}")

# Flatten to 1D vector
feature_vector = df.values.flatten()

# Create dataframe with meaningful names
result = pd.DataFrame([feature_vector], columns=column_names)
result.insert(0, 'sample_id', sample_id)

# Save to temp CSV (for incremental building)
result.to_csv(f"/tmp/{feature}_matrix_temp.csv", index=False)
print(f"Flattened {df.shape} to {len(feature_vector)} features")
PYTHON_FLATTEN_EOF

# Append all other samples
COUNT=1
export FEATURE
for file in $FEATURE_FILES; do
    if [ "$file" == "$FIRST_FILE" ]; then
        continue
    fi
    
    export SAMPLE_ID=$(basename "$file" | sed "s/_${FEATURE}.csv//")
    export file
    
    if [ $((COUNT % 50)) -eq 0 ]; then
        echo "  Processed $COUNT/$FILE_COUNT..."
    fi
    
    python3 << 'PYTHON_APPEND_EOF'
import pandas as pd
import numpy as np
import os

feature = os.environ['FEATURE']
file_path = os.environ['file']
sample_id = os.environ['SAMPLE_ID']

# Load and flatten
df = pd.read_csv(file_path, index_col=0)
feature_vector = df.values.flatten()

# Append to file
with open(f"/tmp/{feature}_matrix_temp.csv", 'a') as f:
    f.write(sample_id + "," + ",".join(map(str, feature_vector)) + "\n")
PYTHON_APPEND_EOF
    
    ((COUNT++))
done

echo "Processed $FILE_COUNT files"

# Convert to HDF5 and add metadata
echo "Converting to HDF5 and adding metadata..."

export OUTPUT_DIR
export MANIFEST

python3 << 'CONVERT_HDF5_EOF'
import pandas as pd
import os

feature = os.environ['FEATURE']
output_dir = os.environ['OUTPUT_DIR']
manifest_path = os.environ['MANIFEST']

print("Loading temp CSV and adding metadata...")

# Load temp CSV directly (skip intermediate HDF5)
temp_csv = f"/tmp/{feature}_matrix_temp.csv"
matrix = pd.read_csv(temp_csv)

# Load manifest and add metadata
manifest = pd.read_csv(manifest_path)

matrix['seqrun_id'] = matrix['sample_id'].str.replace('EE', '').astype(int)
matrix = matrix.merge(
    manifest[['seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']],
    on='seqrun_id',
    how='left'
)

# Reorder: metadata first
metadata_cols = ['sample_id', 'seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']
feature_cols = [c for c in matrix.columns if c not in metadata_cols]
matrix = matrix[metadata_cols + feature_cols]

# Convert string columns to object dtype (HDF5 compatible)
for col in ['sample_id', 'sample_name', 'stage', 'disease', 'tissue']:
    matrix[col] = matrix[col].astype('object')

# Save final HDF5 with 'fixed' format
output_h5 = f"{output_dir}/matrix_{feature}.h5"
matrix.to_hdf(output_h5, key='data', mode='w', complevel=5, format='fixed')

print(f"✓ Saved: {matrix.shape}")
print(f"✓ File: {output_h5}")

# Cleanup temp files
os.remove(temp_csv)
CONVERT_HDF5_EOF

echo "Complete: $OUTPUT_DIR/matrix_${FEATURE}.h5"