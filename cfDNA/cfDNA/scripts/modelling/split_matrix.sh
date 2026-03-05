#!/bin/bash

echo "Building 10 separate matrices by feature type..."

ROWS_FOLDER="./data/rows_sample_temp"
MATRIX_FOLDER="./data/matrix"
MANIFEST="./data/manifest/Cristiano_manifest.csv"

mkdir -p "$MATRIX_FOLDER/by_feature"

# Get first file to identify column indices for each feature type
FIRST_FILE=$(ls "$ROWS_FOLDER"/*.csv | head -1)
HEADER=$(head -1 "$FIRST_FILE")

# Save header for reference
echo "$HEADER" > /tmp/full_header.txt

echo "Extracting feature columns by prefix..."
# Extract column names using bash (strip quotes)
echo "Extracting feature columns by prefix..."
FIRST_FILE=$(ls "$ROWS_FOLDER"/*.csv | head -1)
HEADER=$(head -1 "$FIRST_FILE" | tr ',' '\n' | sed 's/"//g')

# Count columns for each feature type
echo "$HEADER" | grep "^pfe_region_" | nl -nln > /tmp/pfe_cols.txt
echo "$HEADER" | grep "^cov_region_" | nl -nln > /tmp/coverage_cols.txt
echo "$HEADER" | grep "^end_region_" | nl -nln > /tmp/ends_cols.txt
echo "$HEADER" | grep "^ocf_region_" | nl -nln > /tmp/ocf_cols.txt
echo "$HEADER" | grep "^ifs_region_" | nl -nln > /tmp/ifs_cols.txt
echo "$HEADER" | grep "^wps_region_" | nl -nln > /tmp/wps_cols.txt

echo "$HEADER" | grep "^len_chr.*_bin_" | nl -nln > /tmp/length_cols.txt
echo "$HEADER" | grep "^fsd_chr.*_bin_" | nl -nln > /tmp/fsd_cols.txt
echo "$HEADER" | grep "^edm_chr.*_bin_" | nl -nln > /tmp/edm_cols.txt
echo "$HEADER" | grep "^fsr_.*_bin_" | grep -v "_chr" | nl -nln > /tmp/fsr_cols.txt

# Add sample_id column (column 1) to each
for feature in pfe coverage ends ocf ifs wps length fsd edm fsr; do
    echo "1" > /tmp/${feature}_cols_with_id.txt
    cat /tmp/${feature}_cols.txt >> /tmp/${feature}_cols_with_id.txt
    # Convert to comma-separated for cut command
    paste -sd, /tmp/${feature}_cols_with_id.txt > /tmp/${feature}_cols_cut.txt
    
    COUNT=$(wc -l < /tmp/${feature}_cols.txt)
    echo "$feature: $COUNT features"
done
# Build each feature-specific matrix
for FEATURE in length pfe fsr fsd coverage ends ocf ifs wps edm; do
    echo ""
    echo "Building $FEATURE matrix..."
    
    COLS=$(cat /tmp/${FEATURE}_cols.txt)
    OUTPUT="${MATRIX_FOLDER}/by_feature/matrix_${FEATURE}.csv"
    
    # Get header
    head -1 "$FIRST_FILE" | cut -d',' -f"$COLS" > "$OUTPUT"
    
    # Append data from all samples
    COUNT=0
    for file in "$ROWS_FOLDER"/*.csv; do
        tail -n +2 "$file" | cut -d',' -f"$COLS" >> "$OUTPUT"
        ((COUNT++))
        if [ $((COUNT % 50)) -eq 0 ]; then
            echo "  Processed $COUNT files..."
        fi
    done
    
    echo "  ✓ ${FEATURE}: $(wc -l < $OUTPUT) rows"
done

# Add metadata to each matrix using Python (fast - small matrices!)
echo ""
echo "Adding metadata to each matrix..."

python3 << 'ADD_METADATA'
import pandas as pd
import os

manifest = pd.read_csv('./data/manifest/Cristiano_manifest.csv')
matrix_folder = './data/matrix/by_feature'

for feature in ['length', 'pfe', 'fsr', 'fsd', 'coverage', 'ends', 'ocf', 'ifs', 'wps', 'edm']:
    print(f"Processing {feature}...")
    
    # Load feature matrix (fast - only one feature type!)
    df = pd.read_csv(f'{matrix_folder}/matrix_{feature}.csv')
    
    # Add metadata
    df['seqrun_id'] = df['sample_id'].str.replace('EE', '').astype(int)
    df = df.merge(
        manifest[['seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']],
        on='seqrun_id',
        how='left'
    )
    
    # Reorder: metadata first
    metadata_cols = ['sample_id', 'seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']
    feature_cols = [c for c in df.columns if c not in metadata_cols]
    df = df[metadata_cols + feature_cols]
    
    # Save
    df.to_csv(f'{matrix_folder}/matrix_{feature}_with_metadata.csv', index=False)
    print(f"  ✓ Saved: {df.shape}")
ADD_METADATA

# Cleanup
rm /tmp/*_cols.txt /tmp/full_header.txt

echo ""
echo "=========================================="
echo "✓ COMPLETE!"
echo "=========================================="
echo "Created 10 matrices in: $MATRIX_FOLDER/by_feature/"
ls -lh "$MATRIX_FOLDER/by_feature/" | grep "with_metadata"