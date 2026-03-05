#!/bin/bash
# Find which samples failed and need re-downloading

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MANIFEST=./data/manifest/Cristiano_manifest.csv

# Get successful sample IDs
ls ./data/rows_sample_temp/*.csv | sed 's|.*/||' | sed 's|_features.csv||' | sed 's|EE||' | sort > successful_ids.txt

# Get all sample IDs from manifest
tail -n +2 "$MANIFEST" | cut -d',' -f1 | sort > all_ids.txt

# Find missing ones (need re-download)
comm -23 all_ids.txt successful_ids.txt > missing_ids.txt

# Convert to manifest line numbers (for array job)
> corrupted_indices.txt
while read seqrun_id; do
    # Find line number in manifest (add 1 because we skip header)
    LINE_NUM=$(tail -n +2 "$MANIFEST" | grep -n "^${seqrun_id}," | cut -d':' -f1)
    echo "$LINE_NUM" >> corrupted_indices.txt
done < missing_ids.txt

echo "Found $(wc -l < corrupted_indices.txt) corrupted samples"
echo "Indices saved to: corrupted_indices.txt"

# Cleanup
rm successful_ids.txt all_ids.txt missing_ids.txt