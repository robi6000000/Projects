#!/bin/bash
#SBATCH --job-name=downlaod_corrupted
#SBATCH --output=logs/redownload_%A_%a.out
#SBATCH --error=logs/redownload_%A_%a.err
#SBATCH --array=1-160%20          # Adjust first number based on how many corrupted
#SBATCH --time=00:30:00           # 30 min per file
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# Get the manifest line number for this array task
INDEX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" corrupted_indices.txt)

echo "Array task ${SLURM_ARRAY_TASK_ID}: Re-downloading manifest line ${INDEX}"

# Delete existing corrupted files first
MANIFEST=./data/manifest/Cristiano_manifest.csv
LINE=$(tail -n +2 "$MANIFEST" | sed -n "${INDEX}p")
IFS=',' read -r seqrun_id rest <<< "$LINE"

SAMPLE_ID="EE${seqrun_id}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

# Delete corrupted files
if [ -f "$FRAG_FILE" ]; then
    echo "Deleting corrupted frag file: $FRAG_FILE"
    rm "$FRAG_FILE"
fi


# Re-download
bash ./scripts/sample_download_script.sh "$INDEX"

echo "Re-downloaded: $SAMPLE_ID"