#!/bin/bash
#SBATCH --job-name=recompute_ocf
#SBATCH --output=logs/recompute_ocf_%A_%a.out
#SBATCH --error=logs/recompute_ocf_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# Get sample from manifest
MANIFEST="./data/manifest/Cristiano_manifest.csv"
SEQRUN_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | cut -d',' -f1)
SAMPLE_ID="EE${SEQRUN_ID}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

if [ ! -f "$FRAG_FILE" ]; then
    echo "Fragment file missing for $SAMPLE_ID"
    exit 1
fi

bash scripts/current/recompute_ocf.sh "$FRAG_FILE" "$SAMPLE_ID"

echo "Complete"