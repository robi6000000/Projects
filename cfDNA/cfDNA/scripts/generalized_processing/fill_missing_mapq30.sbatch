#!/bin/bash
# Usage:
#   python3 scripts/current/audit_feature_completeness.py
#   ARRAY_IDS=$(paste -sd, log/missing_mapq30_array_indices.txt)
#   sbatch --array="${ARRAY_IDS}%10" scripts/current/fill_missing_mapq30.sbatch

#SBATCH --job-name=fill_mapq30
#SBATCH --output=logs/fill_mapq30_%A_%a.out
#SBATCH --error=logs/fill_mapq30_%A_%a.err
#SBATCH --array=1-1
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MANIFEST="./data/manifest/Cristiano_manifest.csv"
SEQRUN_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | cut -d',' -f1)
SAMPLE_ID="EE${SEQRUN_ID}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

echo "Backfilling mapq30 features for ${SAMPLE_ID}"
bash scripts/current/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" 30
echo "Finished ${SAMPLE_ID}"
