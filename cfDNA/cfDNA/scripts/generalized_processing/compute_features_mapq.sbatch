#!/bin/bash
#SBATCH --job-name=features
#SBATCH --output=logs/features_%j_%A_%a.out
#SBATCH --error=logs/features_%j_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MAPQ_FILTER="${1:-none}"
RERUN="${2:-true}"
METADATA="${3:-./data/manifest/internal_metadata.csv}"
FEATURES_PATH="${4:-./data/internal_features}"
export FEATURES_PATH

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")
SAMPLE_ID=$(echo "$ROW" | cut -d',' -f1)
FRAG_FILE=$(echo "$ROW" | cut -d',' -f5)

echo "processing $SAMPLE_ID with mapq_filter $MAPQ_FILTER rerun $RERUN features_path $FEATURES_PATH"

bash scripts/generalized_processing/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER" "$RERUN"

echo "finished $SAMPLE_ID"