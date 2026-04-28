#!/bin/bash
#
# GS cohort example (38 samples):
#   sbatch --array=1-38%10 scripts/generalized_processing/fill_missing_mapq30.sbatch \
#       30 ./data/manifest/internal_metadata_gs.csv ./data/internal_features_gs false

#SBATCH --job-name=fill_missing
#SBATCH --output=logs/fill_missing_%A_%a.out
#SBATCH --error=logs/fill_missing_%A_%a.err
#SBATCH --array=1-1
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONUNBUFFERED=1

MAPQ_FILTER="${1:-30}"
METADATA="${2:-./data/manifest/internal_metadata_gs.csv}"
FEATURES_PATH="${3:-./data/internal_features_gs}"
RERUN="${4:-false}"

export FEATURES_PATH

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")
SAMPLE_ID=$(echo "$ROW" | cut -d',' -f1)
FRAG_FILE=$(echo "$ROW" | cut -d',' -f5)

echo "Processing $SAMPLE_ID (task $SLURM_ARRAY_TASK_ID)"
echo "  mapq_filter=$MAPQ_FILTER  rerun=$RERUN"
echo "  features_path=$FEATURES_PATH"

bash scripts/generalized_processing/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER" "$RERUN"

echo "Finished $SAMPLE_ID"
