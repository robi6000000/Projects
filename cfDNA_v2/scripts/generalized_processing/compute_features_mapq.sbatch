#!/bin/bash
#SBATCH --job-name=features
#SBATCH --output=logs/features_%j_%A_%a.out
#SBATCH --error=logs/features_%j_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# usecase:
# sbatch --array=1-N%20 compute_features_mapq.sbatch 30 false data/metadata/Cristiano_metadata.csv data/cristiano_features

#   $1 MAPQ_FILTER   (e.g. 30, none)
#   $2 RERUN         (true|false — skip samples already processed if false)
#   $3 METADATA      (path to manifest CSV; column 5 must hold the fragment file path)
#   $4 FEATURES_PATH (output dir for per-sample feature CSVs)

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA_v2

MAPQ_FILTER="${1:-none}"
RERUN="${2:-true}"
METADATA="${3:-./data/metadata/internal_metadata.csv}"
FEATURES_PATH="${4:-./data/internal_features}"
export FEATURES_PATH

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")
SAMPLE_ID=$(echo "$ROW" | cut -d',' -f1)
FRAG_FILE=$(echo "$ROW" | cut -d',' -f5)

echo "processing $SAMPLE_ID with mapq_filter $MAPQ_FILTER rerun $RERUN features_path $FEATURES_PATH"

bash scripts/generalized_processing/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER" "$RERUN"

echo "finished $SAMPLE_ID"