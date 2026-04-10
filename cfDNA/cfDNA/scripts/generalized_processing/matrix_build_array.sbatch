#!/bin/bash
#SBATCH --job-name=build_matrices
#SBATCH --output=logs/build_matrix_%A_%a.out
#SBATCH --error=logs/build_matrix_%A_%a.err
#SBATCH --array=0-13%4
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
MAPQ_FILTER=${1:-None}
FEATURES_FOLDER=${2:-./data/internal_features}
MATRIX_FOLDER=${3:-./data/internal_matrix}
METADATA_PATH=${4:-./data/manifest/internal_metadata.csv}
echo "Using MAPQ filter: $MAPQ_FILTER"
echo "Using features from $FEATURES_FOLDER"
echo "Using metadata from $METADATA_PATH"

# Map array task ID to feature name
FEATURES=(length edm poem_prem ext_poem_prem fsd fsr pfe coverage ends ocf ifs wps wps_compute )

FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}

echo "$SLURM_ARRAY_TASK_ID: building $FEATURE matrix"

python scripts/generalized_processing/matrix_build.py "$FEATURE" "$MAPQ_FILTER" "$FEATURES_FOLDER" "$MATRIX_FOLDER" "$METADATA_PATH"

echo "$FEATURE matrix complete"