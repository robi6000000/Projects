#!/bin/bash
#SBATCH --job-name=build_matrices
#SBATCH --output=logs/build_matrix_%A_%a.out
#SBATCH --error=logs/build_matrix_%A_%a.err
#SBATCH --array=0-13
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
MAPQ_FILTER=$1
FEATURES_FOLDER=$2
MATRIX_FOLDER=$3
METADATA_PATH=$4
echo "Using MAPQ filter: $MAPQ_FILTER"
echo "Using features from $FEATURES_FOLDER"
echo "Using metadata from $METADATA_PATH"

FEATURES=(length edm iedm eedm eoedm cposedm fsd pfe coverage ends ocf ifs wps fsr)

FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}

echo "$SLURM_ARRAY_TASK_ID: building $FEATURE matrix"

python scripts/generalized_processing/matrix_build.py "$FEATURE" "$MAPQ_FILTER" "$FEATURES_FOLDER" "$MATRIX_FOLDER" "$METADATA_PATH"

echo "$FEATURE matrix complete"