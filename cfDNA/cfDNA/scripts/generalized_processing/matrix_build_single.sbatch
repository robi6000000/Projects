#!/bin/bash
#SBATCH --job-name=edm_build
#SBATCH --output=logs/edm_build_%j.out
#SBATCH --error=logs/edm_build_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

FEATURE=$1
MAPQ_FILTER=${2:-None}
FEATURES_FOLDER=${3:-./data/internal_features}
MATRIX_FOLDER=${4:-./data/internal_matrix}
METADATA_PATH=${5:-./data/manifest/internal_metadata.csv}
export FEATURES_FOLDER
echo "Building $FEATURE matrix with MAPQ filter $MAPQ_FILTER using features from $FEATURES_FOLDER"
export MATRIX_FOLDER
echo "Building matrix for $FEATURE with MAPQ filter $MAPQ_FILTER at $MATRIX_FOLDER"
echo "Using metadata from $METADATA_PATH"
python scripts/generalized_processing/matrix_build.py "$FEATURE" "$MAPQ_FILTER" "$FEATURES_FOLDER" "$MATRIX_FOLDER" "$METADATA_PATH"

echo "$FEATURE complete"