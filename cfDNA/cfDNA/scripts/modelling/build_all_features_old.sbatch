#!/bin/bash
#SBATCH --job-name=build_features
#SBATCH --output=logs/build_features_%A_%a.out
#SBATCH --error=logs/build_features_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=01:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# Map array task ID to feature name
FEATURES=(length pfe fsr fsd coverage ends ocf ifs wps edm)
FEATURE=${FEATURES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Array task $SLURM_ARRAY_TASK_ID: Building $FEATURE matrix"

bash scripts/modelling/build_feature_matrix_old.sh "$FEATURE"

echo "$FEATURE complete"