#!/bin/bash
#SBATCH --job-name=build_matrices
#SBATCH --output=logs/build_matrix_%A_%a.out
#SBATCH --error=logs/build_matrix_%A_%a.err
#SBATCH --array=1-10%4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# Map array task ID to feature name
FEATURES=(length pfe fsr fsd coverage ends ocf ifs wps edm)
FEATURE=${FEATURES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "$SLURM_ARRAY_TASK_ID: building $FEATURE matrix"

python scripts/current/matrix_build.py "$FEATURE"

echo "$FEATURE matrix complete"