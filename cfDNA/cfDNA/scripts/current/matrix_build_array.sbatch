#!/bin/bash
#SBATCH --job-name=build_matrices
#SBATCH --output=logs/build_matrix_%A_%a.out
#SBATCH --error=logs/build_matrix_%A_%a.err
#SBATCH --array=0-11%4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
MAPQ_FILTER=${1:-None}
echo "Using MAPQ filter: $MAPQ_FILTER"

# Map array task ID to feature name
FEATURES=(length edm iedm eedm fsd fsr pfe coverage ends ocf ifs wps)

FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}

echo "$SLURM_ARRAY_TASK_ID: building $FEATURE matrix"

python scripts/current/matrix_build.py "$FEATURE" "$MAPQ_FILTER"

echo "$FEATURE matrix complete"