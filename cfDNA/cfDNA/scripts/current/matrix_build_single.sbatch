#!/bin/bash
#SBATCH --job-name=build_single
#SBATCH --output=logs/build_single_%j.out
#SBATCH --error=logs/build_single_%j.err
#SBATCH --time=05:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

FEATURE=$1
MAPQ_FILTER=${2:-None}
echo "Building matrix for $FEATURE with MAPQ filter $MAPQ_FILTER"
python scripts/current/matrix_build.py "$FEATURE" "$MAPQ_FILTER"

echo "$FEATURE complete"