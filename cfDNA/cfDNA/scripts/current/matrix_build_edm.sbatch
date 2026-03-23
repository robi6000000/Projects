#!/bin/bash
#SBATCH --job-name=build_edm
#SBATCH --output=logs/build_edm_%j.out
#SBATCH --error=logs/build_edm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna


python scripts/current/matrix_build.py edm

echo "edm complete"