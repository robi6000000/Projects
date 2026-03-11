#!/bin/bash
#SBATCH --job-name=build_fsr
#SBATCH --output=logs/build_fsr_%j.out
#SBATCH --error=logs/build_fsr_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

# script made solely because fsr matrix computatioin needs more memory and cant be run with matrix_build_array.sbatch
set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

python scripts/current/matrix_build.py fsr

echo "Fsr complete"