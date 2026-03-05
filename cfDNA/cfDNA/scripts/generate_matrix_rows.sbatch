#!/bin/bash
#SBATCH --job-name=cristiano_matrix
#SBATCH --output=logs/matrix_%j.out
#SBATCH --error=logs/matrix_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# conda activate cfdna
python scripts/generate_matrix.py