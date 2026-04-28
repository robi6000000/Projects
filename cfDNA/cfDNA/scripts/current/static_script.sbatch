#!/bin/bash
#SBATCH --job-name=static_script
#SBATCH --output=logs/static_script_%j.out
#SBATCH --error=logs/static_script_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/scripts/current

bash static_script.sh
