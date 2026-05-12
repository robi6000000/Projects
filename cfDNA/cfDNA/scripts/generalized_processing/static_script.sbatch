#!/bin/bash
#SBATCH --job-name=static_script
#SBATCH --output=logs/static_script_%j.out
#SBATCH --error=logs/static_script_%j.err
#SBATCH --mem=16G

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/scripts/generalized_processing

bash static_script.sh
