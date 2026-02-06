#!/bin/bash
#SBATCH --job-name=cristiano_dataset
#SBATCH --output=logs/dl_%A_%a.out
#SBATCH --error=logs/dl_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --array=1-540%20

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

bash scripts/sample_download_script.sh "$SLURM_ARRAY_TASK_ID" 
