#!/bin/bash
#SBATCH --job-name=cristiano_dataset
#SBATCH --output=logs/dl_%A_%a.out
#SBATCH --error=logs/dl_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --array=1-459%20
# array length needs to be adjusted to the number of entries in the manifest file

# slurm array job - calls sample_download_script.sh which downlaods each sample defined in manifest 
set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

bash scripts/sample_download_script.sh "$SLURM_ARRAY_TASK_ID" 

