#!/bin/bash
#SBATCH --job-name=re_ocf
#SBATCH --output=logs/ocf_%A_%a.out
#SBATCH --error=logs/ocf_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

METADATA="${1:-./data/manifest/internal_metadata.csv}"

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")
SAMPLE_ID=$(echo "$ROW" | cut -d',' -f1)
FRAG_FILE=$(echo "$ROW" | cut -d',' -f5)

bash scripts/generalized_processing/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" none true