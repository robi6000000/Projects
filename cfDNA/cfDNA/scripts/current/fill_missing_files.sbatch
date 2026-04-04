#!/bin/bash
#SBATCH --job-name=fill_features
#SBATCH --output=logs/fill_features_%A_%a.out
#SBATCH --error=logs/fill_features_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

MANIFEST="./data/manifest/Cristiano_manifest.csv"
SEQRUN_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | cut -d',' -f1)
SAMPLE_ID="EE${SEQRUN_ID}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

bash scripts/current/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID"