#!/bin/bash
#SBATCH --job-name=sample_prep
#SBATCH --output=logs/sample_%A_%a.out
#SBATCH --error=logs/sample_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

SAMPLE_ID=$1
MAPQ_FILTER=$2
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

bash scripts/current/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER"