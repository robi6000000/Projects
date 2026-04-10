#!/bin/bash
#SBATCH --job-name=sample_prep
#SBATCH --output=logs/sample_%A_%a.out
#SBATCH --error=logs/sample_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

FRAG_FILE=$1
SAMPLE_ID=$2
MAPQ_FILTER=${3:-none}
FEATURES_PATH=${4:-./data/internal_features}
export FEATURES_PATH

echo "processing $SAMPLE_ID with mapq_filter $MAPQ_FILTER features_path $FEATURES_PATH"

bash scripts/generalized_processing/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER" true