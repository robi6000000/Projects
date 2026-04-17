#!/bin/bash
#SBATCH --job-name=rerun_missing
#SBATCH --output=logs/rerun_%A_%a.out
#SBATCH --error=logs/rerun_%A_%a.err
#SBATCH --array=1-1  
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

MISSING_SAMPLES=(
	"EE88325"
)
SAMPLE_ID="${MISSING_SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
MAPQ_FILTER=${1:-30}
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

echo "Processing sample: $SAMPLE_ID (task ${SLURM_ARRAY_TASK_ID})"

bash scripts/current/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER"

echo "Completed: $SAMPLE_ID"