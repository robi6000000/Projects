#!/bin/bash
#SBATCH --job-name=features_mapq
#SBATCH --output=logs/features_mapq%j_%A_%a.out
#SBATCH --error=logs/features_mapq%j_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MAPQ_FILTER="${MAPQ_FILTER:-none}"

MANIFEST="./data/manifest/Cristiano_manifest.csv"
SEQRUN_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | cut -d',' -f1)
SAMPLE_ID="EE${SEQRUN_ID}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

echo "processing $SAMPLE_ID with mapq_filter $MAPQ_FILTER"

bash scripts/current/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER"

echo "finished $SAMPLE_ID"