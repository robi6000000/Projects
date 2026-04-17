#!/bin/bash
#SBATCH --job-name=mapq_fill
#SBATCH --output=logs/mapq%j_%A_%a.out
#SBATCH --error=logs/mapq%j_%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MAPQ_FILTER="${1:-none}"
RERUN="${2:-true}"

# TODO - internal data has form gsca_pca_00125_01_01_pl.hg19.frag.bed.gz
# need to adjust the formatting to be generalized
MANIFEST="./data/manifest/Cristiano_manifest.csv"
SEQRUN_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | cut -d',' -f1)
SAMPLE_ID="EE${SEQRUN_ID}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

echo "processing $SAMPLE_ID with mapq_filter $MAPQ_FILTER rerun $RERUN"

bash scripts/current/sample_prep_script.sh "$FRAG_FILE" "$SAMPLE_ID" "$MAPQ_FILTER" "$RERUN"

echo "finished $SAMPLE_ID"