#!/bin/bash
#SBATCH --job-name=rerun_fsr
#SBATCH --output=logs/rerun_fsr%A_%a.out
#SBATCH --error=logs/rerun_fsr%A_%a.err
#SBATCH --array=1-459%20
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MANIFEST="./data/manifest/Cristiano_manifest.csv"
SEQRUN_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | cut -d',' -f1)
SAMPLE_ID="EE${SEQRUN_ID}"
FRAG_FILE="./data/source/cristiano/${SAMPLE_ID}.hg19.frag.tsv.bgz"

# Recreate only the centroid intersection file needed for FSD/FSR
zcat "$FRAG_FILE" | \
awk -F"\t" 'BEGIN{OFS="\t"} $1 ~ /^[0-9]+$/ && $1 <= 22 {
    chrom=$1; if (chrom !~ /^chr/) chrom="chr"chrom; $1=chrom; print
}' | \
awk 'BEGIN{OFS="\t"} {centroid=int(($2+$3)/2); print $1, centroid, centroid, $2, $3, $4, $5}' | \
bedtools intersect -a - -b ./data/processing/openchrom_with_id.bed -wa -wb \
    > ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed

# Rerun only fsr
python scripts/current/sample_features_script.py \
    "$SAMPLE_ID" \
    "./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed" \
    "./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed" \
    "./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed"

rm ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed
echo "Done: $SAMPLE_ID"