#!/usr/bin/env bash
set -euo pipefail

SAMPLE_FILE="$1"
SAMPLE_ID="$2"

# Regenerate only the files needed for OCF
echo "Regenerating OCF data for $SAMPLE_ID..."

# Create autosomes file
zcat $SAMPLE_FILE | \
awk -F"\t" 'BEGIN{OFS="\t"} $1 ~ /^[0-9]+$/ && $1 <= 22 {
    chrom=$1
    if (chrom !~ /^chr/) chrom="chr"chrom
    $1=chrom
    print
}' > ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed

# Create centroids intersect
awk 'BEGIN{OFS="\t"} {centroid=int(($2+$3)/2); print $1, centroid, centroid, $2, $3, $4, $5}' \
    ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed | \
bedtools intersect -a - -b ./data/processing/openchrom_with_id.bed -wa -wb \
    > ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed

# FIXED: Create fragment ends with region_id from centroid intersection
awk 'BEGIN{OFS="\t"} {
    print $1, $4, $4+1, "U", $11;
    print $1, $5, $5+1, "D", $11
}' ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed | \
sort -k5,5n > ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed

# Join with centroid file on region_id
join -1 5 -2 4 -t $'\t' \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed \
  <(sort -k4,4n ./data/processing/openchrom_with_id_centroid.bed) | \
awk 'BEGIN{OFS="\t"} {
    region_id = $1
    chrom = $2
    end1 = $3
    end2 = $4
    end_type = $5
    centroid = $9
    rel_pos = end1 - centroid
    print chrom, end1, end2, end_type, $7, $8, region_id, centroid, rel_pos
}' > ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed

# Run Python with OCF recompute only
python ./scripts/current/sample_features_script.py \
  $SAMPLE_ID \
  ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed

# Cleanup
rm ./data/sample_temp/${SAMPLE_ID}_*.bed

echo "OCF recomputed for $SAMPLE_ID"