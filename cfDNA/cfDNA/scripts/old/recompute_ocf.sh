#!/usr/bin/env bash
set -euo pipefail

# TODO unused/ old file
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

# Create fragment ends with region_id from centroid intersection
# Create fragment ends from centroid intersection
awk 'BEGIN{OFS="\t"} {
    region_id = $11
    print $1, $4, $4+1, "U", region_id
    print $1, $5, $5+1, "D", region_id
}' ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
> ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed

# Load centroid file into awk array and match
awk 'BEGIN{OFS="\t"}
NR==FNR {
    centroid[$4] = $4
    next
}
{
    region_id = $5
    if (region_id in centroid) {
        chrom = $1
        end1 = $2
        end2 = $3
        end_type = $4
        rel_pos = end1 - centroid[region_id]
        print chrom, end1, end2, end_type, ".", ".", region_id, centroid[region_id], rel_pos
    }
}' ./data/processing/openchrom_with_id_centroid.bed \
   ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed \
> ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed

# Run Python with OCF recompute only
python ./scripts/current/sample_features_script.py \
  $SAMPLE_ID \
  ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed

# Cleanup
rm ./data/sample_temp/${SAMPLE_ID}_*.bed

echo "OCF recomputed for $SAMPLE_ID"