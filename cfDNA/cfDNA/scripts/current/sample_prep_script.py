# #!/usr/bin/env bash
# set -euo pipefail
import sys
import pandas as pd

# # input is sample fragment file and sample id
# SAMPLE_FILE="$1"     
# SAMPLE_ID="$2" 

if __name__ == '__main__':
    SAMPLE_FILE = sys.argv[1]
    SAMPLE_ID = sys.argv[2]

# # filtering autosomes and adding chr preifx (piped)
# zcat $SAMPLE_FILE | \
# awk -F"\t" 'BEGIN{OFS="\t"} $1 ~ /^[0-9]+$/ && $1 <= 22 {
#     chrom=$1
#     if (chrom !~ /^chr/) chrom="chr"chrom
#     $1=chrom
#     print
# }' > ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed

df_sample = pd.read_csv(SAMPLE_FILE, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

# # create fragment centroid file from sample (moving centroids to 2 and 3 positions so intersect works for the centroids not fragment boundaries)
# awk 'BEGIN{OFS="\t"} {centroid=int(($2+$3)/2); print $1, centroid, centroid, $2, $3, $4, $5}' ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed > ./data/sample_temp/${SAMPLE_ID}_fragments_centroids.bed

# bedtools intersect -a ./data/sample_temp/${SAMPLE_ID}_fragments_centroids.bed -b ./data/processing/openchrom_with_id.bed -wa -wb > ./data/processing/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed

# this is the file that will be used for most of the analysis
# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# # intersect fragmet centroids with open chromatin regions with ids
# awk 'BEGIN{OFS="\t"} {centroid=int(($2+$3)/2); print $1, centroid, centroid, $2, $3, $4, $5}' \
#     ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed | \
# bedtools intersect -a - -b ./data/processing/openchrom_with_id.bed -wa -wb \
#     > ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed

# # for fragment ends, we need to create a separate bedfile, intersecting openchrom regions with fragment start and ends, can t use the centroid file
# # create 2 lines from each fragment - one start, one end
# awk 'BEGIN{OFS="\t"} {u=int(($2)); print $1,u,u+1,"U"; d=int($3)-1; print $1,d,d+1,"D"}' \
#     ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed | \
# bedtools intersect -a - -b ./data/processing/openchrom_with_id.bed -wa -wb \
#     > ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed


# awk 'BEGIN{OFS="\t"} {
#     region_id = $11
#     print $1, $4, $4+1, "U", region_id
#     print $1, $5, $5+1, "D", region_id
# }' ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
# > ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed

# # Match with centroids
# awk 'BEGIN{OFS="\t"}
# NR==FNR {
#     centroid[$4] = $4
#     next
# }
# {
#     region_id = $5
#     if (region_id in centroid) {
#         chrom = $1
#         end1 = $2
#         end2 = $3
#         end_type = $4
#         rel_pos = end1 - centroid[region_id]
#         print chrom, end1, end2, end_type, ".", ".", region_id, centroid[region_id], rel_pos
#     }
# }' ./data/processing/openchrom_with_id_centroid.bed \
#    ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed \
# > ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed

# python ./scripts/current/sample_features_script.py \
#   $SAMPLE_ID \
#   ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
#   ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed \
#   ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed

# if [ $? -eq 0 ]; then
#     echo "job finished succesfully, deleting temp files"
#     rm ./data/sample_temp/${SAMPLE_ID}_autosomes_chr.bed
#     rm ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed
#     rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed
#     rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed
# else
#     echo "job failed"
#     exit 1
# fi