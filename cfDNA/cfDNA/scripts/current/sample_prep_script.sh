#!/usr/bin/env bash
set -euo pipefail

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

# input is sample fragment file and sample id
SAMPLE_FILE="$1"     
SAMPLE_ID="$2" 
MAPQ_FILTER="${3:-none}"

# filtering autosomes and adding chr preifx (piped)
zcat "$SAMPLE_FILE" | \
awk -F"\t" 'BEGIN{OFS="\t"} $1 ~ /^[0-9]+$/ && $1 <= 22 {
    chrom=$1
    if (chrom !~ /^chr/) chrom="chr"chrom
    $1 = chrom
    print $0, NR-1
}' > ./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed
# chrom, start, end, score, strand, frag_id

# bedtools intersect -a ./data/sample_temp/${SAMPLE_ID}_fragments_centroids.bed -b ./data/processing/openchrom_with_id.bed -wa -wb > ./data/processing/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed

# this is the file that will be used for most of the analysis
# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# intersect fragmet centroids with open chromatin regions with ids
awk 'BEGIN{OFS="\t"} {centroid=int(($2+$3)/2); print $1, centroid, centroid, $2, $3, $4, $5, $6}' \
    ./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed | \
bedtools intersect -a - -b ./data/processing/openchrom_with_id.bed -wa -wb \
    > ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed
# f_chrom, centroid1, centroid2, f_start, f_end, score, strand, frag_id, oc_chrom, oc_start, oc_end, region_id

# for fragment ends, we need to create a separate bedfile, intersecting openchrom regions with fragment start and ends, can t use the centroid file
# create 2 lines from each fragment - one start, one end
# U for upstream end, D for downastream end. (if +, othewise if -, start is D and end is U)

awk 'BEGIN{OFS="\t"} {
    chrom = $1
    start = $2
    end = $3
    score = $4
    strand = $5
    frag_id = $6

    print chrom, start, start+1, "U", score, strand, frag_id
    print chrom, end-1, end, "D", score, strand, frag_id

}' ./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed > ./data/sample_temp/${SAMPLE_ID}_frag_ends_U_D.bed
# f_chrom, end1, end2, end_type, score, strand, frag_id

bedtools intersect -a ./data/sample_temp/${SAMPLE_ID}_frag_ends_U_D.bed -b ./data/processing/openchrom_with_id.bed -wa -wb \
    > ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed
# f_chrom, end1, end2, end_type, score, strand, frag_id, oc_chrom, oc_start, oc_end, region_id

# WPS - 'For WPS,[44] according to the genomic coordinate position of each
# cfDNA fragment, a window of 120 bp was slid at 1 bp intervals, and the
# likelihood of each base pair being covered at the whole genome level, fully
# covered (+1), and partially covered (−1), was counted. The mean value of
# all loci within each open chromatin region was calculated. 

# we intersect the extended open chromatin regions with the fragments of the sample
# (wps_openchrom_with_id.bed is created once by static_script.sh)
# resulting columns: f_chrom, f_start, f_end, f_score, f_strand, f_id, oc_chrom, oc_start, oc_end, oc_region_id
bedtools intersect \
    -a ./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed \
    -b ./data/processing/wps_openchrom_with_id.bed -wa -wb \
    > ./data/sample_temp/${SAMPLE_ID}_frag_wps_intersect.bed
# f_chrom, f_start, f_end, score, strand, frag_id, oc_chrom, oc_start, oc_end, region_id

# OCF - we need 
awk 'BEGIN{OFS="\t"} {
    f_chrom = $1
    f_start = $4
    f_end = $5
    score = $6
    strand = $7
    frag_id = $8
    region_id = $12

    print f_chrom, f_start, f_start+1, "U", score, strand, frag_id, region_id
    print f_chrom, f_end-1, f_end, "D", score, strand, frag_id, region_id
}' ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
> ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed
# chrom, end1, end2, end_type, score, strand, frag_id, region_id

# Match with open-chromatin centroids
# openchrom_with_id_centroid.bed columns: 1=oc_chrom, 2=oc_start, 3=oc_end, 4=centroid, 5=region_id
awk 'BEGIN{OFS="\t"}
NR==FNR {
    oc_start[$5] = $2
    oc_end[$5] = $3
    centroid[$5] = $4
    next
}
{
    chrom = $1
    end1 = $2
    end2 = $3
    end_type = $4
    score = $5
    strand = $6
    frag_id = $7
    region_id = $8

    if (region_id in centroid) {
        rel_pos = end1 - centroid[region_id]
        print chrom, end1, end2, end_type, score, strand, frag_id, oc_start[region_id], oc_end[region_id], region_id, centroid[region_id], rel_pos
    }
}' ./data/processing/openchrom_with_id_centroid.bed \
   ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed \
> ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed
# chrom, end1, end2, end_type, score, strand, frag_id, oc_start, oc_end, region_id, centroid, rel_pos

python ./scripts/current/sample_features_script.py \
  $SAMPLE_ID \
  ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed \
  ./data/sample_temp/${SAMPLE_ID}_frag_wps_intersect.bed \
  $MAPQ_FILTER

echo "job finished successfully"
# echo "deleting temp files"
# rm ./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_wps_intersect.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_U_D.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region_sorted.bed

