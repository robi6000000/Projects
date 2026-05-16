#!/usr/bin/env bash
set -euo pipefail

cd /data/projects/liquid_biopsy/Projects/cfDNA_v2
source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONUNBUFFERED=1

# input is sample fragment file and sample id
SAMPLE_FILE="$1"     
SAMPLE_ID="$2" 
MAPQ_FILTER="${3:-none}"
RERUN="${4:-true}"

AUTOSOMES_BED="./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed"
CENTROIDS_INTERSECT_BED="./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed"
FRAG_ENDS_BED="./data/sample_temp/${SAMPLE_ID}_frag_ends_U_D.bed"
FRAG_ENDS_INTERSECT_BED="./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed"
FRAG_WPS_INTERSECT_BED="./data/sample_temp/${SAMPLE_ID}_frag_wps_intersect.bed"
FRAG_ENDS_WITH_REGION_BED="./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed"
FRAG_ENDS_OCF_BED="./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed"

# if frag input starts with #, skip those lines (e.g. if there is a header)

# filtering autosomes and adding chr preifx (piped)
if [[ -s "$AUTOSOMES_BED" ]]; then
    echo "${SAMPLE_ID}_autosomes_chr_id.bed - reusing savedautosome-filtered fragment file"
else
    zcat "$SAMPLE_FILE" | \
    awk -F"\t" 'BEGIN{OFS="\t"} $1 !~ /^#/ && $1 ~ /^[0-9]+$/ && $1 <= 22 {
        chrom=$1
        # add chr prefix
        if (chrom !~ /^chr/) chrom="chr"chrom
        # Some samples from finaleDB had 6 columns (chrom,start,end,name,mapq,strand), but normally its (chrom,start,end,mapq,strand).
        if (NF >= 6) {
            score = $5
            strand = $6
        } else {
            score = $4
            strand = $5
        }
        print chrom, $2, $3, score, strand, NR-1
    }' > "$AUTOSOMES_BED"
    echo "${SAMPLE_ID}_autosomes_chr_id.bed - frag file filtered to autosomes and added indexing"
fi
# chrom, start, end, score, strand, frag_id

# bedtools intersect -a ./data/sample_temp/${SAMPLE_ID}_fragments_centroids.bed -b ./data/processing/openchrom_with_id.bed -wa -wb > ./data/processing/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed

# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# this is the file that will be used for most of the analysis:
# intersect fragmet centroids with open chromatin regions with ids
if [[ -s "$CENTROIDS_INTERSECT_BED" ]]; then
    echo "${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed - reusing existing centroid/open chromatin intersect"
else
    awk 'BEGIN{OFS="\t"} {centroid=int(($2+$3)/2); print $1, centroid, centroid, $2, $3, $4, $5, $6}' \
        "$AUTOSOMES_BED" | \
    bedtools intersect -a - -b ./data/processing/openchrom_with_id.bed -wa -wb \
        > "$CENTROIDS_INTERSECT_BED"
fi
# f_chrom, centroid1, centroid2, f_start, f_end, score, strand, frag_id, oc_chrom, oc_start, oc_end, region_id

# for fragment ends feature, we need to create a separate bedfile, 
# intersecting openchrom regions with fragment start and ends, cant use the centroid file
# we create 2 lines from each fragment - one for start coordinate, one for end coordinate
# U for upstream end, D for downastream end. 
if [[ -s "$FRAG_ENDS_BED" ]]; then
    echo "${SAMPLE_ID}_frag_ends_U_D.bed - reusing existing fragment ends file"
else
    awk 'BEGIN{OFS="\t"} {
        chrom = $1
        start = $2
        end = $3
        score = $4
        strand = $5
        frag_id = $6

        print chrom, start, start+1, "U", score, strand, frag_id
        print chrom, end-1, end, "D", score, strand, frag_id

    }' "$AUTOSOMES_BED" > "$FRAG_ENDS_BED"
    echo "${SAMPLE_ID}_frag_ends_U_D.bed - created fragment ends bed file with U/D annotation"
fi
# f_chrom, end1, end2, end_type, score, strand, frag_id

if [[ -s "$FRAG_ENDS_INTERSECT_BED" ]]; then
    echo "${SAMPLE_ID}_frag_ends_openchrom_intersect.bed - reusing existing fragment-end/open chromatin intersect"
else
    bedtools intersect -a "$FRAG_ENDS_BED" -b ./data/processing/openchrom_with_id.bed -wa -wb \
        > "$FRAG_ENDS_INTERSECT_BED"
    echo "${SAMPLE_ID}_frag_ends_openchrom_intersect.bed - intersected fragment ends with open chromatin regions"
fi
# f_chrom, end1, end2, end_type, score, strand, frag_id, oc_chrom, oc_start, oc_end, region_id

# WPS - 'For WPS,[44] according to the genomic coordinate position of each
# cfDNA fragment, a window of 120 bp was slid at 1 bp intervals, and the
# likelihood of each base pair being covered at the whole genome level, fully
# covered (+1), and partially covered (−1), was counted. The mean value of
# all loci within each open chromatin region was calculated. 

# we intersect the extended open chromatin regions with the fragments of the sample
# (wps_openchrom_with_id.bed is created once by static_script.sh)
# resulting columns: f_chrom, f_start, f_end, f_score, f_strand, f_id, oc_chrom, oc_start, oc_end, oc_region_id
if [[ -s "$FRAG_WPS_INTERSECT_BED" ]]; then
    echo "${SAMPLE_ID}_frag_wps_intersect.bed - reusing existing WPS/open chromatin intersect"
else
    bedtools intersect \
        -a "$AUTOSOMES_BED" \
        -b ./data/processing/wps_openchrom_with_id.bed -wa -wb \
        > "$FRAG_WPS_INTERSECT_BED"
    echo "${SAMPLE_ID}_frag_wps_intersect.bed - intersected fragments with wps padded open chromatin regions"
fi
# f_chrom, f_start, f_end, score, strand, frag_id, oc_chrom, oc_start, oc_end, region_id

# OCF 
if [[ -s "$FRAG_ENDS_WITH_REGION_BED" ]]; then
    echo "${SAMPLE_ID}_frag_ends_with_region.bed - reusing existing fragment ends with region annotation"
else
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
    }' "$CENTROIDS_INTERSECT_BED" \
    > "$FRAG_ENDS_WITH_REGION_BED"
    echo "${SAMPLE_ID}_frag_ends_with_region.bed - created fragment ends bed file with U/D annotation and region id from centroid intersect"
fi
# chrom, end1, end2, end_type, score, strand, frag_id, region_id

# Match with open-chromatin centroids
# openchrom_with_id_centroid.bed columns: 1=oc_chrom, 2=oc_start, 3=oc_end, 4=centroid, 5=region_id
if [[ -s "$FRAG_ENDS_OCF_BED" ]]; then
    echo "${SAMPLE_ID}_frag_ends_ocf.bed - reusing existing OCF-annotated fragment ends"
else
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
       "$FRAG_ENDS_WITH_REGION_BED" \
    > "$FRAG_ENDS_OCF_BED"
    echo "${SAMPLE_ID}_frag_ends_ocf.bed - fragment ends bed file with OCF annotation relative position column"
fi
# chrom, end1, end2, end_type, score, strand, frag_id, oc_start, oc_end, region_id, centroid, rel_pos

python -u ./scripts/generalized_processing/sample_features_script.py \
  $SAMPLE_ID \
  "$CENTROIDS_INTERSECT_BED" \
  "$FRAG_ENDS_INTERSECT_BED" \
  "$FRAG_ENDS_OCF_BED" \
  "$FRAG_WPS_INTERSECT_BED" \
    $MAPQ_FILTER \
    $RERUN \
    "${FEATURES_PATH:-./data/internal_features}"

echo "job finished successfully"
echo "deleting temp files"
rm ./data/sample_temp/${SAMPLE_ID}_autosomes_chr_id.bed
rm ./data/sample_temp/${SAMPLE_ID}_frag_centroids_openchrom_intersect.bed
rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_openchrom_intersect.bed
rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_ocf.bed
rm ./data/sample_temp/${SAMPLE_ID}_frag_wps_intersect.bed
rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_U_D.bed
rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region.bed
# rm ./data/sample_temp/${SAMPLE_ID}_frag_ends_with_region_sorted.bed

