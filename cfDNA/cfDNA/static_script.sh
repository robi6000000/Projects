#!/usr/bin/env bash


# liftover pancancer peaks from hg38 to hg19
# download liftover file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O data/source_data/hg38ToHg19.over.chain.gz
gunzip data/source_data/hg38ToHg19.over.chain.gz

# TODO ADD WGET COMMANDS FOR ATAC AND DNAS SeTS
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E029-DNase.hotspot.fdr0.01.broad.bed.gz -O data/source_data/E029-DNase.hotspot.fdr0.01.broad.bed.gz
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E032-DNase.hotspot.fdr0.01.broad.bed.gz -O data/source_data/E032-DNase.hotspot.fdr0.01.broad.bed.gz
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E034-DNase.hotspot.fdr0.01.broad.bed.gz -O data/source_data/E034-DNase.hotspot.fdr0.01.broad.bed.gz
wget https://api.gdc.cancer.gov/data/116ebba2-d284-485b-9121-faf73ce0a4ec -O data/source_data/TCGA-ATAC_PanCancer_PeakSet.txt

# remove header from pancancer peaks file and select only thee 3 columns
tail -n +2 data/source_data/TCGA-ATAC_PanCancer_PeakSet.txt | cut -f1-3 > data/source_data/TCGA-ATAC_PanCancer_PeakSet.noheader.txt

# liftover pancancer peaks
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O liftOver
chmod +x liftOver
liftOver \
  data/source_data/TCGA-ATAC_PanCancer_PeakSet.noheader.txt \
  data/source_data/hg38ToHg19.over.chain \
  data/processing/mat.pancancer.hg19.bed \
  data/processing/mat.unmapped.hg19.bed

# dnaseq files - unzip and select only 3 columns (pancancer also has header)
( zcat data/source_data/E029-DNase.hotspot.fdr0.01.broad.bed.gz | cut -f1-3 ; 
  zcat data/source_data/E032-DNase.hotspot.fdr0.01.broad.bed.gz | cut -f1-3 ;
  zcat data/source_data/E034-DNase.hotspot.fdr0.01.broad.bed.gz | cut -f1-3 ;
  cat data/processing/mat.pancancer.hg19.bed ) \
  > data/processing/all_openchrom_regions.bed
# sort and merge to get open chromatin union
sort -k1,1V -k2,2n data/processing/all_openchrom_regions.bed > data/processing/all_openchrom_regions.sorted.bed
bedtools merge -i data/processing/all_openchrom_regions.sorted.bed > data/processing/openchrom_union.bed

# Autosomes = chromosomes 1–22 (exclude x,y) etc
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' data/processing/openchrom_union.bed > data/processing/openchrom_union_autosomes.bed
wc -l data/processing/openchrom_union_autosomes.bed
# bin the open chromatin regions
# 200bp centered regions - calculate centroid and add 100bp each side
awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); s=c-100; if(s<0)s=0; e=c+100; print $1,s,e}' data/processing/openchrom_union_autosomes.bed > data/processing/openchrom_200bp.bed
# add region ids
awk 'BEGIN{OFS="\t"} {print $0, NR-1}' \
    data/processing/openchrom_200bp.bed \
    > data/processing/openchrom_with_id.bed

# print head
# head data/processing/openchrom_with_id.bed
# get reference genome chromosome sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O data/source_data/hg19.genome


# get TSS positions from gencode gtf
awk 'BEGIN{OFS="\t"} 
$3=="transcript" {
    if ($7 == "+") print $1, $4, $4, ".", ".", $7;
    else if ($7 == "-") print $1, $5, $5, ".", ".", $7;
}' data/source_data/gencode.v30lift37.annotation.gtf > data/processing/tss_raw.bed


# check chromosomes and filter:
cut -f1 data/processing/tss_raw.bed | sort | uniq -c
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' data/processing/tss_raw.bed > data/processing/tss_autosomes.bed


bedtools slop -i data/processing/tss_autosomes.bed -g data/source_data/hg19.genome -l 150 -r 50 -s > data/processing/tss_150_50.bed
bedtools slop -i data/processing/tss_autosomes.bed -g data/source_data/hg19.genome -l 1000 -r 1000 -s > data/processing/tss_1000_1000.bed

# get reference genome fasta for end motif:
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O data/source_data/hg19.fa.gz
gunzip data/source_data/hg19.fa.gz
mkdir data/fullsample

# here ends the static processing
