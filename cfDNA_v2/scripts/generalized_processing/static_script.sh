#!/usr/bin/env bash
set -euo pipefail

mkdir -p ./data/processing
mkdir -p ./data/sample_temp
mkdir -p ./data/rows_sample_temp
mkdir -p ./data/source

echo "Getting liftover, ATAC-seq, and DNase-seq files"
# download liftover file
if [ ! -f "./data/source/hg38ToHg19.over.chain" ]; then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O ./data/source/hg38ToHg19.over.chain.gz
    gunzip ./data/source/hg38ToHg19.over.chain.gz
fi

# WGET COMMANDS FOR ATAC AND DNASE peaks
if [ ! -f "./data/source/E029-DNase.hotspot.fdr0.01.broad.bed.gz" ]; then
    wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E029-DNase.hotspot.fdr0.01.broad.bed.gz -O ./data/source/E029-DNase.hotspot.fdr0.01.broad.bed.gz
fi
if [ ! -f "./data/source/E032-DNase.hotspot.fdr0.01.broad.bed.gz" ]; then
    wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E032-DNase.hotspot.fdr0.01.broad.bed.gz -O ./data/source/E032-DNase.hotspot.fdr0.01.broad.bed.gz
fi
if [ ! -f "./data/source/E034-DNase.hotspot.fdr0.01.broad.bed.gz" ]; then
    wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E034-DNase.hotspot.fdr0.01.broad.bed.gz -O ./data/source/E034-DNase.hotspot.fdr0.01.broad.bed.gz
fi
if [ ! -f "./data/source/TCGA-ATAC_PanCancer_PeakSet.txt" ]; then
    wget https://api.gdc.cancer.gov/data/116ebba2-d284-485b-9121-faf73ce0a4ec -O ./data/source/TCGA-ATAC_PanCancer_PeakSet.txt
fi    
echo "Liftover and ATAC-seq, DNase-seq files done"

# pancancer peaks file has header, which is skipped, and select only the 3 relevant columns
tail -n +2 ./data/source/TCGA-ATAC_PanCancer_PeakSet.txt | cut -f1-3 > ./data/source/TCGA-ATAC_PanCancer_PeakSet.noheader.txt

echo "Performing liftover of pancancer peaks from hg38 to hg19"
# liftover pancancer peaks
if [ ! -f "liftOver" ]; then
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O liftOver
    chmod +x liftOver
fi
./liftOver \
  ./data/source/TCGA-ATAC_PanCancer_PeakSet.noheader.txt \
  ./data/source/hg38ToHg19.over.chain \
  ./data/processing/mat.pancancer.hg19.bed \
  ./data/processing/mat.unmapped.hg19.bed
echo "Liftover complete"

echo "creating open chromatin region pool"
# pipe version of processing dna and atacseq
( zcat ./data/source/E029-DNase.hotspot.fdr0.01.broad.bed.gz | cut -f1-3 ; 
  zcat ./data/source/E032-DNase.hotspot.fdr0.01.broad.bed.gz | cut -f1-3 ;
  zcat ./data/source/E034-DNase.hotspot.fdr0.01.broad.bed.gz | cut -f1-3 ;
  cat ./data/processing/mat.pancancer.hg19.bed ) | \
  sort -k1,1V -k2,2n | \
  bedtools merge | \
  awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' | \
#   compute centroid and adjust regions to 200bp centered on the centroid
  awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); start=c-100; if(start<0)start=0; end=c+100; print $1,start,end}' | \
#   add region ids, $0 means the whole line, NR is the line number
  awk 'BEGIN{OFS="\t"} {print $0, NR-1}' \
  > ./data/processing/openchrom_with_id.bed 
echo "Open chromatin regions file created - ./data/processing/openchrom_with_id.bed"


# TSS processing — not used by current pipeline features, kept for reference
# # get reference genome chromosome sizes
# if [ ! -f "./data/source/hg19.genome" ]; then
#     wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O ./data/source/hg19.genome
# fi
# if [ ! -f "./data/source/gencode.v30lift37.annotation.gtf.gz" ] && [ ! -f "./data/source/gencode.v30lift37.annotation.gtf" ]; then
#     wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gtf.gz -O ./data/source/gencode.v30lift37.annotation.gtf.gz
# fi
# if [ ! -f "./data/source/gencode.v30lift37.annotation.gtf" ]; then
#     gunzip ./data/source/gencode.v30lift37.annotation.gtf.gz
# fi
# awk 'BEGIN{OFS="\t"}
# $3=="transcript" {
#     if ($7 == "+") print $1, $4, $4, ".", ".", $7;
#     else if ($7 == "-") print $1, $5, $5, ".", ".", $7;
# }' ./data/source/gencode.v30lift37.annotation.gtf | \
# awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' > ./data/processing/tss_autosomes.bed
# bedtools slop -i ./data/processing/tss_autosomes.bed -g ./data/source/hg19.genome -l 150 -r 50 -s > ./data/processing/tss_150_50.bed
# bedtools slop -i ./data/processing/tss_autosomes.bed -g ./data/source/hg19.genome -l 1000 -r 1000 -s > ./data/processing/tss_1000_1000.bed

# for OCF we also need fragment end information intersected with openchrom
# add a centroid column
awk 'BEGIN{OFS="\t"} {
    c=int(($2+$3)/2);
    print $1, $2, $3, c, $4
}' ./data/processing/openchrom_with_id.bed \
> ./data/processing/openchrom_with_id_centroid.bed

# extend openchrom regions by 60bp on each side for WPS computation
awk 'BEGIN{OFS="\t"} {
    chrom = $1
    start = $2
    end = $3
    region_id = $4

    wps_start = start - 60
    if (wps_start < 0) wps_start = 0
    wps_end = end + 60
    print chrom, wps_start, wps_end, region_id
}' ./data/processing/openchrom_with_id.bed > ./data/processing/wps_openchrom_with_id.bed

# get reference genome fasta for end motif:
if [ ! -f "./data/hg19/hg19.fa.gz" ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O ./data/hg19/hg19.fa.gz
fi
if [ ! -f "./data/hg19/hg19.fa" ]; then
    gunzip ./data/hg19/hg19.fa.gz
fi

# get v37 reference genome fasta for finaledb pipeline:
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
echo "Downloading reference genome for finaledb pipeline"
if [ ! -f "./data/hg19/human_g1k_v37.fasta.gz" ] && [ ! -f "./data/hg19/human_g1k_v37.fasta" ]; then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -O ./data/hg19/human_g1k_v37.fasta.gz
fi
if [ ! -f "./data/hg19/human_g1k_v37.fasta" ]; then
    gunzip ./data/hg19/human_g1k_v37.fasta.gz
fi
echo "Reference genome download complete"


# Run GC content calculation 
echo "Calculating GC content per region"
python scripts/generalized_processing/calculate_gc_content.py ./data/processing/openchrom_with_id.bed ./data/processing/gc_content_per_region.csv ./data/hg19/hg19.fa
echo "GC content calculation complete"