#!/usr/bin/env bash
set -euo pipefail

# script to download a sample of the cristiano dataset from the manifest file
MANIFEST=data/manifest/Cristiano_manifest.csv
# MANIFEST=data/manifest/cristiano_manifest.csv     # old
OUTDIR=data/source/cristiano/
OUTDIR_WPS=data/source/cristiano_features/cristiano_wps_download/
INDEX=$1   
BASE=https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/

line=$(tail -n +2 "$MANIFEST" | sed -n "${INDEX}p")
IFS=, read -r seqrun_id sample_name age gender tissue disease read_length instrument frag_num frag_tsv_hg19 wps_bw_hg19 stage <<< "$line"
filename=$(basename "$frag_tsv_hg19")
filename_wps=$(basename "$wps_bw_hg19")
outfile="$OUTDIR/$filename"
outfile_wps="$OUTDIR_WPS/$filename_wps"
url="$BASE$frag_tsv_hg19"
url_wps="$BASE$wps_bw_hg19"

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR_WPS"
mkdir -p logs

# download fragment tsv
if [ -f "$outfile" ]; then
    echo "frag file already exists: $filename"
else
    wget -c -O "$outfile" "$url"
fi

# download WPS bigwig
if [ -f "$outfile_wps" ]; then
    echo "WPS file already exists: $filename_wps"
else
    wget -c -O "$outfile_wps" "$url_wps"
fi

# wget -c -O data/source/cristiano_features/cristiano_wps_download/ https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/entries/EE87786/hg19/EE87786.hg19.wps.mapq30.bw

