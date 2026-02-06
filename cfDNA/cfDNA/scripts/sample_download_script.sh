#!/usr/bin/env bash
set -euo pipefail

MANIFEST=data/manifest/cristiano_manifest.csv
OUTDIR=data/source/cristiano/
INDEX=$1   
BASE=https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/

line=$(tail -n +2 "$MANIFEST" | sed -n "${INDEX}p")
IFS=, read -r seqrun sample disease tissue readlen instrument file_path <<< "$line"

filename=$(basename "$file_path")
outfile="$OUTDIR/$filename"
url="$BASE$file_path"
mkdir -p "$OUTDIR"
# mkdir -p "logs"

if [ -f "$outfile" ]; then
    echo "Already exists, skipping: $filename"
    exit 0
fi

echo "Downloading: $filename"
wget -c -O "$outfile" "$url"