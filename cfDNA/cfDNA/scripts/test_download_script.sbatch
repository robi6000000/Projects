#!/bin/bash
#SBATCH --job-name=sample_test
#SBATCH --output=logs/download_test_%j.out
#SBATCH --error=logs/download_test_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=1G

set -euo pipefail

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

MANIFEST=data/manifest/cristiano_manifest.csv
OUTDIR=data/source/cristiano/
INDEX=2   
BASE=https://finaledb.research.cchmc.org/

mkdir -p "$OUTDIR"

url=$(sed -n "${INDEX}p" "$MANIFEST" | cut -d, -f7)
url="$BASE$url"

echo "Downloading:"
echo "$url"

cd "$OUTDIR"
wget -c "$url"
