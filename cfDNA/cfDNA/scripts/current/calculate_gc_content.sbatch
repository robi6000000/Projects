#!/bin/bash
#SBATCH --job-name=gc_content
#SBATCH --output=logs/gc_content_%j.out
#SBATCH --error=logs/gc_content_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# Edit these paths as needed
REGIONS=./data/processing/openchrom_with_id.bed
OUTPUT=./data/processing/gc_content_per_region.csv
GENOME=./data/hg19/hg19.fa

python scripts/modelling/calculate_gc_content.py "$REGIONS" "$OUTPUT" "$GENOME"

echo "GC content calculation complete"
