#!/bin/bash
#SBATCH --job-name=gc_content
#SBATCH --output=logs/gc_content_%j.out
#SBATCH --error=logs/gc_content_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

python scripts/modelling/calculate_gc_content.py

echo "GC content calculation complete"