#!/bin/bash
#SBATCH --job-name=cfDNA_ref
#SBATCH --output=logs/ref_%j.out
#SBATCH --error=logs/ref_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
bash scripts/static_script.sh
