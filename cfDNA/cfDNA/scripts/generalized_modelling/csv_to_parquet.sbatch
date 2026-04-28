#!/bin/bash
#SBATCH --job-name=csv_to_parquet
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# Convert a single feature matrix CSV to Parquet (float32).
# Usage:
#   sbatch scripts/generalized_modelling/csv_to_parquet.sbatch <csv_path>
#
# Example:
#   sbatch scripts/generalized_modelling/csv_to_parquet.sbatch data/matrix/by_feature/matrix_pfe_mapq30.csv

CSV_PATH=$1

if [ -z "$CSV_PATH" ]; then
    echo "Error: no CSV path provided."
    echo "Usage: sbatch csv_to_parquet.sbatch <csv_path>"
    exit 1
fi

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

python -u scripts/generalized_modelling/csv_to_parquet.py $CSV_PATH
