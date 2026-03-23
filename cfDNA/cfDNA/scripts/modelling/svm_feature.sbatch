#!/bin/bash
#SBATCH --job-name=cfdna_svm
#SBATCH --output=logs/svm_%x_%j.out
#SBATCH --error=logs/svm_%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

FEATURE=$1
KERNEL=$2
GC_CORRECTION=$3
PCA=$4

MATRIX=data/matrix/by_feature/matrix_${FEATURE}.csv
GC=data/processing/gc_content_per_region.csv


# mkdir -p data/matrix/svm_by_feature
mkdir -p data/matrix/svm_linear_by_feature
mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

python -u scripts/modelling/svm_feature.py $MATRIX $GC $FEATURE $KERNEL $GC_CORRECTION $PCA