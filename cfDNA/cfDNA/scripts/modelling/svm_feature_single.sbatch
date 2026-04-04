#!/bin/bash
#SBATCH --job-name=svm
#SBATCH --output=logs/svm_%x_%j.out
#SBATCH --error=logs/svm_%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

FEATURE=$1
KERNEL=$2
GC_CORRECTION=$3
PCA=$4
PCA_COMPONENTS=${5:-0.95}
CV_REPEATS=${6:-10}
MAPQ_FILTER=${7:-None}
MATRIX=data/matrix/by_feature/matrix_${FEATURE}.csv
if [ "$MAPQ_FILTER" != "None" ]; then
    MATRIX=data/matrix/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi
GC=data/processing/gc_content_per_region.csv

python -u scripts/modelling/svm_feature.py $MATRIX $GC $FEATURE $KERNEL $GC_CORRECTION $PCA $PCA_COMPONENTS $CV_REPEATS $MAPQ_FILTER  