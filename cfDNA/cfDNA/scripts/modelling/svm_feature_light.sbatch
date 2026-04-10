#!/bin/bash
#SBATCH --job-name=cfdna_svm
#SBATCH --output=logs/svm_%x_%j.out
#SBATCH --error=logs/svm_%x_%j.err
#SBATCH --mem=24G
#SBATCH --array=0-11%4

echo "Modelling SVM on lighter features using 24G, this svm is using standardization and pca inside the CV loop"
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# FEATURE=$1
# KERNEL=$2
# GC_CORRECTION=$3
# PCA=$4
# CV_REPEATS=$5 (optional, default 10)
FEATURES=(length edm poem_prem ext_poem_prem fsd pfe coverage ends ocf ifs wps wps_compute)
FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}
KERNEL=$1
GC_CORRECTION=$2
PCA=$3
PCA_COMPONENTS=$4
CV_REPEATS=${5:-10}
MAPQ_FILTER=${6:-None}

MATRIX=data/matrix/by_feature/matrix_${FEATURE}.csv
if [ "$MAPQ_FILTER" != "None" ]; then
    MATRIX=data/matrix/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi
GC=data/processing/gc_content_per_region.csv


mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

python -u scripts/modelling/svm_feature.py $MATRIX $GC $FEATURE $KERNEL $GC_CORRECTION $PCA $PCA_COMPONENTS $CV_REPEATS $MAPQ_FILTER