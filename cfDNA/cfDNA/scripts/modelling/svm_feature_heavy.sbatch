#!/bin/bash
#SBATCH --job-name=cfdna_svm
#SBATCH --output=logs/svm_%x_%j.out
#SBATCH --error=logs/svm_%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --array=0-0

echo "Modelling SVM on heavier features using 48G, this svm is using standardization before CV"
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# FEATURE=$1
# KERNEL=$2
# GC_CORRECTION=$3
# PCA=$4
# CV_REPEATS=$5 (optional, default 10)
FEATURES=(fsr)
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