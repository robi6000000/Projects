#!/bin/bash
#SBATCH --job-name=sig_wps
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mem=24G

echo "Single-feature generalized SVM run (cv/train/test)"
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

# Unified sbatch argument order:
#   $1 OPTION (cv|train|test)
#   $2 MATRIX_FOLDER
#   $3 FEATURE
#   $4 KERNEL
#   $5 GC_CORRECTION (true|false)
#   $6 PCA (true|false)
#   $7 PCA_COMPONENTS
#   $8 CV_REPEATS
#   $9 MAPQ_FILTER
#   $10 RUN_TAG
OPTION=$1
MATRIX_FOLDER=$2
FEATURE=$3
KERNEL=$4
GC_CORRECTION=$5
PCA=$6
PCA_COMPONENTS=$7
CV_REPEATS=$8
MAPQ_FILTER=$9
RUN_TAG=${10}

MATRIX=$MATRIX_FOLDER/by_feature/matrix_${FEATURE}.csv
if [[ "${MAPQ_FILTER,,}" != "none" && "$MAPQ_FILTER" != "0" ]]; then
    MATRIX=$MATRIX_FOLDER/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi
GC=data/processing/gc_content_per_region.csv


mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONWARNINGS="once:Mean of empty slice:RuntimeWarning,once:Degrees of freedom <= 0 for slice.:RuntimeWarning"
python -u scripts/generalized_modelling/gen_svm_feature.py $OPTION $MATRIX $MATRIX_FOLDER $GC \
                                                             $FEATURE $KERNEL $GC_CORRECTION \
                                                             $PCA $PCA_COMPONENTS $CV_REPEATS $MAPQ_FILTER $RUN_TAG