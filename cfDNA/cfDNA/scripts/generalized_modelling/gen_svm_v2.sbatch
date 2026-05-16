#!/bin/bash
#SBATCH --job-name=cfdna_svm_v2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mem=24G

#   $1 OPTION (cv, train)
#   $2 MATRIX_FOLDER
#   $3 KERNEL
#   $4 GC_CORRECTION (true/false)
#   $5 PCA (true/false)
#   $6 PCA_COMPONENTS
#   $7 CV_REPEATS
#   $8 MAPQ_FILTER
#   $9 RUN_TAG (optional) — also feeds into the auto-derived config name; pass e.g. "v2" to keep outputs separate from v1.

FEATURES=(length edm iedm eedm eoedm cposedm fsd pfe coverage ends ocf ifs wps fsr)
FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}
OPTION=$1
MATRIX_FOLDER=$2
KERNEL=$3
GC_CORRECTION=$4
PCA=$5
PCA_COMPONENTS=$6
CV_REPEATS=$7
MAPQ_FILTER=$8
RUN_TAG=$9

MATRIX=$MATRIX_FOLDER/by_feature/matrix_${FEATURE}.csv
if [[ "${MAPQ_FILTER,,}" != "none" && "$MAPQ_FILTER" != "0" ]]; then
    MATRIX=$MATRIX_FOLDER/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi
GC=data/processing/gc_content_per_region.csv

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONWARNINGS="once:Mean of empty slice:RuntimeWarning,once:Degrees of freedom <= 0 for slice.:RuntimeWarning"
python -u scripts/generalized_modelling/gen_svm_feature_v2.py $OPTION $MATRIX $MATRIX_FOLDER $GC $FEATURE $KERNEL $GC_CORRECTION $PCA $PCA_COMPONENTS $CV_REPEATS $MAPQ_FILTER $RUN_TAG
