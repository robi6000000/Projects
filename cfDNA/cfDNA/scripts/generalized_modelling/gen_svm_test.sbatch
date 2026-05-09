#!/bin/bash
#SBATCH --job-name=cfdna_svm_test
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mem=64G

#   $1 TEST_MATRIX_FOLDER  (e.g. data/matrix_internal_gs)
#   $2 MODEL_CONFIG_DIR    (e.g. data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30)
#   $3 DATASET_DIR         (e.g. cristiano, ly_test) — output subfolder under test/
#   $4 MAPQ_FILTER         (e.g. 30, none)

# FEATURES index reference:
#   light (8G)  -- 0-6:  length edm iedm eedm eoedm cposedm fsd
#   med   (16G) -- 7-12: pfe coverage ends ocf ifs wps
#   heavy (32G) -- 13:   fsr
FEATURES=(length edm iedm eedm eoedm cposedm fsd pfe coverage ends ocf ifs wps fsr)
FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}

TEST_MATRIX_FOLDER=$1
MODEL_CONFIG_DIR=$2
DATASET_DIR=$3
MAPQ_FILTER=$4

CONFIG_NAME=$(basename "$MODEL_CONFIG_DIR")

TEST_MATRIX=$TEST_MATRIX_FOLDER/by_feature/matrix_${FEATURE}.csv
if [[ "${MAPQ_FILTER,,}" != "none" && "$MAPQ_FILTER" != "0" ]]; then
    TEST_MATRIX=$TEST_MATRIX_FOLDER/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi

MODEL_PATH=$MODEL_CONFIG_DIR/models/${CONFIG_NAME}_${FEATURE}.pkl
OUTPUT_PATH=$MODEL_CONFIG_DIR/test/${DATASET_DIR}/${FEATURE}_probs.csv

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONWARNINGS="once:Mean of empty slice:RuntimeWarning,once:Degrees of freedom <= 0 for slice.:RuntimeWarning"

python -u scripts/generalized_modelling/gen_svm_test.py "$TEST_MATRIX" "$MODEL_PATH" "$OUTPUT_PATH"
