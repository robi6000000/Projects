#!/bin/bash
#SBATCH --job-name=cfdna_svm_test
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mem=24G

# $1 TEST_MATRIX_PATH
# $2 MODEL_PATH
# $3 OUTPUT_PROBS_PATH

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONWARNINGS="once:Mean of empty slice:RuntimeWarning,once:Degrees of freedom <= 0 for slice.:RuntimeWarning"

python -u scripts/generalized_modelling/gen_svm_test.py "$1" "$2" "$3"
