#!/bin/bash
#SBATCH --job-name=ensemble_svm
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mem=4G

# usecase for cv:
# sbatch gen_ensemble_svm.sbatch cv <train_meta_matrix.csv> linear 10 none none none nostd false
# usecase for train:
# sbatch gen_ensemble_svm.sbatch train <train_meta_matrix.csv> linear 10 none none none nostd false
# usecase for test (NEEDS dataset_tag and trained model):
# sbatch gen_ensemble_svm.sbatch test <test_meta_matrix.csv> linear 10 <dataset_tag> <trained_model.pkl> none nostd false

#   $1 OPTION           (cv/train/test)
#   $2 META_MATRIX_PATH
#   $3 KERNEL
#   $4 CV_REPEATS
#   $5 DATASET_TAG      (none if not applicable, e.g. cristiano)
#   $6 MODEL_PATH       (none to use default ensemble/models/ensemble.pkl)
#   $7 TEST_OUTPUT_PATH (none to use default ensemble/test/<dataset_tag>_preds.csv)
#   $8 RUN_TAG          (none, or suffix for ensemble dir, e.g. v2 -> ensemble_v2/)
#   $9 STANDARDIZE      (optional, true|false, default false)

cd /data/projects/liquid_biopsy/Projects/cfDNA_v2
mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna

python -u scripts/generalized_modelling/gen_ensemble_svm.py "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "${9:-true}"
