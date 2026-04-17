#!/bin/bash
#SBATCH --job-name=cfdna_svm
#SBATCH --output=logs/svm_%x_%j.out
#SBATCH --error=logs/svm_%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --array=0-0

echo "Modelling SVM on heavier features using 64G"
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

#   $1 OPTION (cv,train,test)
#   $2 MATRIX_FOLDER
#   $3 KERNEL
#   $4 GC_CORRECTION (true/false)
#   $5 PCA (true/false)
#   $6 PCA_COMPONENTS
#   $7 CV_REPEATS (optional, default 10)
#   $8 MAPQ_FILTER (optional, default None)

# usecase: sbatch scripts/generalized_modelling/gen_svm_feature_heavy.sbatch cv data/matrix rbf false false 150 1 30
FEATURES=(fsr)
FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}
OPTION=${1:-cv}
MATRIX_FOLDER=$2
KERNEL=$3
GC_CORRECTION=$4
PCA=$5
PCA_COMPONENTS=$6
CV_REPEATS=${7:-10}
MAPQ_FILTER=${8:-None}

MATRIX=data/matrix/by_feature/matrix_${FEATURE}.csv
if [ "$MAPQ_FILTER" != "None" ]; then
    MATRIX=data/matrix/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi
GC=data/processing/gc_content_per_region.csv


mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONWARNINGS="once:Mean of empty slice:RuntimeWarning,once:Degrees of freedom <= 0 for slice.:RuntimeWarning"
python -u scripts/generalized_modelling/gen_svm_feature.py $OPTION $MATRIX $MATRIX_FOLDER $GC $FEATURE $KERNEL $GC_CORRECTION $PCA $PCA_COMPONENTS $CV_REPEATS $MAPQ_FILTER