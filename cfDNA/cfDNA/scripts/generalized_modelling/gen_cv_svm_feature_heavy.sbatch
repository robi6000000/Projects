#!/bin/bash
#SBATCH --job-name=cfdna_svm
#SBATCH --output=logs/svm_%x_%j.out
#SBATCH --error=logs/svm_%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --array=0-0

echo "Modelling SVM on heavier features using 64G, this svm is using standardization before CV (fsr is the largest matrix)"
cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA

#   $1 OPTION (cv,train,test)
#   $2 KERNEL
#   $3 GC_CORRECTION (true/false)
#   $4 PCA (true/false)
#   $5 PCA_COMPONENTS
#   $6 CV_REPEATS (default 10)
#   $7 MAPQ_FILTER (default None)
FEATURES=(fsr)
FEATURE=${FEATURES[$SLURM_ARRAY_TASK_ID]}
OPTION=${1:-cv}
MATRIX_FOLDER=data/matrix
KERNEL=$2
GC_CORRECTION=$3
PCA=$4
PCA_COMPONENTS=$5
CV_REPEATS=${6:-10}
MAPQ_FILTER=${7:-None}

MATRIX=data/matrix/by_feature/matrix_${FEATURE}.csv
if [ "$MAPQ_FILTER" != "None" ]; then
    MATRIX=data/matrix/by_feature/matrix_${FEATURE}_mapq${MAPQ_FILTER}.csv
fi
GC=data/processing/gc_content_per_region.csv


mkdir -p logs

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate cfdna
export PYTHONWARNINGS="once:Mean of empty slice:RuntimeWarning,once:Degrees of freedom <= 0 for slice.:RuntimeWarning"
python -u scripts/generalized_modelling/gen_cv_svm_feature.py $OPTION $MATRIX $MATRIX_FOLDER $GC $FEATURE $KERNEL $GC_CORRECTION $PCA $PCA_COMPONENTS $CV_REPEATS $MAPQ_FILTER