# /data/projects/liquid_biopsy/run_finaledb_gs.sh
#!/bin/bash
#SBATCH --job-name=snakemake_finaledb
#SBATCH --output=/data/projects/liquid_biopsy/slurm_logs/snakemake_master.%j.log
#SBATCH --error=/data/projects/liquid_biopsy/slurm_logs/snakemake_master.%j.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=RM-shared

# Submit this controller job with sbatch. It runs Snakemake, and Snakemake
# submits the per-rule worker jobs to SLURM using the repo-local profile.
#
# The worker-job partition is controlled by PARTITION below.

source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

cd /data/projects/liquid_biopsy

RUN_DIR="${RUN_DIR:-/data/projects/liquid_biopsy}"
PROFILE_DIR="${PROFILE_DIR:-/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/finaledb_workflow/slurm/profile}"
PARTITION="${PARTITION:-RM-shared}"

mkdir -p "$RUN_DIR/slurm_logs"
export PARTITION

TARGETS=$(cat gs_targets.txt | tr '\n' ' ')

snakemake \
  --profile "$PROFILE_DIR" \
  --directory "$RUN_DIR" \
  -s finaledb_workflow/slurm/workflow.slurm.smk \
  --reason \
  --rerun-incomplete \
  $TARGETS