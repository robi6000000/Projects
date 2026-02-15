#!/bin/bash
#SBATCH --job-name=finaledb_pipeline_test
#SBATCH --output=logs/finaledb_%j.out
#SBATCH --error=logs/finaledb_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

cd /data/projects/liquid_biopsy/Projects/cfDNA/cfDNA
mkdir -p logs
conda activate snakemake


snakemake \
  -s finaledb_workflow/finaledb.smk \
  --use-conda \
  --cores 1 \
  frag/gsca_pca_00001_01_01_pl.hg19.frag.bed.gz