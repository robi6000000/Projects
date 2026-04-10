cat > run_finaledb.sh << 'EOF'
#!/bin/bash

export TMPDIR=/data/users/kenderes/tmp
mkdir -p slurm_logs

# input path to file with target file names, in the formm "frag/samplename.GRCh37.frag.bed.gz"
TARGETS_PATH=$1
TARGETS=$(cat $TARGETS_PATH)

snakemake \
  --profile finaledb_workflow/slurm/profile \
  -s finaledb_workflow/slurm/workflow.slurm.smk \
  --reason \
  $TARGETS
EOF
chmod +x run_finaledb.sh