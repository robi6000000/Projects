conda create -n snakemake python=3.12 snakemake=7.32.4 -c conda-forge -c bioconda
conda activate snakemake
conda install -c conda-forge -c bioconda \
  bwa samtools samblaster bedtools bioawk htslib trimmomatic

<!-- dry run: -->
snakemake   -s finaledb_workflow/slurm/workflow.slurm.smk   -n -p   frag/gsca_pca_00001_01_01_pl.GRCh37.frag.bed.gz                                

<!-- added the lines: -->
wildcard_constraints:
    ref_genome="hg19|hg38|hg19_mm10"

<!-- run: -->
snakemake --cores 16 \
  -s finaledb_workflow/slurm/workflow.slurm.smk -p \
  frag/gsca_pca_00001_01_01_pl.GRCh37.frag.bed.gz

number of sample files (these are pairs):
kenderes@gen-manager:/data/projects/liquid_biopsy$ find -L reads/original -maxdepth 1 -type f -prin
tf '%f\n' | cut -c1-4 | sort | uniq -c | sort -nr
    954 lypp
    438 lycc
    114 gspp
    114 gsca
      1 ln.s

<!-- slurm setup -->
We suggest invoking snakemake workflows with --profile (https://snakemake.readthedocs.io/en/stable/executing/cli.html?highlight=profile#profiles). For more information, please refer to this excellent blog post. Below is a profile configiuration for example:
'''jobs: 1000
cluster: "sbatch -t 4320 -p {params.partition} --mem={resources.mem_mb} --time-min {resources.time_min} --ntasks-per-node {resources.cpus} --job-name snakemake.{rule}.{params.label} -o slurm_logs/slurm.jobid_%j.{rule}.{params.label}.log"
default-resources: [cpus=1, time_min=1440, mem_mb=4096]
singularity-prefix: "direcotry for keeping singularity images"
singularity-args: "--bind $LOCAL:/local"
directory: "your working directory"
keep-going: true'''

<!-- creating symlink folder -->
<!-- the slurm finaledb workflow needs fq.gz format so we also need to rename the files -->
for f in /data/projects/liquid_biopsy/reads/original/gs*.fastq.gz; do
    base=$(basename "$f" .fastq.gz)
    ln -s "$f" "fastq/${base}.fq.gz"
done

