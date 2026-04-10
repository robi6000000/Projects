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
    read=$(echo "$base" | grep -o 'R[12]$')
    sample="${base%_${read}}"
    if [[ "$sample" == *"_pl" ]]; then
        ln -s "$f" "fastq/${sample}.${read}.fq.gz"
    fi
done

delete frag/ files which are not pl or are not gs_tartget_names.txt:
find frag/ -type f ! -name "*_pl*" ! -name "gs_target_names.txt" -delete


```Considering that feature dimensions of different
fragmentation patterns vary considerably, we performed a principal
component analysis dimensionality reduction on all fragmentation
patterns to the same dimensions before constructing
the classificationmodel.Model performances before and after dimensionality
reduction were found to be similar (Table 1).```

sed -i 's|HTTP.remote("https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa")|"/gen-nas01/active_projects/gen-manager/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/supplementary/TruSeq3-PE-2.fa"|' \
  finaledb_workflow/slurm/workflow.slurm.smk

snakemake \
  --profile finaledb_workflow/slurm/profile \
  -s finaledb_workflow/slurm/workflow.slurm.smk \
  --reason \
  frag/gsca_pca_00001_01_01_pl.GRCh37.frag.bed.gz


  <!-- generate file with all target names -->
ls fastq/gs*.R1.fq.gz | while read f; do
    base=$(basename "$f" .R1.fq.gz)
    echo "frag/${base}.GRCh37.frag.bed.gz"
done > frag/gs_target_names.txt

cat frag/gs_target_names.txt | head -5
wc -l frag/gs_target_names.txt  

<!-- activate screen env -->
screen -S lynch_preprocess

snakemake \
  --profile finaledb_workflow/slurm/profile \
  -s finaledb_workflow/slurm/workflow.slurm.smk \
  --reason \
  $(cat frag/ly_target_names.txt)
