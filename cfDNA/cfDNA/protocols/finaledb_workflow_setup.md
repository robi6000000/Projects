# finaledb workflow — local setup notes

## Source and isolation strategy

The finaledb workflow was originally a separate git repository (submodule). To avoid pushing local changes back to the upstream finaledb repo, the `.git` directory was removed and the files vendored directly into this project repo (`robi6000000/Projects`). The original `.git` folder was preserved as `.git.backup` inside `finaledb_workflow/` in case the upstream remote config is ever needed.

Commit: `879d63b` — "vendor finaledb workflow into repo" (2026-04-04)

---

## Changes made to workflow.slurm.smk

### 1. Trimmomatic adapter — replace HTTP remote with local file path
**Commit:** `02dcbd1` — "update 10_4_2026"

The original workflow fetched the Trimmomatic adapter file at runtime over HTTP:
```python
# original
adapter=HTTP.remote("https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa"),
```

This was replaced with a local absolute path (the file was downloaded once and placed in `supplementary/`):
```python
# changed
adapter="/gen-nas01/active_projects/gen-manager/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/supplementary/TruSeq3-PE-2.fa",
```

**Why:** HTTP remote fetching was unreliable in the SLURM environment and introduced an external dependency at job runtime.

### 2. File permissions
**Commit:** `124aead` — file mode changed from `100644` → `100755` (made executable). No content change.

---

## MAPQ filtering — workflow does NOT filter by mapping quality

The fragment extraction step uses:
```bash
samtools view -h -f 3 -F 3852 -G 48 --incl-flags 48 {input.bam}
```

There is **no `-q` flag**. All properly-paired, non-duplicate, primary alignments are written to the `.frag.bed.gz` file regardless of mapping quality. The mapq score is stored in column 5 of the fragment BED (`"#chrom", "start", "end", "name", "mapq", "strand"`).

MAPQ filtering is applied post-hoc in `compute_features_mapq.sbatch` when computing features (parameter 1: `mapq_filter`). The differences between mapq0/15/30/45 CV results are real — there is no hidden floor baked into the frag files.

---

## How to run the workflow

Activate the screen session and run snakemake with the slurm profile:

```bash
screen -S lynch_preprocess

snakemake \
  --profile finaledb_workflow/slurm/profile \
  -s finaledb_workflow/slurm/workflow.slurm.smk \
  --reason \
  $(cat frag/gs_target_names.txt)
```

For a single test sample:
```bash
snakemake \
  --profile finaledb_workflow/slurm/profile \
  -s finaledb_workflow/slurm/workflow.slurm.smk \
  --reason \
  frag/gsca_pca_00001_01_01_pl.GRCh37.frag.bed.gz
```

### Preparing input fastq symlinks
The workflow expects `fastq/{sample}.R1.fq.gz` / `fastq/{sample}.R2.fq.gz`. Symlinks from the original `.fastq.gz` files:
```bash
for f in /data/projects/liquid_biopsy/reads/original/gs*.fastq.gz; do
    base=$(basename "$f" .fastq.gz)
    read=$(echo "$base" | grep -o 'R[12]$')
    sample="${base%_${read}}"
    if [[ "$sample" == *"_pl" ]]; then
        ln -s "$f" "fastq/${sample}.${read}.fq.gz"
    fi
done
```

### Generating the target names file
```bash
ls fastq/gs*.R1.fq.gz | while read f; do
    base=$(basename "$f" .R1.fq.gz)
    echo "frag/${base}.GRCh37.frag.bed.gz"
done > frag/gs_target_names.txt
```

### Cleaning non-target frag files
```bash
find frag/ -type f ! -name "*_pl*" ! -name "gs_target_names.txt" -delete
```
