# cfDNA_v2 pipeline

Working directory: `/data/projects/liquid_biopsy/Projects/cfDNA_v2/`

## Setup

- Conda env `cfdna` (`source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh && conda activate cfdna`)
- SLURM cluster (gen-compute01–04)
- Tools on PATH: `wget`, `bedtools`, `liftOver`, `samtools`, `python` (pandas, numpy, scikit-learn, statsmodels, pysam)
- One-time prerequisite: `bash scripts/generalized_processing/static_script.sh` — downloads reference files and produces `data/processing/openchrom_with_id.bed` (~560k open chromatin regions) and `data/processing/gc_content_per_region.csv`
- Input: gzipped fragment BED files (chrom, start, end, MAPQ, strand). Run `scripts/snakemake_finaledb/` first if input is paired FASTQ.

## Directory structure

```
cfDNA_v2/
├── data/
│   ├── metadata/                       # sample manifest CSVs
│   │   ├── Cristiano_metadata.csv
│   │   ├── internal_metadata_ly_filtered_train.csv
│   │   ├── internal_metadata_ly_filtered_test.csv
│   │   └── internal_metadata_gs.csv
│   ├── processing/                     # static reference files (built once by static_script.sh)
│   │   ├── openchrom_with_id.bed
│   │   ├── openchrom_with_id_centroid.bed
│   │   ├── wps_openchrom_with_id.bed
│   │   └── gc_content_per_region.csv
│   ├── source/                         # downloaded reference files
│   ├── hg19/                           # reference genome FASTA
│   ├── sample_temp/                    # per-sample intermediate BED files (auto-deleted)
│   ├── cristiano_features/             # per-sample feature CSVs, Cristiano cohort
│   ├── internal_features_ly/           # per-sample feature CSVs, PreveLynch cohort
│   ├── internal_features_gs/           # per-sample feature CSVs, GenoScan cohort
│   ├── cristiano_runs/                 # matrices and model runs for Cristiano
│   │   ├── feature_matrices/           # per-feature matrices
│   │   │   └── matrix_{feature}_mapq{N}.csv
│   │   └── model_runs/
│   │       └── svm_{kernel}_{pca}_{gc}_{mapq}/
│   │           ├── cv/                 # per-feature CV probabilities
│   │           ├── models/             # per-feature trained .pkl files
│   │           ├── test/               # per-feature test predictions
│   │           │   └── {dataset}/
│   │           ├── meta_matrix/        # merged CV probabilities for ensemble input
│   │           └── ensemble_{tag}/     # ensemble model run
│   │               ├── cv/
│   │               ├── models/
│   │               └── test/
│   ├── internal_ly_train_runs/         # same structure as cristiano_runs, PreveLynch train
│   ├── internal_ly_test_runs/          # same structure, PreveLynch test
│   ├── internal_gs_runs/
│   └── crc_runs/
└── scripts/
    ├── generalized_processing/
    │   ├── static_script.sh / .sbatch  # one-time reference file builder
    │   ├── generate_metadata_file.py   # generate manifest CSV for internal cohorts
    │   ├── compute_features_mapq.sbatch# SLURM array: per-sample feature extraction
    │   ├── sample_prep_script.sh       # called by above; creates intersection BEDs
    │   ├── sample_features_script.py   # called by above; computes 14 features per sample
    │   ├── matrix_build_array.sbatch   # SLURM array: aggregate per-sample CSVs into matrices
    │   ├── matrix_build.py             # called by above; one feature per invocation
    │   └── meta_matrix_build.py        # merge per-feature CV probs into ensemble input
    └── generalized_modelling/
        ├── gen_svm.sbatch              # SLURM array: per-feature CV or training
        ├── gen_svm_feature.py
        ├── gen_svm_test.sbatch         # SLURM array: per-feature inference on test cohort
        ├── gen_svm_test.py
        ├── gen_ensemble_svm.sbatch     # SLURM: ensemble CV, training, or test inference
        ├── gen_ensemble_svm.py
        └── modules/gen_cfdna_model.py  # CFDNAModel class (SVM, GC correction, CV logic)
```

## File naming conventions

**Manifest CSV** — one row per sample, columns: `sample_id, disease, dataset, material, frag_path, stage, cancer_true`

- Auto-generated for internal cohorts: `python scripts/generalized_processing/generate_metadata_file.py ly` → `data/metadata/internal_metadata_ly.csv`
- After generation, manually filter/split into train/test and fill the `stage` column

**Fragment files (internal cohorts)**
`{dataset}_{disease}_{...optional}_{material}.GRCh37.frag.bed.gz`
- `dataset` suffix `cc`/`ca` ⇒ `cancer_true=1`; anything else ⇒ 0
- Only `material == pl` (plasma) is kept

**Feature matrix**: `{run_dir}/feature_matrices/matrix_{feature}_mapq{N}.csv`

**Config name** (auto-derived): `svm_{kernel}[_pca{N}][_gc][_cv{R}][_mapq{M}][_{tag}]` — `_cv{R}` omitted when `cv_repeats=10`

## Assumptions

- Fragment BEDs: autosomes only (chr1–chr22); `chr` prefix optional, normalised by `sample_prep_script.sh`
- MAPQ filter is applied at fragment-prep time — thread the same value through every step
- 14 features computed per sample. Feature array index reference (used by every sbatch array job):
  `0 length, 1 edm, 2 iedm, 3 eedm, 4 eoedm, 5 cposedm, 6 fsd, 7 pfe, 8 coverage, 9 ends, 10 ocf, 11 ifs, 12 wps, 13 fsr`
- GC correction is only meaningful for `pfe, coverage, ends, ocf, ifs, wps`; `gen_svm_feature.py` auto-disables it for others with a warning
- Memory groups (empirical at MAPQ 30, ~459 samples):
  - light (0–6): 16G for CV; 8G for test
  - med (7–12): 32G for CV/train; 16G for test
  - heavy (13, fsr ~1.7M cols): 64G for CV/train; 40G for test
- Ensemble runs are cheap (4G) — operates on a 14-column meta-matrix

## Best config

`svm_linear_pca150.0_gc_mapq30` — linear kernel, PCA=150, GC correction on, MAPQ≥30, 10×10 CV.

## Command order

Steps 1–3 run **once per cohort**. Steps 4–8 run once per **training cohort**. Steps 9–10 run once per **(training cohort × test cohort)** pair.

```bash
# 1. one-time prerequisites (run from project root)
bash scripts/generalized_processing/static_script.sh

# 2. (internal cohorts only) generate manifest, then manually filter/split and fill stage column
python scripts/generalized_processing/generate_metadata_file.py ly
# → produces data/metadata/internal_metadata_ly.csv
# → manually split into internal_metadata_ly_filtered_train.csv / _test.csv

# 3. per-sample feature computation (N = sample count in manifest)
sbatch --array=1-N%20 scripts/generalized_processing/compute_features_mapq.sbatch \
    30 false data/metadata/Cristiano_metadata.csv data/cristiano_features

# 4. build per-feature matrices
sbatch --array=0-13 --mem=32G scripts/generalized_processing/matrix_build_array.sbatch \
    30 data/cristiano_features data/cristiano_runs/feature_matrices data/metadata/Cristiano_metadata.csv

# 5. per-feature 10×10 CV
sbatch --array=0-6  --mem=16G scripts/generalized_modelling/gen_svm.sbatch cv data/cristiano_runs linear true true 150 10 30
sbatch --array=7-12 --mem=32G scripts/generalized_modelling/gen_svm.sbatch cv data/cristiano_runs linear true true 150 10 30
sbatch --array=13   --mem=64G scripts/generalized_modelling/gen_svm.sbatch cv data/cristiano_runs linear true true 150 10 30

# 6. per-feature training on the full matrix
sbatch --array=0-6  --mem=16G scripts/generalized_modelling/gen_svm.sbatch train data/cristiano_runs linear true true 150 10 30
sbatch --array=7-12 --mem=32G scripts/generalized_modelling/gen_svm.sbatch train data/cristiano_runs linear true true 150 10 30
sbatch --array=13   --mem=64G scripts/generalized_modelling/gen_svm.sbatch train data/cristiano_runs linear true true 150 10 30

# 7. build CV meta-matrix
python scripts/generalized_processing/meta_matrix_build.py \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/cv/ \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv \
    data/metadata/Cristiano_metadata.csv

# 8. ensemble CV + train
# run_tag=nostd, standardize=false — the ensemble operates directly on per-feature probabilities without standardization
sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch cv    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv linear 10 none none none nostd false
sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch train data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv linear 10 none none none nostd false

# 9. per-feature external test (repeat per test cohort)
sbatch --array=0-6  --mem=8G  scripts/generalized_modelling/gen_svm_test.sbatch data/internal_ly_test_runs data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30
sbatch --array=7-12 --mem=16G scripts/generalized_modelling/gen_svm_test.sbatch data/internal_ly_test_runs data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30
sbatch --array=13   --mem=40G scripts/generalized_modelling/gen_svm_test.sbatch data/internal_ly_test_runs data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30

# 10. ensemble external test
# build test meta-matrix from per-feature test probs ($2 must be full output file path, not a directory)
python scripts/generalized_processing/meta_matrix_build.py \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/test/prevelynch_test/ \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/meta_matrix.csv \
    data/metadata/internal_metadata_ly_filtered_test.csv

# ensemble inference (pass $7 output path explicitly — default derivation is wrong for test case)
sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch test \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/meta_matrix.csv \
    linear 10 prevelynch_test \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/models/ensemble.pkl \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_preds.csv \
    nostd false
```
