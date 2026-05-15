# cfDNA pipeline — runbook

Concise companion to [protocol_pipeline.md](protocol_pipeline.md). Use this as a quick reference; the longer doc explains the science.

Working dir: `/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/`

## 1. Setup / requirements

- Conda env `cfdna` (`source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh && conda activate cfdna`)
- SLURM cluster (gen-compute01–04)
- Tools available on PATH: `wget`, `bedtools`, `liftOver`, `samtools`, `python` (pandas, numpy, scikit-learn, statsmodels, pysam)
- One-time prerequisite step: `bash scripts/generalized_processing/static_script.sh` — downloads hg19 FASTA, hg38→hg19 chain, Roadmap DNase peaks (E029/E032/E034), TCGA-ATAC pan-cancer peak set; produces `data/processing/openchrom_with_id.bed` (~560k regions), the centroid and WPS-padded variants, and `gc_content_per_region.csv`
- Input: gzipped fragment BED files (chrom, start, end, MAPQ, strand). If only paired FASTQ available, run `scripts/snakemake_finaledb/` first (out of scope for thesis)

## 2. File naming conventions

### Fragment files (internal cohorts)
`{dataset}_{disease}_{...optional}_{material}.GRCh37.frag.bed.gz`
- `dataset` suffix `cc`/`ca` ⇒ `cancer_true=1`; anything else ⇒ 0
- Only `material == pl` (plasma) is kept; everything else excluded
- Example: `lycc_Colorectal_LY00123_pl.GRCh37.frag.bed.gz`

### Manifest CSV (one row per sample)
Columns: `sample_id, disease, dataset, material, frag_path, stage, cancer_true`
- Auto-generated for internal cohorts: `python scripts/generalized_processing/generate_metadata_file.py ly` → `data/manifest/internal_metadata_ly.csv`
- After generation, **manually filter / split** into train/test (e.g. `internal_metadata_ly_filtered_train.csv`, `_test.csv`) and fill the `stage` column where data exists
- External cohorts (e.g. Cristiano) supply the CSV directly — column names must match

### Per-feature outputs
- Per-sample feature CSV: `{features_path}/{feature}_mapq{N}/{sample_id}_{feature}.csv`
- Feature matrix: `{matrix_folder}/by_feature/matrix_{feature}_mapq{N}.csv`
- Model config name (auto-derived): `svm_{kernel}[_pca{N}][_gc][_cv{R}][_mapq{M}][_{tag}]` — `_cv{R}` omitted when `cv_repeats=10`
- Trained pickle: `{config_dir}/models/{config_name}_{feature}.pkl`
- CV probs: `{config_dir}/cv/{feature}_cv_probs.csv`
- Test probs: `{config_dir}/test/{dataset_dir}/{feature}_probs.csv`
- Meta-matrix: `{config_dir}/meta_matrix/meta_matrix.csv`
- Ensemble outputs: `{config_dir}/ensemble_{run_tag}/{cv,models,test}/...`

## 3. Assumptions

- Fragment BEDs are autosomes-relevant only (chr1–chr22 retained; chrom names with or without `chr` prefix both accepted)
- MAPQ filter is applied at fragment-prep time, not at feature time — the same `mapq_filter` value must be threaded through every step or the matrix file names won't resolve
- All 14 features are computed per sample. Feature array index reference (used by every sbatch):
  `0 length, 1 edm, 2 iedm, 3 eedm, 4 eoedm, 5 cposedm, 6 fsd, 7 pfe, 8 coverage, 9 ends, 10 ocf, 11 ifs, 12 wps, 13 fsr`
- GC correction is *only* meaningful for per-region features (`pfe, coverage, ends, ocf, ifs, wps`); the modelling wrapper auto-disables it for the other 8 with a warning
- GC correction is per-sample LOWESS of coverage vs GC, applied independently to train and test (not refit-on-train-and-replay) — see [2026_05_15_protocol.md](2026_05_15_protocol.md)
- Per-region features merge against the full ~560k region list; missing regions are zero-filled
- Trained pickles store: feature name, kernel, gc_correction flag, pca/pca_components, feature_columns (used to align test matrices), scaler (mean/std), PCA model, SVC. No GC state is stored — test-time GC correction is per-sample
- Memory groups (empirical at MAPQ 30, ~459 samples):
  - light (0–6, length / edm-family / fsd): 16G for CV; 8G plausibly enough for train/test (unverified)
  - med (7–12, per-region single-column): 32G for CV/train (peaks ~22.5G — a 24G cap swaps heavily); 16G safe for test
  - heavy (13, fsr ~1.7M cols): 64G for CV/train (safe margin); 40G for test (under monitoring)
- Ensemble runs are cheap (4G) — operates on a 14-column meta-matrix

## 4. Study flow — command order

Steps 1–3 (sample → matrix) run **once per cohort**. Steps 4–8 run once per **training cohort**. Steps 9–10 run once per **(training cohort × test cohort)** pair.

```bash
# 1. one-time prerequisites
bash scripts/generalized_processing/static_script.sh

# 2. (internal cohorts only) generate manifest, then manually filter/split
python scripts/generalized_processing/generate_metadata_file.py ly
# → manually produce internal_metadata_ly_filtered_train.csv / _test.csv

# 3. per-sample feature computation (N = sample count from manifest)
sbatch --array=1-N%20 scripts/generalized_processing/compute_features_mapq.sbatch \
    30 false data/manifest/Cristiano_metadata.csv data/cristiano_features

# 4. build per-feature matrices
sbatch --array=0-6  --mem=8G  scripts/generalized_processing/matrix_build_array.sbatch 30 data/cristiano_features data/matrix/by_feature data/manifest/Cristiano_metadata.csv
sbatch --array=7-12 --mem=32G scripts/generalized_processing/matrix_build_array.sbatch 30 data/cristiano_features data/matrix/by_feature data/manifest/Cristiano_metadata.csv
sbatch --array=13   --mem=64G scripts/generalized_processing/matrix_build_array.sbatch 30 data/cristiano_features data/matrix/by_feature data/manifest/Cristiano_metadata.csv

# 5. per-feature 10×10 CV
sbatch --array=0-6  --mem=16G scripts/generalized_modelling/gen_svm.sbatch cv data/matrix linear true true 150 10 30
sbatch --array=7-12 --mem=32G scripts/generalized_modelling/gen_svm.sbatch cv data/matrix linear true true 150 10 30
sbatch --array=13   --mem=64G scripts/generalized_modelling/gen_svm.sbatch cv data/matrix linear true true 150 10 30

# 6. per-feature training on the full matrix
sbatch --array=0-6  --mem=16G scripts/generalized_modelling/gen_svm.sbatch train data/matrix linear true true 150 10 30
sbatch --array=7-12 --mem=32G scripts/generalized_modelling/gen_svm.sbatch train data/matrix linear true true 150 10 30
sbatch --array=13   --mem=64G scripts/generalized_modelling/gen_svm.sbatch train data/matrix linear true true 150 10 30

# 7. build CV meta-matrix
python scripts/generalized_processing/meta_matrix_build.py \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/cv/ \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/ \
    data/manifest/Cristiano_metadata.csv

# 8. ensemble CV + train
sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch cv    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv linear 10 none none none nostd false
sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch train data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv linear 10 none none none nostd false

# 9. per-feature external test (repeat per test cohort)
sbatch --array=0-6  --mem=8G  scripts/generalized_modelling/gen_svm_test.sbatch data/matrix_internal_ly_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30
sbatch --array=7-12 --mem=16G scripts/generalized_modelling/gen_svm_test.sbatch data/matrix_internal_ly_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30
sbatch --array=13   --mem=40G scripts/generalized_modelling/gen_svm_test.sbatch data/matrix_internal_ly_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30

# 10. ensemble external test (build test meta-matrix from per-feature test probs, then run ensemble inference)
python scripts/generalized_processing/meta_matrix_build.py \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/test/prevelynch_test/ \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/ \
    data/manifest/internal_metadata_ly_filtered_test.csv

sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch test \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/meta_matrix.csv \
    linear 10 prevelynch_test \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/models/ensemble.pkl \
    none nostd false
```

## 5. Best config (used throughout thesis)

`svm_linear_pca150.0_gc_mapq30` — linear kernel, PCA=150, GC correction on, MAPQ≥30, 10×10 CV, ensemble run-tag `nostd` (no standardisation at the ensemble layer).
