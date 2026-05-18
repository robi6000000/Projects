# cfDNA_v2 pipeline

Working directory: `/data/projects/liquid_biopsy/Projects/cfDNA_v2/`

## Setup

- Conda env `cfdna` (`source /data/users/kenderes/miniconda3/etc/profile.d/conda.sh && conda activate cfdna`)
- SLURM cluster (gen-compute01-04)
- Tools on PATH: `wget`, `bedtools`, `liftOver`, `samtools`, `python` (pandas, numpy, scikit-learn, statsmodels, pysam)
- One-time prerequisite: `bash scripts/generalized_processing/static_script.sh` - downloads reference files and produces `data/processing/openchrom_with_id.bed` (~560k open chromatin regions) and `data/processing/gc_content_per_region.csv`
- Input: gzipped fragment BED files (chrom, start, end, MAPQ, strand). Run `scripts/snakemake_finaledb/` first if input is paired FASTQ.

## Directory structure

```
cfDNA_v2/
    data/
        metadata/                           # sample metadata CSVs
            Cristiano_metadata.csv
            internal_metadata_ly_filtered_train.csv
            internal_metadata_ly_filtered_test.csv
            internal_metadata_gs.csv
        processing/                         # static reference files
            openchrom_with_id.bed
            openchrom_with_id_centroid.bed
            wps_openchrom_with_id.bed
            gc_content_per_region.csv
        source/                             # downloaded reference files
        hg19/                               # downloaded reference genome
        sample_temp/                        # intermediate BED files 
        cristiano_features/                 # per-sample feature vectors for Cristiano 
        internal_features_ly/               # feature vectors for PreveLynch cohort
        internal_features_gs/               # feature vectors for GenoScan cohort
        cristiano_runs/                     # constructed feature matrices and model runs for Cristiano (same for all _runs folders)
            feature_matrices/             
                matrix_{feature}_mapq{N}.csv
            model_runs/
                svm_{kernel}_{pca}_{gc}_{mapq}/
                    cv/                     # per-feature CV probabilities
                    models/                 # per-feature trained pickle models
                    test/                   # per-feature test predictions
                        {dataset}/
                    meta_matrix/            # merged CV probabilities
                    ensemble_{tag}/         # ensemble model run
                        cv/
                        models/             # ensemble trained model
                        test/
        internal_ly_train_runs/            
        internal_gs_runs/
        crc_runs/
    scripts/
        generalized_processing/
            static_script.sh / .sbatch      # one-time reference file builder
            generate_metadata_file.py       # generate manifest CSV for internal cohorts
            compute_features_mapq.sbatch    
            sample_prep_script.sh           # called by above; creates intersection BEDs
            sample_features_script.py       # called by above; computes 14 features per sample
            matrix_build_array.sbatch      
            matrix_build.py                 # called by above; one feature per invocation
            meta_matrix_build.py            # merge per-feature CV probs into ensemble input
        generalized_modelling/
            gen_svm.sbatch                  # SLURM array job script ,per-feature CV or training
            gen_svm_feature.py
            gen_svm_test.sbatch             # SLURM array: per-feature inference on test cohort
            gen_svm_test.py
            gen_ensemble_svm.sbatch         # SLURM: ensemble CV, training, or test inference
            gen_ensemble_svm.py
            modules/gen_cfdna_model.py      # CFDNAModel class (SVM, GC correction, pca, standardization, CV logic)
```

## Naming conventions

**Metadata CSV** - one row per sample, columns: `sample_id, disease, dataset, material, frag_path, stage, cancer_true`

- generated for internal data from fastq filenames: `python scripts/generalized_processing/generate_metadata_file.py ly`  outputs  `data/metadata/internal_metadata_ly.csv`
- After generation, manually filter/split into train/test and fill the `stage` column

**Fragment files (internal cohorts)**
`{dataset}_{disease}_{...optional}_{material}.GRCh37.frag.bed.gz`
- `dataset` suffix `cc`/`ca` -> `cancer_true=1`; anything else -> 0
- Only `material == pl` (plasma) is kept

**Feature matrix**: `{run_dir}/feature_matrices/matrix_{feature}_mapq{N}.csv`

**Config name** (auto-derived): `svm_{kernel}_pca{N}_gc_cv{R}_mapq{M}_{tag}` - `_cv{R}` omitted when `cv_repeats=10`

## Requirments 

- all commands are ran from the project root 
- Fragment BEDs: autosomes only (chr1-chr22); chr prefix is added by `sample_prep_script.sh` if missing 
- 14 features computed per sample. The features are stored in this order, and can be referenced by slurm array index:
`[length, edm, iedm, eedm, eoedm, cposedm, fsd, pfe, coverage, ends, ocf, ifs, wps, fsr]`
- GC correction is only applied to gc_correct_features = ['pfe', 'coverage', 'ends', 'ocf', 'ifs', 'wps'] - results are still stored under the same feature name regardless of GC correction status, to keep them together with other features in a given RUN
- Memory groups (empirical mostly):
  - light (0-6): 16G for CV; 8G for test
  - med (7-12): 32G for CV/train; 16G for test
  - heavy (13, fsr 1.6m cols): 64G for CV/train; 40G for test (48G currently running, most likely will be ok)
- Ensemble runs are cheap (4G) - operates on a 14-column meta-matrix

- we noticed after switching from pandas to PyArrow for matrix loading, the loading step sped up from hours to a couple of minutes, some modules may not be fully adjusted for this as this was done late in our study

## Best config

`svm_linear_pca150.0_gc_mapq30` - linear kernel, PCA=150, GC correction on, MAPQ>=30, 10x10 CV.

## Command order

Steps 1-3 run once per dataset, steps 4-8 run once per training subset. Steps 9-10 run once for each train/test pair 

```bash
# 1. one-time prerequisites (run from project root)
bash scripts/generalized_processing/static_script.sh

# 2. (internal cohorts only) generates metadata files, then manually filter/split and fill stage column
python scripts/generalized_processing/generate_metadata_file.py ly
# -> produces data/metadata/internal_metadata_ly.csv
# -> manually split into internal_metadata_ly_filtered_train.csv / _test.csv

# 3. per-sample feature computation (N = sample count in manifest)
sbatch --array=1-N%20 scripts/generalized_processing/compute_features_mapq.sbatch \
    30 false data/metadata/Cristiano_metadata.csv data/cristiano_features

# 4. build per-feature matrices
sbatch --array=0-13 --mem=32G scripts/generalized_processing/matrix_build_array.sbatch \
    30 data/cristiano_features data/cristiano_runs/feature_matrices data/metadata/Cristiano_metadata.csv

# 5. per-feature 10x10 CV
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
# run_tag=nostd, standardize=false - the ensemble operates directly on per-feature probabilities without standardization
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

# ensemble inference (pass $7 output path explicitly - default derivation is wrong for test case)
sbatch scripts/generalized_modelling/gen_ensemble_svm.sbatch test \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/meta_matrix.csv \
    linear 10 prevelynch_test \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/models/ensemble.pkl \
    data/cristiano_runs/model_runs/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_preds.csv \
    nostd false
```
